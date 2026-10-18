import numpy as np
import pandas as pd
import pytest

from find50 import Control, analyze, analyze_file, parse_sample_id, read_octet
from find50.octet import Readout, auc, default_readouts, plate_column, value_at
from find50.synthetic import synthetic_plate, write_demo_file, write_octet_xlsx


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    path, truth = write_demo_file(tmp_path_factory.mktemp("octet") / "plate.xlsx")
    return path, truth


@pytest.mark.parametrize("text, analyte, assoc, conc, compound, vehicle", [
    ("Sample ID: hEphB2-His (80 sec) + 0.1 uM EDTA", "hEphB2-His", 80, 0.1, "EDTA", ""),
    (" Sample ID: 25nM hEphB2-His (800 sec) + 0.05 uM BQPB4.3xHCl in PBSTD",
     "25nM hEphB2-His", 800, 0.05, "BQPB4.3xHCl", "PBSTD"),
    ("25nM hEphB2-His (800 sec) in PBSTD", "25nM hEphB2-His", 800, None, "", "PBSTD"),
    ("hTarget + 500 nM Cmpd-7", "hTarget", None, 0.5, "Cmpd-7", ""),
    ("hTarget (60 s) + 2 mM Cmpd-8", "hTarget", 60, 2000.0, "Cmpd-8", ""),
])
def test_parse_sample_id(text, analyte, assoc, conc, compound, vehicle):
    s = parse_sample_id(text)
    assert (s.analyte, s.association_s, s.compound, s.vehicle) == (analyte, assoc, compound, vehicle)
    assert s.concentration_um == (pytest.approx(conc) if conc is not None else None)


def test_plate_column():
    assert plate_column("C5") == 5 and plate_column("h12") == 12 and plate_column("x") is None


def test_readouts_on_a_linear_trace():
    t = np.arange(10.0, 210.25, 0.5)             # 10-210 s: traces need not start at 0
    y = 0.01 * (t - 10.0)
    assert value_at(t - 10, y, 75.0) == pytest.approx(0.75)
    assert Readout("p", "p", "point", 75.0)(t, y) == pytest.approx(0.75)
    assert Readout("a", "a", "auc")(t, y) == pytest.approx(0.01 * 200 ** 2 / 2)
    assert Readout("w", "w", "auc", 0.0, 80.0)(t, y) == pytest.approx(0.01 * 80 ** 2 / 2)
    assert auc(np.array([0, 1.0]), np.array([1, 1.0]), 0.25, 0.75) == pytest.approx(0.5)
    assert np.isnan(value_at(t, y, 1000.0))


def test_default_readouts_follow_the_original_times():
    r = {x.name: x for x in default_readouts(80.0)}
    assert r["association"].start == 75.0 and r["dissociation"].start == 130.0
    assert (r["auc_association"].start, r["auc_association"].end) == (0.0, 80.0)
    assert [x.name for x in default_readouts(None)] == ["auc_full"]


def test_control_parsing():
    assert Control.parse("A20-I@2") == Control("A20-I", 2.0)
    assert Control.parse("A20-I@500nM") == Control("A20-I", 0.5)
    c = Control.parse("A20-I@2")
    assert c.matches(parse_sample_id("x + 2 uM A20-I.2xHCl"))
    assert not c.matches(parse_sample_id("x + 0.2 uM A20-I.2xHCl"))     # the test series
    with pytest.raises(ValueError):
        Control.parse("A20-I@two")


def test_read_octet_layout(demo):
    path, truth = demo
    wells = read_octet(path)
    assert len(wells) == 8 * len(truth)
    assert wells[0].loc == "A1" and wells[0].column == 1
    assert wells[0].time[0] == 0 and len(wells[0].time) == 401


def test_csv_export_is_read_the_same(demo, tmp_path):
    path, _ = demo
    csv = tmp_path / "plate.csv"
    pd.read_excel(path, header=None).to_csv(csv, header=False, index=False)
    a, b = read_octet(path), read_octet(csv)
    assert [w.sample.text for w in a] == [w.sample.text for w in b]
    assert np.allclose(a[5].response, b[5].response)


def test_recovers_known_ic50s(demo):
    path, truth = demo
    an = analyze_file(path, control="CTRL-1@2")
    r = an.results.set_index(["compound", "readout", "normalization"])
    roles = an.well_table["role"].value_counts()
    assert (roles["reference"], roles["control"], roles["dose"]) == (4, 4, 24)
    for name in ("Cmpd-1", "Cmpd-2"):                      # inside the tested 0.1-3.2 uM range
        for readout in ("association", "dissociation", "auc_full", "auc_association"):
            row = r.loc[(name, readout, "relative")]
            assert row["method"] == "4PL"
            assert row["ic50"] == pytest.approx(truth[name], rel=0.1)
            # absolute normalisation keeps the inhibitor-independent signal -> higher IC50
            assert r.loc[(name, readout, "absolute")]["ic50"] > row["ic50"] * 0.97
    assert set(r.xs("Cmpd-3")["label"]) == {"> 3.2"}
    assert set(r.xs("Cmpd-4")["label"]) == {"< 0.1"}
    wide = an.summary()
    assert wide.shape == (8, 4) and wide.index[1] == "Relative association at 75 s"


def test_without_control_only_absolute(demo):
    path, _ = demo
    an = analyze_file(path)
    assert set(an.results["normalization"]) == {"absolute"}
    assert set(an.well_table["role"]) == {"reference", "dose"}
    # the control is now an ordinary dose well, but one concentration is not a curve
    assert len(an.curves) == 4 and not any(c.startswith("CTRL-1") for c in an.curves)
    ctrl = an.well_table[an.well_table["compound"] == "CTRL-1"]
    assert (ctrl["curve"] == "").all()


def test_exclude_wells(demo):
    path, _ = demo
    an = analyze_file(path, control="CTRL-1@2", exclude=["c1", "D1"])
    n = an.results.set_index(["compound", "readout", "normalization"]).loc[("Cmpd-1", "auc_full", "absolute")]
    assert n["n_points"] == 4
    assert (an.well_table.set_index("loc").loc[["C1", "D1"], "role"] == "excluded").all()


def test_reference_scope_column_vs_file():
    wells, _ = synthetic_plate({"Cmpd-1": 0.3, "Cmpd-2": 1.0}, seed=2)
    loc, sid, t, y = wells[8]                                 # make column 2's reference brighter
    wells[8] = (loc, sid, t, y * 1.5)
    from find50.octet import Well
    ws = [Well(i, w[0], parse_sample_id(w[1]), w[2], w[3]) for i, w in enumerate(wells)]
    col = analyze(ws, control="CTRL-1@2").points
    fil = analyze(ws, control="CTRL-1@2", reference_scope="file").points
    pick = lambda p: p[(p.curve.str.startswith("Cmpd-2")) & (p.readout == "auc_full")  # noqa: E731
                       & (p.normalization == "absolute")]
    assert set(pick(col)["reference_scope"]) == {"column"}
    assert (pick(col)["pct_activity"].to_numpy() < pick(fil)["pct_activity"].to_numpy()).all()


def test_replicate_series_are_split_by_plate_column(tmp_path):
    wells, _ = synthetic_plate({"Cmpd-1": 0.3}, seed=4)
    wells2, _ = synthetic_plate({"Cmpd-1": 0.3}, seed=5)
    wells += [(loc.replace("1", "2"), sid, t, y) for loc, sid, t, y in wells2]
    p = write_octet_xlsx(wells, tmp_path / "reps.xlsx")
    an = analyze_file(p, control="CTRL-1@2")
    assert len(an.curves) == 2 and all("rep" in c for c in an.curves)


def test_not_an_octet_file(tmp_path):
    p = tmp_path / "x.csv"
    p.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError, match="Sample ID"):
        read_octet(p)
