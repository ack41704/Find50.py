import openpyxl
import pandas as pd

from find50 import cli
from find50.synthetic import write_demo_file


def test_demo(tmp_path, capsys):
    assert cli.main(["demo", "--out", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "Relative full AUC" in out and "true IC50" in out
    assert (tmp_path / "synthetic_curves.pdf").stat().st_size > 10_000


def test_octet_command_writes_workbook(tmp_path):
    path, truth = write_demo_file(tmp_path / "plate.xlsx")
    assert cli.main(["octet", str(path), "--control", "CTRL-1@2", "--no-plots",
                     "--out", str(tmp_path / "res")]) == 0
    wb = openpyxl.load_workbook(tmp_path / "res" / "plate_results.xlsx")
    assert wb.sheetnames == ["IC50 summary", "IC50 details", "Dose-response", "Wells", "Settings"]
    details = pd.read_csv(tmp_path / "res" / "plate_ic50.csv")
    assert len(details) == 4 * 4 * 2                    # compounds x readouts x normalisations
    assert not (tmp_path / "res" / "plate_curves.pdf").exists()


def test_custom_readout_times(tmp_path):
    path, _ = write_demo_file(tmp_path / "plate.xlsx")
    cli.main(["octet", str(path), "--assoc-time", "60", "--dissoc-time", "150", "--association-s", "80",
              "--no-plots", "--out", str(tmp_path / "r")])
    d = pd.read_csv(tmp_path / "r" / "plate_ic50.csv")
    assert {"Association at 60 s", "Dissociation at 150 s"} <= set(d["readout_label"])


def test_curve_command(capsys, tmp_path):
    assert cli.main(["curve", "--conc", "0.1", "0.3", "1", "3", "10", "30",
                     "--response", "98", "90", "70", "40", "15", "5", "--top", "100",
                     "--plot", str(tmp_path / "c.png")]) == 0
    out = capsys.readouterr().out
    assert out.startswith("IC50 = ") and "(4PL)" in out
    assert (tmp_path / "c.png").exists()
