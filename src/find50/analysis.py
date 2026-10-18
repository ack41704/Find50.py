"""
Plate analysis: wells -> roles -> readouts -> normalised dose-response -> IC50 table.

Roles of the wells:

* **reference**  no inhibitor in the Sample ID: 100 % activity
* **control**    the positive (maximal-inhibition) control given with ``control=``,
                 e.g. ``"A20-I@2"`` = compound A20-I at 2 uM: defines 0 % for the
                 *relative* normalisation
* **dose**       every other well with an inhibitor dose
* **blank**      buffer-only wells (``blank_regex``), ignored
* **excluded**   wells listed in ``exclude`` (well location, e.g. "C5"), ignored

Each readout r of a dose well is normalised two ways:

    absolute   % activity = 100 * r / r_ref
    relative   % activity = 100 * (r - r_ctrl) / (r_ref - r_ctrl)

where r_ref and r_ctrl are the means of the reference and control wells of the same
analyte. With ``reference_scope="column"`` (default) the references are taken from the
same plate column as the dose well when that column has any (one compound per column is
the usual Octet plate layout), otherwise from the whole file.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
import re

import numpy as np
import pandas as pd

from .curve import ic50
from .octet import DEFAULT_BLANK, default_readouts, read_octet, UNITS_UM


@dataclass(frozen=True)
class Control:
    """Positive-control definition: compound name (case-insensitive prefix, so "A20-I"
    matches "A20-I.2xHCl") and optionally its concentration in uM."""
    compound: str
    concentration_um: float | None = None

    @classmethod
    def parse(cls, text):
        if text is None or isinstance(text, Control):
            return text
        name, _, conc = str(text).partition("@")
        if not conc:
            return cls(name.strip())
        m = re.fullmatch(r"\s*(\d*\.?\d+)\s*(pM|nM|uM|µM|μM|mM|M)?\s*", conc)
        if not m:
            raise ValueError(f"control must look like NAME or NAME@CONC_uM, not {text!r}")
        return cls(name.strip(), float(m.group(1)) * UNITS_UM[(m.group(2) or "uM").lower()])

    def matches(self, sample) -> bool:
        if not sample.has_inhibitor or not sample.compound.casefold().startswith(self.compound.casefold()):
            return False
        return self.concentration_um is None or np.isclose(sample.concentration_um, self.concentration_um,
                                                           rtol=1e-6)


@dataclass
class Analysis:
    wells: list
    readouts: list
    well_table: pd.DataFrame
    points: pd.DataFrame
    results: pd.DataFrame
    curves: dict = field(default_factory=dict)      # curve id -> list of dose-well indices
    fits: dict = field(default_factory=dict, repr=False)   # (curve, readout, norm) -> IC50Result
    settings: dict = field(default_factory=dict)

    def summary(self) -> pd.DataFrame:
        """Wide table like the original output: one row per readout x normalisation,
        one column per curve, cells = IC50 labels (uM)."""
        if self.results.empty:
            return pd.DataFrame()
        r = self.results.copy()
        r["row"] = np.where(r["normalization"] == "relative", "Relative " + r["readout_label"].str[0].str.lower()
                            + r["readout_label"].str[1:], r["readout_label"])
        order = list(dict.fromkeys(r["row"]))
        wide = r.pivot(index="row", columns="curve", values="label").reindex(order)
        wide = wide[list(dict.fromkeys(r["curve"]))]
        wide.index.name = "IC50 (uM)"
        return wide


def _role(w, control, blank, exclude):
    if w.loc.upper() in exclude or str(w.index) in exclude:
        return "excluded"
    if blank and re.search(blank, w.sample.text):
        return "blank"
    if control is not None and control.matches(w.sample):
        return "control"
    return "dose" if w.sample.has_inhibitor else "reference"


def _curve_ids(dose_wells):
    """Group dose wells into curves by (analyte, compound). Replicate series of the same
    compound (repeated concentrations) are split by plate column when every column holds
    a series of its own. Groups with a single concentration are not curves."""
    groups = defaultdict(list)
    for w in dose_wells:
        groups[(w.sample.analyte, w.sample.compound)].append(w)
    out = {}
    for (analyte, compound), ws in groups.items():
        name = f"{compound} | {analyte}" if analyte else compound
        if len({w.sample.concentration_um for w in ws}) < 2:
            continue
        concs = Counter(w.sample.concentration_um for w in ws)
        by_col = defaultdict(list)
        for w in ws:
            by_col[w.column].append(w)
        if (max(concs.values()) > 1 and len(by_col) > 1 and None not in by_col
                and all(len({w.sample.concentration_um for w in v}) >= 2 for v in by_col.values())):
            for k, col in enumerate(sorted(by_col), 1):
                out[f"{name} (rep {k}, col {col})"] = by_col[col]
        else:
            out[name] = ws
    return out


def analyze(wells, control=None, readouts=None, association_s=None, exclude=(),
            reference_scope: str = "column", r2_min: float = 0.8,
            blank_regex: str | None = DEFAULT_BLANK) -> Analysis:
    """IC50s for every inhibitor series in a list of wells (see module docstring)."""
    control = Control.parse(control)
    excl = {str(e).strip().upper() for e in exclude}
    if association_s is None:
        found = [w.sample.association_s for w in wells if w.sample.association_s]
        association_s = Counter(found).most_common(1)[0][0] if found else None
    if readouts is None:
        readouts = default_readouts(association_s)
    roles = {w.index: _role(w, control, blank_regex, excl) for w in wells}
    values = {w.index: {r.name: r(w.time, w.response) for r in readouts} for w in wells}

    rows = []
    for w in wells:
        s = w.sample
        rows.append({"well": w.index, "loc": w.loc, "plate_column": w.column, "role": roles[w.index],
                     "sample_id": s.text, "analyte": s.analyte, "compound": s.compound,
                     "concentration_uM": s.concentration_um, "association_s": s.association_s,
                     "vehicle": s.vehicle, "n_timepoints": len(w.time),
                     **{f"readout:{k}": v for k, v in values[w.index].items()}})
    well_table = pd.DataFrame(rows)

    def refs(w, role):
        """Wells of `role` for dose well `w` (same analyte; same plate column if any)."""
        pool = [v for v in wells if roles[v.index] == role and v.sample.analyte == w.sample.analyte]
        if reference_scope == "column" and w.column is not None:
            same = [v for v in pool if v.column == w.column]
            if same:
                return same, "column"
        return pool, ("file" if pool else "")

    def mean_of(ws, readout):
        return float(np.nanmean([values[v.index][readout] for v in ws])) if ws else np.nan

    curves = _curve_ids([w for w in wells if roles[w.index] == "dose"])
    member = {w.index: cid for cid, ws in curves.items() for w in ws}
    well_table.insert(4, "curve", [member.get(i, "") for i in well_table["well"]])
    points, results, fits = [], [], {}
    for cid, ws in curves.items():
        for r in readouts:
            for norm in ("absolute", "relative"):
                ctrl = {w.index: refs(w, "control") for w in ws}
                if norm == "relative" and not all(c[0] for c in ctrl.values()):
                    continue                      # no positive control for this curve
                conc, act, flags = [], [], set()
                for w in ws:
                    val = values[w.index][r.name]
                    ref_w, scope = refs(w, "reference")
                    top = mean_of(ref_w, r.name)
                    if not np.isfinite(top):
                        top = float(np.nanmax([values[v.index][r.name] for v in ws]))
                        scope = "max dose well"
                        flags.add("no reference wells: 100 % = highest dose-well response")
                    bottom = mean_of(ctrl[w.index][0], r.name) if norm == "relative" else 0.0
                    pct = 100.0 * (val - bottom) / (top - bottom) if top != bottom else np.nan
                    conc.append(w.sample.concentration_um)
                    act.append(pct)
                    points.append({"curve": cid, "readout": r.name, "readout_label": r.label,
                                   "normalization": norm, "well": w.index, "loc": w.loc,
                                   "concentration_uM": w.sample.concentration_um, "value": val,
                                   "reference_value": top,
                                   "control_value": bottom if norm == "relative" else np.nan,
                                   "reference_scope": scope,
                                   "control_scope": ctrl[w.index][1] if norm == "relative" else "",
                                   "pct_activity": pct, "pct_inhibition": 100.0 - pct})
                res = ic50(conc, act, r2_min=r2_min)
                res.flags = sorted(flags) + res.flags
                fits[(cid, r.name, norm)] = res
                s = ws[0].sample
                results.append({"curve": cid, "analyte": s.analyte, "compound": s.compound,
                                "readout": r.name, "readout_label": r.label, "normalization": norm,
                                "cmin_uM": float(np.nanmin(conc)), "cmax_uM": float(np.nanmax(conc)),
                                **res.as_dict()})
    settings = {"control": f"{control.compound} @ {control.concentration_um} uM" if control else "none",
                "association_s": association_s, "reference_scope": reference_scope,
                "r2_min": r2_min, "blank_regex": blank_regex, "excluded": sorted(excl),
                "readouts": "; ".join(f"{r.name} = {r.label}" for r in readouts)}
    return Analysis(wells, readouts, well_table, pd.DataFrame(points), pd.DataFrame(results),
                    {k: [w.index for w in v] for k, v in curves.items()}, fits, settings)


def analyze_file(path, sheet=0, **kwargs) -> Analysis:
    a = analyze(read_octet(path, sheet), **kwargs)
    a.settings["file"] = str(path)
    return a
