"""
Reading Octet (ForteBio / Sartorius BLI) exports and turning each sensor trace into a
single binding readout.

Export layout (Octet Data Analysis "export to Excel/CSV"): one column pair per well,
time (X) then response (Y). Somewhere above the data, the X column carries
``Sample Loc: <well>`` and the Y column ``Sample ID: <text>``. Numeric rows below are
the trace. Nothing else about the sheet is assumed.

The Sample ID text carries the experiment, for example::

    hEphB2-His (80 sec) + 0.1 uM A20-I.2xHCl in PBSTD
    \\_______/  \\_____/   \\____/ \\__________/    \\_____/
     analyte  association  dose   compound     vehicle

A well without ``+ <dose> <unit> <compound>`` is a no-inhibitor (100 %) reference.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re

import numpy as np
import pandas as pd

UNITS_UM = {"pm": 1e-6, "nm": 1e-3, "um": 1.0, "µm": 1.0, "μm": 1.0, "mm": 1e3, "m": 1e6}
_TREAT = re.compile(r"\+\s*(\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(pM|nM|uM|µM|μM|mM|M)\b\s*(.*)$")
_ASSOC = re.compile(r"\(\s*(\d+(?:\.\d+)?)\s*s(?:ec(?:onds?)?)?\s*\)", re.I)
_WELL = re.compile(r"^\s*([A-Pa-p])\s*0*(\d{1,2})\s*$")
DEFAULT_BLANK = r"(?i)^\s*(association|buffer|blank|baseline)\b"


@dataclass
class SampleInfo:
    text: str
    analyte: str
    association_s: float | None
    concentration_um: float | None
    compound: str
    vehicle: str

    @property
    def has_inhibitor(self) -> bool:
        return self.concentration_um is not None


def parse_sample_id(text) -> SampleInfo:
    """Split a Sample ID into analyte, association time, inhibitor dose and compound."""
    t = re.sub(r"^\s*sample\s*id\s*:?\s*", "", str(text or ""), flags=re.I).strip()
    m = _ASSOC.search(t)
    assoc = float(m.group(1)) if m else None
    tm = _TREAT.search(t)
    conc, compound, vehicle = None, "", ""
    head = t
    if tm:
        conc = float(tm.group(1)) * UNITS_UM[tm.group(2).lower()]
        rest = tm.group(3).strip()
        vm = re.search(r"\s+in\s+(.+)$", rest)
        if vm:
            vehicle, rest = vm.group(1).strip(), rest[:vm.start()].strip()
        compound = rest
        head = t[:tm.start()]
    else:
        vm = re.search(r"\s+in\s+(.+)$", t)
        if vm:
            vehicle, head = vm.group(1).strip(), t[:vm.start()]
    analyte = _ASSOC.sub("", head).strip(" +,;")
    return SampleInfo(t, analyte, assoc, conc, compound, vehicle)


def plate_column(loc: str):
    """Plate column number of a well location ("C5" -> 5), or None."""
    m = _WELL.match(str(loc or ""))
    return int(m.group(2)) if m else None


@dataclass
class Well:
    index: int
    loc: str
    sample: SampleInfo
    time: np.ndarray = field(repr=False)
    response: np.ndarray = field(repr=False)

    @property
    def column(self):
        return plate_column(self.loc)


def _read_table(path, sheet=0) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in (".csv", ".tsv", ".txt"):
        sep = "\t" if path.suffix.lower() == ".tsv" else ","
        return pd.read_csv(path, header=None, sep=sep, dtype=object, encoding="utf-8-sig",
                           engine="python", on_bad_lines="skip")
    return pd.read_excel(path, sheet_name=sheet, header=None, dtype=object)


def read_octet(path, sheet=0) -> list:
    """All wells (time/response traces with parsed Sample IDs) of an Octet export."""
    raw = _read_table(path, sheet)
    cells = raw.astype(str)
    is_id = cells.apply(lambda col: col.str.contains(r"(?i)^\s*sample\s*id", regex=True))
    if not is_id.values.any():
        raise ValueError(f"{Path(path).name}: no 'Sample ID' header found - is this an Octet export?")
    wells = []
    for j in np.nonzero(is_id.values.any(axis=0))[0]:
        if j == 0:
            continue
        r = int(np.nonzero(is_id.iloc[:, j].values)[0][0])
        loc = str(raw.iat[r, j - 1]) if pd.notna(raw.iat[r, j - 1]) else ""
        loc = re.sub(r"(?i)^\s*sample\s*loc\w*\s*:?\s*", "", loc).strip()
        t = pd.to_numeric(raw.iloc[r + 1:, j - 1], errors="coerce").to_numpy(float)
        y = pd.to_numeric(raw.iloc[r + 1:, j], errors="coerce").to_numpy(float)
        ok = np.isfinite(t) & np.isfinite(y)
        if ok.sum() < 2:
            continue
        wells.append(Well(len(wells), loc, parse_sample_id(raw.iat[r, j]), t[ok], y[ok]))
    if not wells:
        raise ValueError(f"{Path(path).name}: Sample ID headers found but no numeric traces")
    return wells


# ---------------------------------------------------------------- readouts

def value_at(t, y, time):
    """Response at `time` (linear interpolation; nan outside the trace)."""
    if time < t[0] or time > t[-1]:
        return np.nan
    return float(np.interp(time, t, y))


def auc(t, y, t0=None, t1=None):
    """Trapezoidal area under the trace between t0 and t1 (trace ends by default)."""
    t0 = t[0] if t0 is None else max(t0, t[0])
    t1 = t[-1] if t1 is None else min(t1, t[-1])
    if t1 <= t0:
        return np.nan
    inside = (t > t0) & (t < t1)
    tt = np.concatenate([[t0], t[inside], [t1]])
    yy = np.concatenate([[np.interp(t0, t, y)], y[inside], [np.interp(t1, t, y)]])
    return float(np.sum((yy[1:] + yy[:-1]) * np.diff(tt)) / 2.0)


@dataclass(frozen=True)
class Readout:
    """How a trace becomes one number. kind: "point" (value at `start`) or "auc"
    (area from `start` to `end`; None = trace start / end). Times are seconds from the
    first time point of each trace."""
    name: str
    label: str
    kind: str
    start: float | None = None
    end: float | None = None

    def __call__(self, t, y):
        t = np.asarray(t, float) - float(t[0])
        y = np.asarray(y, float)
        if self.kind == "point":
            return value_at(t, y, self.start)
        return auc(t, y, self.start, self.end)


def default_readouts(association_s=None, assoc_time=None, dissoc_time=None) -> list:
    """The four readouts of the original analysis: response 5 s before the end of
    association, response 50 s into dissociation, AUC of the whole trace and AUC of the
    association phase (75 s, 130 s, 0-end and 0-80 s for an 80 s association)."""
    out = []
    if association_s is None and (assoc_time is None or dissoc_time is None):
        return [Readout("auc_full", "Full AUC", "auc")]
    ta = assoc_time if assoc_time is not None else association_s - 5.0
    td = dissoc_time if dissoc_time is not None else association_s + 50.0
    out.append(Readout("association", f"Association at {ta:g} s", "point", ta))
    out.append(Readout("dissociation", f"Dissociation at {td:g} s", "point", td))
    out.append(Readout("auc_full", "Full AUC", "auc"))
    if association_s is not None:
        out.append(Readout("auc_association", f"Association AUC (0-{association_s:g} s)", "auc",
                           0.0, association_s))
    return out
