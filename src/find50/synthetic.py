"""
Synthetic Octet export with known IC50s (for tests, the demo and the README figures).

Each compound occupies one plate column (the layout used by the original notebook):
row A = no-inhibitor reference, row B = positive control, rows C-H = six doses.
Binding of the analyte to the sensor is inhibited as

    f(c) = 1 / (1 + (c / IC50) ** hill)                      fraction still binding

    association (0 <= t <= ta):   R(t) = Rmax f (1 - exp(-kobs t)) + N(t)
    dissociation (t > ta):        R(t) = R(ta) specific * exp(-koff (t - ta)) + N(t)

where N(t) is a fast, inhibitor-independent non-specific signal. After subtracting
the positive control (relative normalisation) every readout gives back the true IC50.
The absolute normalisation still contains N (unaffected by the inhibitor), so absolute
IC50s come out higher than the true value.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_DOSES = (0.1, 0.2, 0.4, 0.8, 1.6, 3.2)
DEFAULT_COMPOUNDS = {"Cmpd-1": 0.30, "Cmpd-2": 1.00, "Cmpd-3": 12.0, "Cmpd-4": 0.02}


def trace(conc_um, ic50_um, hill=1.0, association_s=80.0, dissociation_s=120.0, dt=0.5,
          rmax=0.30, kobs=0.06, koff=0.012, nonspecific=0.04, noise=0.002, rng=None):
    rng = rng or np.random.default_rng(0)
    t = np.arange(0.0, association_s + dissociation_s + dt / 2, dt)
    f = 1.0 / (1.0 + (conc_um / ic50_um) ** hill) if conc_um > 0 else 1.0
    spec_a = rmax * f * (1 - np.exp(-kobs * np.minimum(t, association_s)))
    spec = np.where(t <= association_s, spec_a,
                    spec_a * np.exp(-koff * np.clip(t - association_s, 0, None)))
    ns_a = nonspecific * (1 - np.exp(-0.25 * np.minimum(t, association_s)))
    ns = np.where(t <= association_s, ns_a, ns_a * np.exp(-0.08 * np.clip(t - association_s, 0, None)))
    return t, spec + ns + rng.normal(0, noise, t.size)


def synthetic_plate(compounds=None, doses=DEFAULT_DOSES, analyte="hTarget-His",
                    association_s=80.0, control=("CTRL-1", 2.0, 0.001), seed=0, **kw):
    """List of (sample_loc, sample_id, t, response) and the truth {compound: IC50}."""
    compounds = dict(compounds or DEFAULT_COMPOUNDS)
    rng = np.random.default_rng(seed)
    rows = "ABCDEFGH"
    wells = []
    for col, (name, true_ic50) in enumerate(compounds.items(), start=1):
        layout = [(f"{analyte} ({association_s:g} sec) in PBST", 0.0, 1.0),
                  (f"{analyte} ({association_s:g} sec) + {control[1]:g} uM {control[0]} in PBST",
                   control[1], control[2])]
        layout += [(f"{analyte} ({association_s:g} sec) + {d:g} uM {name} in PBST", d, true_ic50)
                   for d in doses]
        for r, (sid, c, k) in zip(rows, layout):
            t, y = trace(c, k, association_s=association_s, rng=rng, **kw)
            wells.append((f"{r}{col}", sid, t, y))
    return wells, compounds


def write_octet_xlsx(wells, path) -> Path:
    """Write wells in the Octet export layout: per well an X (time) and Y (response)
    column; header rows 1-5 (sensor, sample location / Sample ID, concentration,
    X/Y labels), data from row 6."""
    n = max(len(t) for _, _, t, _ in wells)
    cols = {}
    for i, (loc, sid, t, y) in enumerate(wells):
        conc = sid.split("+")[1].split("uM")[0].strip() if "+" in sid else "0"
        head_x = [f"Sensor Loc: {loc}", "Sensor Type: SA", f"Sample Loc: {loc}", "Sample Conc.:", "X Values"]
        head_y = ["", "", f"Sample ID: {sid}", conc, "Y Values"]
        pad = [None] * (n - len(t))
        cols[2 * i] = head_x + list(t) + pad
        cols[2 * i + 1] = head_y + list(y) + pad
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(cols).to_excel(path, header=False, index=False)
    return path


def write_demo_file(path, seed=0, **kw):
    wells, truth = synthetic_plate(seed=seed, **kw)
    return write_octet_xlsx(wells, path), truth
