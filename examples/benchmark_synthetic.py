"""
Synthetic benchmark: does Find50 recover a known IC50 from simulated Octet traces?

Each simulated compound has a random true IC50 inside the tested range (0.1-3.2 uM,
six 2-fold doses) and a random Hill slope; traces include an inhibitor-independent
non-specific signal and measurement noise (see find50.synthetic). Writes
benchmark_synthetic.csv and prints the accuracy per readout and normalisation.
Synthetic data checks the implementation; it does not replace orthogonal assays.

    python examples/benchmark_synthetic.py [--n 60] [--noise 0.002] [--out docs/benchmark]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from find50 import analyze, parse_sample_id
from find50.octet import Well
from find50.synthetic import synthetic_plate


def run(n=60, noise=0.002, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        true = float(10 ** rng.uniform(np.log10(0.15), np.log10(2.5)))
        hill = float(rng.uniform(0.7, 1.8))
        wells, _ = synthetic_plate({"X": true}, seed=seed * 10_000 + i, hill=hill, noise=noise)
        ws = [Well(k, w[0], parse_sample_id(w[1]), w[2], w[3]) for k, w in enumerate(wells)]
        # the control well uses the same Hill slope; its IC50 (1 nM) makes it ~0 % activity
        an = analyze(ws, control="CTRL-1@2")
        for r in an.results.itertuples():
            rows.append({"compound": i, "true_ic50_uM": true, "hill_true": hill, "readout": r.readout,
                         "normalization": r.normalization, "ic50_uM": r.ic50, "method": r.method,
                         "qualifier": r.qualifier, "r2": r.r2})
    df = pd.DataFrame(rows)
    df["fold_error"] = 10 ** np.abs(np.log10(df.ic50_uM / df.true_ic50_uM))
    return df


def summary(df):
    g = df.groupby(["readout", "normalization"], sort=False)
    return pd.DataFrame({
        "curves": g.size(),
        "4PL reported (%)": (100 * g["method"].apply(lambda m: (m == "4PL").mean())).round(1),
        "median fold error": g["fold_error"].median().round(3),
        "max fold error": g["fold_error"].max().round(3),
        "within 1.25-fold (%)": (100 * g["fold_error"].apply(lambda f: (f <= 1.25).mean())).round(1),
    })


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, default=60)
    ap.add_argument("--noise", type=float, default=0.002)
    ap.add_argument("--out", type=Path, default=Path("docs/benchmark"))
    a = ap.parse_args()
    df = run(a.n, a.noise)
    a.out.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out / "benchmark_synthetic.csv", index=False)
    s = summary(df)
    s.to_csv(a.out / "benchmark_summary.csv")
    print(s.to_string())
