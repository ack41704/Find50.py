"""
Regenerate the README figures from synthetic data.

    python docs/make_figures.py

Writes docs/figures/method_overview.png and docs/benchmark/*.csv.
"""
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from find50 import analyze, parse_sample_id  # noqa: E402
from find50.octet import Well  # noqa: E402
from find50.plotting import ABS_C, REL_C, plot_dose_response  # noqa: E402
from find50.synthetic import synthetic_plate  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "examples"))
from benchmark_synthetic import run, summary  # noqa: E402

plt.rcParams.update({"font.family": "sans-serif", "font.size": 8.5, "axes.titlesize": 9.5,
                     "axes.titleweight": "bold", "axes.spines.top": False,
                     "axes.spines.right": False, "savefig.dpi": 200})


def traces_panel(ax, an, curve):
    by = {w.index: w for w in an.wells}
    doses = sorted((by[i] for i in an.curves[curve]), key=lambda w: w.sample.concentration_um)
    col = doses[0].column
    roles = dict(zip(an.well_table["well"], an.well_table["role"]))
    ref = next(w for w in an.wells if roles[w.index] == "reference" and w.column == col)
    ctl = next(w for w in an.wells if roles[w.index] == "control" and w.column == col)
    cmap = plt.get_cmap("viridis")
    ax.plot(ref.time, ref.response, color="black", lw=1.3, label="no inhibitor")
    for k, w in enumerate(doses):
        ax.plot(w.time, w.response, color=cmap(0.9 * k / (len(doses) - 1)), lw=1,
                label=f"{w.sample.concentration_um:g} µM")
    ax.plot(ctl.time, ctl.response, color="0.5", lw=1.2, ls="--", label="positive control")
    a_s = an.settings["association_s"]
    ax.axvspan(0, a_s, color="0.93", zorder=0)
    ax.text(a_s / 2, 0.355, "association", ha="center", fontsize=7, color="0.35")
    ax.text(a_s + 75, 0.355, "dissociation", ha="center", fontsize=7, color="0.35")
    for r in an.readouts:
        if r.kind == "point":
            ax.axvline(r.start, color="0.25", lw=0.8, ls=":")
            ax.text(r.start + 2, 0.003, f"{r.start:g} s", fontsize=6.5, color="0.25")
    ax.set_ylim(-0.01, 0.37)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("BLI response (nm)")
    ax.set_title("a  Sensor traces")
    ax.legend(fontsize=6, frameon=False, loc="center right", bbox_to_anchor=(1.0, 0.55))


def benchmark_panel(ax, df):
    sel = df[(df.readout == "auc_full")]
    for norm, c, m in (("absolute", ABS_C, "o"), ("relative", REL_C, "s")):
        d = sel[sel.normalization == norm]
        ax.scatter(d.true_ic50_uM, d.ic50_uM, s=12, color=c, marker=m, alpha=0.75, lw=0,
                   label=f"{norm} (median {d.fold_error.median():.2f}-fold)")
    lim = np.array([0.08, 5])
    ax.plot(lim, lim, color="0.4", lw=0.8)
    ax.fill_between(lim, lim / 1.25, lim * 1.25, color="0.9", zorder=0, label="±1.25-fold")
    ax.set_xscale("log"), ax.set_yscale("log")
    ax.set_xlim(*lim), ax.set_ylim(*lim)
    ax.set_xlabel("true IC50 (µM)")
    ax.set_ylabel("Find50 IC50 (µM)")
    ax.set_title(f"c  Synthetic benchmark (full AUC, n = {sel.compound.nunique()})")
    ax.legend(fontsize=6.5, frameon=False, loc="upper left")


def main():
    (HERE / "figures").mkdir(exist_ok=True)
    (HERE / "benchmark").mkdir(exist_ok=True)
    wells, _ = synthetic_plate({"Cmpd-1": 0.3}, seed=11)
    ws = [Well(i, w[0], parse_sample_id(w[1]), w[2], w[3]) for i, w in enumerate(wells)]
    an = analyze(ws, control="CTRL-1@2")
    curve = next(iter(an.curves))
    df = run(60)
    df.to_csv(HERE / "benchmark" / "benchmark_synthetic.csv", index=False)
    s = summary(df)
    s.to_csv(HERE / "benchmark" / "benchmark_summary.csv")

    fig, ax = plt.subplots(1, 3, figsize=(11, 3.1), gridspec_kw={"width_ratios": [1.25, 1, 1]})
    traces_panel(ax[0], an, curve)
    plot_dose_response(ax[1], an, curve, "auc_full")
    ax[1].set_title("b  Dose-response (full AUC), true IC50 0.3 µM")
    benchmark_panel(ax[2], df)
    fig.tight_layout()
    fig.savefig(HERE / "figures" / "method_overview.png", bbox_inches="tight")
    print(s.to_string())


if __name__ == "__main__":
    main()
