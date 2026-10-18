"""Figures: raw sensor traces and dose-response fits (matplotlib)."""
from __future__ import annotations

import numpy as np

ABS_C = "#2c6fbb"
REL_C = "#c0392b"


def _style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_traces(ax, analysis, curve, show_readouts=True):
    """Raw response vs time of a curve's dose wells (coloured by dose), with its
    reference (black) and positive-control (grey, dashed) wells."""
    import matplotlib.pyplot as plt

    by_index = {w.index: w for w in analysis.wells}
    doses = [by_index[i] for i in analysis.curves[curve]]
    cols = {w.column for w in doses}
    roles = dict(zip(analysis.well_table["well"], analysis.well_table["role"]))
    analyte = doses[0].sample.analyte
    refs = [w for w in analysis.wells if roles[w.index] in ("reference", "control")
            and w.sample.analyte == analyte and (w.column in cols or None in cols)]
    concs = sorted({w.sample.concentration_um for w in doses})
    cmap = plt.get_cmap("viridis")
    for w in refs:
        t = w.time - w.time[0]
        if roles[w.index] == "reference":
            ax.plot(t, w.response, color="black", lw=1.2, label="no inhibitor")
        else:
            ax.plot(t, w.response, color="0.55", lw=1.2, ls="--", label="positive control")
    for w in sorted(doses, key=lambda v: v.sample.concentration_um):
        k = concs.index(w.sample.concentration_um) / max(1, len(concs) - 1)
        ax.plot(w.time - w.time[0], w.response, color=cmap(0.9 * k), lw=1,
                label=f"{w.sample.concentration_um:g} µM")
    if show_readouts:
        for r in analysis.readouts:
            if r.kind == "point":
                ax.axvline(r.start, color="0.3", lw=0.7, ls=":")
    a_s = analysis.settings.get("association_s")
    if a_s:
        ax.axvspan(0, a_s, color="0.95", zorder=0)
    handles, labels = ax.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    ax.legend(uniq.values(), uniq.keys(), fontsize=6, frameon=False, ncol=1, loc="upper right")
    ax.set_xlabel("time (s)")
    ax.set_ylabel("response (nm)")
    ax.set_title("sensor traces", fontsize=9)
    _style(ax)


def plot_dose_response(ax, analysis, curve, readout, normalizations=("absolute", "relative")):
    """% activity vs concentration with the 4PL fit and the reported IC50."""
    pts = analysis.points
    lo, hi = np.inf, -np.inf
    for norm, colour, marker in (("absolute", ABS_C, "o"), ("relative", REL_C, "s")):
        if norm not in normalizations or (curve, readout, norm) not in analysis.fits:
            continue
        p = pts[(pts.curve == curve) & (pts.readout == readout) & (pts.normalization == norm)]
        res = analysis.fits[(curve, readout, norm)]
        c = p["concentration_uM"].to_numpy(float)
        lo, hi = min(lo, c.min()), max(hi, c.max())
        ax.scatter(c, p["pct_activity"], s=16, color=colour, marker=marker, zorder=3,
                   facecolor=colour if norm == "absolute" else "white", lw=1.1)
        if res.fit is not None:
            xs = np.linspace(np.log10(c.min()) - 0.3, np.log10(c.max()) + 0.3, 200)
            ax.plot(10 ** xs, res.fit(xs), color=colour, lw=1.2)
        if res.qualifier == "" and np.isfinite(res.ic50):
            ax.plot([res.ic50, res.ic50], [-5, 50], color=colour, lw=0.8, ls="--")
        ax.plot([], [], color=colour, marker=marker, mfc=colour if norm == "absolute" else "white",
                label=f"{norm}: IC50 {res.label} µM")
    ax.axhline(50, color="0.6", lw=0.7, ls=":")
    ax.set_xscale("log")
    if np.isfinite(lo):
        ax.set_xlim(lo / 2.5, hi * 2.5)
    ax.set_ylim(-10, 125)
    ax.set_xlabel("concentration (µM)")
    ax.set_ylabel("% activity")
    label = next((r.label for r in analysis.readouts if r.name == readout), readout)
    ax.set_title(label, fontsize=9)
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    _style(ax)


def curve_figure(analysis, curve):
    """One figure per curve: traces + one dose-response panel per readout."""
    import matplotlib.pyplot as plt

    n = len(analysis.readouts)
    fig, axes = plt.subplots(1, n + 1, figsize=(3.0 * (n + 1), 2.9))
    plot_traces(axes[0], analysis, curve)
    for ax, r in zip(axes[1:], analysis.readouts):
        plot_dose_response(ax, analysis, curve, r.name)
    fig.suptitle(curve, fontsize=10, fontweight="bold", x=0.01, ha="left")
    fig.tight_layout()
    return fig


def write_pdf(analysis, path):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        for curve in analysis.curves:
            fig = curve_figure(analysis, curve)
            pdf.savefig(fig)
            plt.close(fig)
    return path
