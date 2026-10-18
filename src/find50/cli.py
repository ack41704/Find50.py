"""Command line: ``find50 octet | curve | demo``."""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from . import __version__


def build_parser():
    ap = argparse.ArgumentParser(prog="find50", description="IC50 values from dose-response and Octet BLI data")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = ap.add_subparsers(dest="command", required=True)

    o = sub.add_parser("octet", help="analyse an Octet export (.xlsx / .csv)")
    o.add_argument("file", type=Path)
    o.add_argument("--sheet", default=0, help="sheet name or index (0)")
    o.add_argument("--control", help="positive control NAME@CONC_uM, e.g. A20-I@2 (enables relative IC50s)")
    o.add_argument("--association-s", type=float,
                   help="association length in s (default: read from the Sample IDs, e.g. '(80 sec)')")
    o.add_argument("--assoc-time", type=float, help="association readout time in s (default association - 5)")
    o.add_argument("--dissoc-time", type=float, help="dissociation readout time in s (default association + 50)")
    o.add_argument("--exclude", nargs="*", default=[], help="well locations to leave out, e.g. C3 H5")
    o.add_argument("--reference-scope", choices=["column", "file"], default="column",
                   help="take reference / control wells from the same plate column (default) or the whole file")
    o.add_argument("--r2-min", type=float, default=0.8, help="minimum R^2 to report the 4PL IC50 (0.8)")
    o.add_argument("--out", type=Path, default=None, help="output folder (default: <file>_find50)")
    o.add_argument("--no-plots", action="store_true")

    c = sub.add_parser("curve", help="IC50 of one concentration series typed on the command line")
    c.add_argument("--conc", type=float, nargs="+", required=True, help="concentrations (> 0)")
    c.add_argument("--response", type=float, nargs="+", required=True, help="responses, same order")
    c.add_argument("--top", type=float, help="response of the no-inhibitor reference (100 %%)")
    c.add_argument("--bottom", type=float, default=0.0, help="response defined as 0 %% (default 0)")
    c.add_argument("--percent", action="store_true", help="responses are already %% activity")
    c.add_argument("--plot", type=Path, help="save the fit as an image")

    d = sub.add_parser("demo", help="analyse a synthetic Octet file with known IC50s")
    d.add_argument("--out", type=Path, default=Path("find50_demo"))
    return ap


def _curve(a):
    from .curve import ic50

    conc, y = np.asarray(a.conc, float), np.asarray(a.response, float)
    if conc.size != y.size:
        raise SystemExit("--conc and --response need the same number of values")
    if a.percent:
        pct = y
    else:
        top = a.top if a.top is not None else float(np.max(y))     # original: normalise to the max
        pct = 100.0 * (y - a.bottom) / (top - a.bottom)
    res = ic50(conc, pct)
    print(f"IC50 = {res.label}  ({res.method})")
    print(f"  4PL: IC50 {res.ic50_fit:.4g}, inflection {res.ec50_fit:.4g}, Hill {res.hill:.3g}, "
          f"top {res.top:.1f} %, bottom {res.bottom:.1f} %, R2 {res.r2:.4f}")
    print(f"  interpolated 50 % crossing: {res.ic50_interp:.4g}")
    for f in res.flags:
        print(f"  note: {f}")
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(3.6, 3.0))
        ax.scatter(conc, pct, color="#2c6fbb", zorder=3)
        if res.fit is not None:
            xs = np.linspace(np.log10(conc.min()) - 0.3, np.log10(conc.max()) + 0.3, 200)
            ax.plot(10 ** xs, res.fit(xs), color="#2c6fbb")
        ax.axhline(50, color="0.6", ls=":")
        ax.set_xscale("log")
        ax.set_xlabel("concentration")
        ax.set_ylabel("% activity")
        ax.set_title(f"IC50 {res.label}")
        fig.tight_layout()
        fig.savefig(a.plot, dpi=200)
    return 0


def main(argv=None) -> int:
    a = build_parser().parse_args(argv)
    if a.command == "curve":
        return _curve(a)

    import matplotlib
    matplotlib.use("Agg")
    from .analysis import analyze_file
    from .octet import default_readouts
    from .report import write_results

    if a.command == "demo":
        from .synthetic import write_demo_file
        path, truth = write_demo_file(a.out / "synthetic_octet.xlsx")
        an = analyze_file(path, control="CTRL-1@2")
        write_results(an, a.out, "synthetic")
        print(an.summary().to_string())
        print("\ntrue IC50 (uM):", ", ".join(f"{k} {v:g}" for k, v in truth.items()))
        print(f"\nwrote {a.out}/")
        return 0

    sheet = int(a.sheet) if str(a.sheet).isdigit() else a.sheet
    readouts = None
    if a.assoc_time is not None or a.dissoc_time is not None:
        readouts = default_readouts(a.association_s, a.assoc_time, a.dissoc_time)
    an = analyze_file(a.file, sheet=sheet, control=a.control, readouts=readouts,
                      association_s=a.association_s, exclude=a.exclude,
                      reference_scope=a.reference_scope, r2_min=a.r2_min)
    if an.settings["association_s"] is None and readouts is None:
        print("note: no association time found in the Sample IDs; only the full-trace AUC is used "
              "(pass --association-s)")
    out = a.out or a.file.with_name(a.file.stem + "_find50")
    paths = write_results(an, out, a.file.stem, plots=not a.no_plots)
    print(an.summary().to_string())
    print(f"\nwrote {paths['xlsx']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
