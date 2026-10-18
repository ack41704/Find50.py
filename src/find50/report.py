"""Writing an analysis to Excel / CSV / PDF."""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import __version__


def write_results(analysis, out_dir, stem="find50", plots=True) -> dict:
    """Write <stem>_results.xlsx (sheets: IC50 summary, IC50 details, Dose-response,
    Wells, Settings), the two main tables as CSV and, with `plots`, <stem>_curves.pdf."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {"xlsx": out / f"{stem}_results.xlsx", "results_csv": out / f"{stem}_ic50.csv",
             "points_csv": out / f"{stem}_dose_response.csv"}
    settings = pd.DataFrame({"setting": ["software", *analysis.settings],
                             "value": [f"Find50 {__version__}",
                                       *map(str, analysis.settings.values())]})
    with pd.ExcelWriter(paths["xlsx"], engine="openpyxl") as xw:
        analysis.summary().to_excel(xw, sheet_name="IC50 summary")
        analysis.results.to_excel(xw, sheet_name="IC50 details", index=False)
        analysis.points.to_excel(xw, sheet_name="Dose-response", index=False)
        analysis.well_table.to_excel(xw, sheet_name="Wells", index=False)
        settings.to_excel(xw, sheet_name="Settings", index=False)
        for ws in xw.book.worksheets:
            ws.column_dimensions["A"].width = 34
            ws.freeze_panes = "B2"
    analysis.results.to_csv(paths["results_csv"], index=False)
    analysis.points.to_csv(paths["points_csv"], index=False)
    if plots and analysis.curves:
        from .plotting import write_pdf
        paths["pdf"] = write_pdf(analysis, out / f"{stem}_curves.pdf")
    return paths
