# Changelog

## 1.0.0 (2026-10-18)

Packaged rewrite of the original notebook, for publication.

- `find50` Python package and command line (`octet`, `curve`, `demo`).
- Octet export reader for `.xlsx` and `.csv`. Plate layout, doses, association time and
  controls are parsed from the Sample IDs; nothing is hard-coded.
- Four readouts per trace on the true time axis: association point, dissociation point,
  full AUC, association AUC. Each gets an absolute and a positive-control-relative
  normalisation, with references taken from the same plate column.
- 4PL fit with a free Hill slope: absolute IC50 (the 50 % crossing), inflection point,
  R² and 95 % CI. An exact interpolated crossing (the original definition) is also
  reported. Out-of-range results are reported as `< c_min` / `> c_max`. Non-monotonic
  curves and poor fits are flagged.
- Excel workbook (summary in the original layout, details, points, wells, settings), CSV, PDF of traces and fits.
- Tests, including reproduction of the original notebook's recorded output; synthetic benchmark; figures; citation file.

## 0.1 (2025-05-30)

Original notebook `IC50 FINAL PROGRAM-5.ipynb`, now in `legacy/` with its outputs cleared.
