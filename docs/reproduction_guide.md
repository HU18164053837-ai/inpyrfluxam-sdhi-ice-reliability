# Reproduction guide

Use Python 3.12 (the tested interpreter version is recorded in the validation log).
From the repository root:

```bash
python -m pip install -r requirements.txt
python data_and_code/scripts/test_reproduce_review.py
python data_and_code/scripts/reproduce_review.py
python scripts/check_checksums.py
python data_and_code/scripts/figures.py
python data_and_code/scripts/supplementary_figures.py
```

Statistics and comparator tests use only the standard library; plotting uses
NumPy, pandas and matplotlib. Reproduction writes only to the ignored directory
data_and_code/reproduction_check. It compares all 12 generated result tables with
frozen inputs using rel_tol=1e-12 and abs_tol=1e-12 for finite floating fields.
Identity fields, category fields, frozen integer fields, row order and column
schema remain exact. The identity_median_FE metric is numerical, not an identity.

The plotting commands overwrite figures/Fig1–3 and FigS1–S2 in their PDF, SVG and
PNG variants. Run plotting in a disposable clone/copy when preserving the exact
published file hashes. Plotting-library versions and PDF metadata can change
rendered bytes without changing scientific data. checksums must pass immediately
after cloning; after regenerated figures replace frozen assets they can differ.
The supplied PDF/SVG/PNG files are the frozen research figures, not newly edited
versions. Tests were run in a clean copy; no source/frozen files were rewritten.

The preserved directory layout is intentional: data_and_code/inputs is the
authoritative frozen input location, supplementary/tables holds archived research
tables and figures contains the frozen visuals. Original relative-path imports
to optional project tools are inert when absent; installed dependencies suffice.
The comparator and source scripts retain v1.0.2 logic without formula/path edits.

