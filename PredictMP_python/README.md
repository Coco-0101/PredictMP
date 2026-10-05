# PredictMP (Python)

Percentile rank of PrediXcan-predicted gene expression against simulated
gnomAD v4.1 reference populations (afr, ami, amr, asj, eas, fin, mid, nfe, sas)
across 49 GTEx v8 tissues.

PredictMP is the multi-population extension of PredictAP (Chan et al. 2024,
Brief Bioinform 25(6) bbae549), which provided an East Asian (EAS) reference only.

```bash
pip install -e .
```

```python
from PredictMP import rankcal
res = rankcal(data, tissue="x32", population="eas")
```

Reference data location: `PREDICTMP_DATA_DIR` (default `PredictMP/data/`),
built by `code/08_build_package_data.py`. Full pipeline and output columns:
see `code/README.md`.
