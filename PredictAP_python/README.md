# PredictAP (Python)

Percentile rank of PrediXcan-predicted gene expression against simulated
gnomAD v4.1 reference populations (afr, ami, amr, asj, eas, fin, mid, nfe, sas)
across 49 GTEx v8 tissues (Chan et al. 2024, Brief Bioinform 25(6) bbae549).

```bash
pip install -e .
```

```python
from PredictAP import rankcal
res = rankcal(data, tissue="x32", population="eas")
```

Reference data location: `PREDICTAP_DATA_DIR` (default `PredictMP/data/`),
built by `code/08_build_package_data.py`. Full pipeline and output columns:
see `code/README.md`.
