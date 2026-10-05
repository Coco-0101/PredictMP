"""
PredictMP: Gene expression percentile ranking based on PrediXcan predictions.
"""

from .functions import ecdf_fn, rankcal
from .loader import load_tissue, load_result_table, load_ks_table

__version__ = "1.1.0"
__all__ = ["ecdf_fn", "rankcal", "load_tissue", "load_result_table", "load_ks_table"]
