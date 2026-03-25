"""Fusion concordance analysis for fusql.

Provides merging and analysis of Ariba and StarFusion fusion calls.
"""

from fusql.concordance.merger import (
    FusionConcordance,
    merge_fusions,
    normalize_fusion_id,
)
from fusql.concordance.analyzer import (
    calculate_concordance_stats,
    generate_summary_report,
)

__all__ = [
    "FusionConcordance",
    "merge_fusions",
    "normalize_fusion_id",
    "calculate_concordance_stats",
    "generate_summary_report",
]
