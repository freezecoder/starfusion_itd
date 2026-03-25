"""Concordance analyzer for fusion results."""

from typing import Any, Dict, List

from fusql.concordance.merger import FusionConcordance


def calculate_concordance_stats(merged: List[FusionConcordance]) -> Dict[str, Any]:
    """Calculate concordance statistics from merged fusion results.

    Args:
        merged: List of FusionConcordance records.

    Returns:
        Dictionary containing concordance statistics:
        - total: total number of unique fusions
        - shared: fusions found by both callers
        - unique_ariba: fusions only found by Ariba
        - unique_starfusion: fusions only found by StarFusion
        - splice_concordant: number with splice concordance
        - splice_concordance_rate: proportion with splice concordance
        - ariba_total: total fusions in Ariba
        - starfusion_total: total fusions in StarFusion
        - overlap_rate: proportion of Ariba fusions also in StarFusion
    """
    if not merged:
        return {
            "total": 0,
            "shared": 0,
            "unique_ariba": 0,
            "unique_starfusion": 0,
            "splice_concordant": 0,
            "splice_concordance_rate": 0.0,
            "ariba_total": 0,
            "starfusion_total": 0,
            "overlap_rate": 0.0,
        }

    total = len(merged)
    shared = sum(1 for f in merged if f.concordance_status == "shared")
    unique_ariba = sum(1 for f in merged if f.concordance_status == "unique_ariba")
    unique_starfusion = sum(
        1 for f in merged if f.concordance_status == "unique_starfusion"
    )
    splice_concordant = sum(1 for f in merged if f.splice_concordant)

    # Count totals per caller
    ariba_total = sum(1 for f in merged if f.ariba_found)
    starfusion_total = sum(1 for f in merged if f.starfusion_found)

    # Calculate rates
    splice_concordance_rate = (
        splice_concordant / shared if shared > 0 else 0.0
    )
    overlap_rate = shared / ariba_total if ariba_total > 0 else 0.0

    return {
        "total": total,
        "shared": shared,
        "unique_ariba": unique_ariba,
        "unique_starfusion": unique_starfusion,
        "splice_concordant": splice_concordant,
        "splice_concordance_rate": round(splice_concordance_rate, 4),
        "ariba_total": ariba_total,
        "starfusion_total": starfusion_total,
        "overlap_rate": round(overlap_rate, 4),
    }


def generate_summary_report(
    merged: List[FusionConcordance], stats: Dict[str, Any]
) -> str:
    """Generate human-readable summary report.

    Args:
        merged: List of FusionConcordance records.
        stats: Statistics dictionary from calculate_concordance_stats.

    Returns:
        Formatted summary string.
    """
    lines = [
        "=" * 60,
        "FUSION CONCORDANCE SUMMARY",
        "=" * 60,
        "",
        f"Total Unique Fusions:     {stats['total']}",
        "",
        "Breakdown by Caller:",
        f"  - Shared (both):        {stats['shared']}",
        f"  - Unique to Ariba:      {stats['unique_ariba']}",
        f"  - Unique to StarFusion: {stats['unique_starfusion']}",
        "",
        "Caller Totals:",
        f"  - Ariba fusions:        {stats['ariba_total']}",
        f"  - StarFusion fusions:   {stats['starfusion_total']}",
        "",
        "Concordance Metrics:",
        f"  - Splice concordance:    {stats['splice_concordant']} of shared",
        f"  - Splice rate:          {stats['splice_concordance_rate']:.1%}",
        f"  - Overlap rate:         {stats['overlap_rate']:.1%}",
        "",
    ]

    # Add details on shared fusions if any
    if stats["shared"] > 0:
        shared_fusions = [
            f for f in merged if f.concordance_status == "shared"
        ]
        lines.append("Shared Fusions (both callers):")
        lines.append("-" * 40)
        for fusion in shared_fusions[:20]:  # Limit to first 20
            lines.append(
                f"  {fusion.fusion_id} "
                f"(Ariba: {fusion.ariba_exons}, "
                f"StarFusion: {fusion.starfusion_exons})"
            )
        if len(shared_fusions) > 20:
            lines.append(f"  ... and {len(shared_fusions) - 20} more")
        lines.append("")

    # Add unique Ariba fusions
    if stats["unique_ariba"] > 0:
        unique_ariba = [
            f for f in merged if f.concordance_status == "unique_ariba"
        ]
        lines.append(f"Unique to Ariba ({len(unique_ariba)}):")
        lines.append("-" * 40)
        for fusion in unique_ariba[:10]:  # Limit to first 10
            lines.append(f"  {fusion.fusion_id} ({fusion.ariba_exons})")
        if len(unique_ariba) > 10:
            lines.append(f"  ... and {len(unique_ariba) - 10} more")
        lines.append("")

    # Add unique StarFusion fusions
    if stats["unique_starfusion"] > 0:
        unique_sf = [
            f for f in merged if f.concordance_status == "unique_starfusion"
        ]
        lines.append(f"Unique to StarFusion ({len(unique_sf)}):")
        lines.append("-" * 40)
        for fusion in unique_sf[:10]:  # Limit to first 10
            lines.append(f"  {fusion.fusion_id} ({fusion.starfusion_exons})")
        if len(unique_sf) > 10:
            lines.append(f"  ... and {len(unique_sf) - 10} more")
        lines.append("")

    lines.append("=" * 60)
    return "\n".join(lines)
