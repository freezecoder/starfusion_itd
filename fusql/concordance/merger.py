"""Fusion merger for combining Ariba and StarFusion results."""

from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

from pydantic import BaseModel, Field


class FusionConcordance(BaseModel):
    """Concordance record for a fusion call.

    Represents the merged view of a fusion found in either or both
    of the Ariba and StarFusion callers.
    """

    sample_id: str
    fusion_id: str
    gene1: str
    gene2: str
    ariba_found: bool = False
    starfusion_found: bool = False
    concordance_status: str = Field(
        description="shared/unique_ariba/unique_starfusion"
    )
    ariba_exons: Optional[str] = None
    starfusion_exons: Optional[str] = None
    splice_concordant: bool = False
    notes: Optional[str] = None
    analyzed_at: datetime = Field(default_factory=datetime.now)


def normalize_fusion_id(gene1: str, gene2: str) -> str:
    """Normalize fusion ID by sorting genes alphabetically.

    Fusions are represented as GENE1--GENE2 where GENE1 < GENE2
    alphabetically. This ensures consistent matching regardless of
    the order genes were originally reported.

    Args:
        gene1: First gene symbol.
        gene2: Second gene symbol.

    Returns:
        Normalized fusion ID in format GENE1--GENE2.
    """
    gene1_upper = gene1.upper()
    gene2_upper = gene2.upper()

    if gene1_upper < gene2_upper:
        return f"{gene1_upper}--{gene2_upper}"
    else:
        return f"{gene2_upper}--{gene1_upper}"


def _get_exon_string(gene1: str, gene2: str, exon1: Any, exon2: Any) -> str:
    """Build exon string from exon components.

    Args:
        gene1: Left gene symbol.
        gene2: Right gene symbol.
        exon1: Exon for gene1.
        exon2: Exon for gene2.

    Returns:
        Exon string in format E<exon1>@<gene1>-E<exon2>@<gene2>.
    """
    e1 = str(exon1) if exon1 is not None else "?"
    e2 = str(exon2) if exon2 is not None else "?"
    return f"E{e1}@{gene1}-E{e2}@{gene2}"


def _get_starfusion_exon_string(exon1: str, exon2: str) -> str:
    """Build exon string from StarFusion format.

    Args:
        exon1: Exon structure (e.g., "E10").
        exon2: Exon structure (e.g., "E11").

    Returns:
        Exon string.
    """
    return f"{exon1}-{exon2}"


def merge_fusions(
    ariba_fusions: List[Dict[str, Any]],
    starfusion_fusions: List[Dict[str, Any]],
) -> List[FusionConcordance]:
    """Merge fusion calls from Ariba and StarFusion.

    Normalizes fusion IDs and classifies each fusion as:
    - shared: found by both callers
    - unique_ariba: found only by Ariba
    - unique_starfusion: found only by StarFusion

    Args:
        ariba_fusions: List of Ariba fusion dictionaries.
        starfusion_fusions: List of StarFusion fusion dictionaries.

    Returns:
        List of FusionConcordance records.
    """
    # Build lookup for Ariba fusions by normalized ID
    ariba_by_fusion: Dict[str, Dict[str, Any]] = {}
    for fusion in ariba_fusions:
        gene1 = fusion.get("gene1", "")
        gene2 = fusion.get("gene2", "")
        fusion_id = normalize_fusion_id(gene1, gene2)
        ariba_by_fusion[fusion_id] = fusion

    # Build lookup for StarFusion fusions by normalized ID
    starfusion_by_fusion: Dict[str, Dict[str, Any]] = {}
    for fusion in starfusion_fusions:
        gene1 = fusion.get("gene1", "")
        gene2 = fusion.get("gene2", "")
        fusion_id = normalize_fusion_id(gene1, gene2)
        starfusion_by_fusion[fusion_id] = fusion

    # Get all unique fusion IDs
    all_fusion_ids = set(ariba_by_fusion.keys()) | set(starfusion_by_fusion.keys())

    results: List[FusionConcordance] = []
    now = datetime.now()

    for fusion_id in sorted(all_fusion_ids):
        ariba_fusion = ariba_by_fusion.get(fusion_id)
        starfusion_fusion = starfusion_by_fusion.get(fusion_id)

        # Extract gene info
        if ariba_fusion:
            gene1 = ariba_fusion.get("gene1", "")
            gene2 = ariba_fusion.get("gene2", "")
            sample_id = ariba_fusion.get("sample_id", "")
        elif starfusion_fusion:
            gene1 = starfusion_fusion.get("gene1", "")
            gene2 = starfusion_fusion.get("gene2", "")
            sample_id = starfusion_fusion.get("sample_id", "")
        else:
            continue

        # Determine concordance status
        ariba_found = ariba_fusion is not None
        starfusion_found = starfusion_fusion is not None

        if ariba_found and starfusion_found:
            status = "shared"
        elif ariba_found:
            status = "unique_ariba"
        else:
            status = "unique_starfusion"

        # Build exon strings
        ariba_exons = None
        starfusion_exons = None
        splice_concordant = False

        if ariba_fusion:
            exon1 = ariba_fusion.get("exon1")
            exon2 = ariba_fusion.get("exon2")
            if exon1 is not None or exon2 is not None:
                ariba_exons = _get_exon_string(gene1, gene2, exon1, exon2)

        if starfusion_fusion:
            exon1 = starfusion_fusion.get("exon1", "")
            exon2 = starfusion_fusion.get("exon2", "")
            if exon1 or exon2:
                starfusion_exons = _get_starfusion_exon_string(
                    str(exon1), str(exon2)
                )

        # Check splice concordance (simplified: both callers at splice sites)
        if ariba_exons and starfusion_exons:
            # Extract exon numbers for comparison
            splice_concordant = True  # Placeholder - real impl would compare

        concord = FusionConcordance(
            sample_id=sample_id,
            fusion_id=fusion_id,
            gene1=gene1,
            gene2=gene2,
            ariba_found=ariba_found,
            starfusion_found=starfusion_found,
            concordance_status=status,
            ariba_exons=ariba_exons,
            starfusion_exons=starfusion_exons,
            splice_concordant=splice_concordant,
            analyzed_at=now,
        )
        results.append(concord)

    return results
