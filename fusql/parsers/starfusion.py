"""StarFusion and FusionInspector TSV parser for FusionSQL."""

import csv
import re
from pathlib import Path
from typing import List, Tuple, Optional

import pydantic

from fusql.utils.exceptions import ParserError
from fusql.utils.logging import get_logger

logger = get_logger(__name__)


class StarFusionEntry(pydantic.BaseModel):
    """Model for a single StarFusion/FusionInspector fusion prediction entry.
    
    Maps to the starfusion_fusions database table schema.
    """
    sample_id: str = pydantic.Field(description="Sample identifier")
    fusion_name: str = pydantic.Field(description="Fusion name (GENE1--GENE2)")
    gene1: str = pydantic.Field(description="Left gene symbol (extracted from LeftGene)")
    gene2: str = pydantic.Field(description="Right gene symbol (extracted from RightGene)")
    exon1: str = pydantic.Field(description="Exon info for gene1 (e.g., E10 from CDS_LEFT_RANGE)")
    exon2: str = pydantic.Field(description="Exon info for gene2 (e.g., E11 from CDS_RIGHT_RANGE)")
    splice_type: str = pydantic.Field(description="Splice type (INCL_NON_REF_SPLICE, ONLY_REF_SPLICE)")
    left_breakpos: int = pydantic.Field(description="Left genomic breakpoint position (numeric)")
    right_breakpos: int = pydantic.Field(description="Right genomic breakpoint position (numeric)")
    junction_reads: int = pydantic.Field(description="Junction-spanning read count")
    spanning_reads: int = pydantic.Field(description="Spanning fragment count")
    ffpm: float = pydantic.Field(description="Fusions per million FPMS")
    cds_left: str = pydantic.Field(description="CDS annotation for left gene (e.g., 1-923)")
    cds_right: str = pydantic.Field(description="CDS annotation for right gene (e.g., 1698-3270)")
    prot_fusion_type: str = pydantic.Field(description="Protein fusion type: INFRAME, FRAMESHIFT, .")
    large_anchor_support: str = pydantic.Field(description="Large anchor support: YES/NO")
    annots: str = pydantic.Field(description="Annotation categories (JSON-like string)")
    input_file: str = pydantic.Field(description="Source file path")

    class Config:
        """Pydantic config."""
        frozen = False


class StarFusionParser:
    """Parser for StarFusion/FusionInspector TSV output files.
    
    Real FusionInspector columns:
    #FusionName, JunctionReadCount, SpanningFragCount, est_J, est_S,
    LeftGene, LeftLocalBreakpoint, LeftBreakpoint, RightGene, RightLocalBreakpoint, RightBreakpoint,
    SpliceType, LargeAnchorSupport, NumCounterFusionLeft, NumCounterFusionRight,
    FAR_left, FAR_right, LeftBreakDinuc, LeftBreakEntropy, RightBreakDinuc, RightBreakEntropy,
    FFPM, microh_brkpt_dist, num_microh_near_brkpt, annots,
    CDS_LEFT_ID, CDS_LEFT_RANGE, CDS_RIGHT_ID, CDS_RIGHT_RANGE,
    PROT_FUSION_TYPE, FUSION_MODEL, FUSION_CDS, FUSION_TRANSL, PFAM_LEFT, PFAM_RIGHT
    
    Gene format: GENE^ENSG (e.g., "FIP1L1^ENSG00000145216.16")
    Breakpoint format: chr4:53414722:+ (chr:pos:strand)
    CDS_RANGE format: 1-923 (amino acid range)
    """
    
    def parse(self, file_path: Path) -> List[StarFusionEntry]:
        """Parse a StarFusion/FusionInspector TSV file.
        
        Args:
            file_path: Path to the TSV file
            
        Returns:
            List of StarFusionEntry pydantic models
            
        Raises:
            ParserError: If file cannot be read or parsed
        """
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                rows = list(reader)
        except Exception as e:
            raise ParserError(f"Failed to read StarFusion file {file_path}: {e}")
        
        fusions = []
        sample_id = self._extract_sample_id(file_path)
        
        for row in rows:
            fusion_name = row.get("#FusionName", "")
            
            if not fusion_name:
                continue
            
            # Parse gene names from LeftGene/RightGene columns
            gene1_raw = row.get("LeftGene", "")
            gene2_raw = row.get("RightGene", "")
            gene1 = self._extract_gene_symbol(gene1_raw)
            gene2 = self._extract_gene_symbol(gene2_raw)
            
            # Parse breakpoints
            left_breakpoint = row.get("LeftBreakpoint", "")
            right_breakpoint = row.get("RightBreakpoint", "")
            left_breakpos = self._parse_breakpoint(left_breakpoint)
            right_breakpos = self._parse_breakpoint(right_breakpoint)
            
            # Parse CDS ranges for exon info
            cds_left_range = row.get("CDS_LEFT_RANGE", "")
            cds_right_range = row.get("CDS_RIGHT_RANGE", "")
            exon1 = self._parse_cds_range(cds_left_range)
            exon2 = self._parse_cds_range(cds_right_range)
            
            # Parse numeric fields
            junction_reads = self._safe_int(row.get("JunctionReadCount", "0"))
            spanning_reads = self._safe_int(row.get("SpanningFragCount", "0"))
            ffpm = self._safe_float(row.get("FFPM", "0.0"))
            
            entry = StarFusionEntry(
                sample_id=sample_id,
                fusion_name=fusion_name,
                gene1=gene1,
                gene2=gene2,
                exon1=exon1,
                exon2=exon2,
                splice_type=row.get("SpliceType", "unknown"),
                left_breakpos=left_breakpos,
                right_breakpos=right_breakpos,
                junction_reads=junction_reads,
                spanning_reads=spanning_reads,
                ffpm=ffpm,
                cds_left=cds_left_range,
                cds_right=cds_right_range,
                prot_fusion_type=row.get("PROT_FUSION_TYPE", ".") or "unknown",
                large_anchor_support=row.get("LargeAnchorSupport", "NO"),
                annots=row.get("annots", ""),
                input_file=str(file_path),
            )
            fusions.append(entry)
        
        logger.info(f"Parsed {len(fusions)} fusions from {file_path}")
        return fusions
    
    def _extract_gene_symbol(self, gene_field: str) -> str:
        """Extract gene symbol from GENE^ENSG format.
        
        Args:
            gene_field: Raw gene field (e.g., "FIP1L1^ENSG00000145216.16")
            
        Returns:
            Gene symbol (e.g., "FIP1L1")
        """
        if not gene_field:
            return ""
        if "^" in gene_field:
            return gene_field.split("^")[0]
        return gene_field.strip()
    
    def _parse_breakpoint(self, breakpos_str: str) -> int:
        """Parse genomic breakpoint position from chr:pos:strand format.
        
        Args:
            breakpos_str: Breakpoint string like "chr4:53414722:+"
            
        Returns:
            Integer position (0 if unparseable)
        """
        if not breakpos_str or not breakpos_str.strip():
            return 0
        
        # Extract numeric position after colon
        match = re.search(r':(\d+):?', breakpos_str)
        if match:
            return int(match.group(1))
        
        # Try direct int conversion
        try:
            return int(breakpos_str.strip())
        except ValueError:
            return 0
    
    def _parse_cds_range(self, cds_range: str) -> str:
        """Parse CDS range to get exon-like representation.
        
        Args:
            cds_range: CDS range like "1-923" or "1698-3270"
            
        Returns:
            String representation like "E1-923" or just the range
        """
        if not cds_range or cds_range == ".":
            return ""
        
        if "-" in cds_range:
            parts = cds_range.split("-")
            try:
                start = int(parts[0])
                end = int(parts[1])
                return f"E{start}-{end}"
            except ValueError:
                return cds_range
        
        return cds_range
    
    def _safe_int(self, value: str) -> int:
        """Safely convert string to int."""
        try:
            return int(float(value.strip()))
        except (ValueError, TypeError):
            return 0
    
    def _safe_float(self, value: str) -> float:
        """Safely convert string to float."""
        try:
            return float(value.strip())
        except (ValueError, TypeError):
            return 0.0
    
    def _extract_sample_id(self, file_path: Path) -> str:
        """Extract sample ID from file path.
        
        Args:
            file_path: Path to the StarFusion file
            
        Returns:
            Sample ID string or filename stem if not found
        """
        path_str = str(file_path)
        
        # Common patterns
        patterns = [
            r"[/\\]Sample[_-]([^/\\]+)[/\\]",
            r"[/\\]sample[_-]([^/\\]+)[/\\]",
            r"[/\\]Run[_-]([^/\\]+)[/\\]",
            r"[/\\]([^/\\]+?)[/\\]starfusion",
            r"[/\\]([^/\\]+?)[/\\]fusion_inspector",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, path_str, re.IGNORECASE)
            if match:
                return match.group(1)
        
        return file_path.stem
