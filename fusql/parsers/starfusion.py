"""StarFusion and FusionInspector TSV parser."""

import csv
import re
from pathlib import Path
from typing import List, Tuple

import pydantic


class StarFusionEntry(pydantic.BaseModel):
    """Model for a single StarFusion fusion prediction entry."""
    
    sample_id: str
    fusion_name: str
    gene1: str
    gene2: str
    exon1: str
    exon2: str
    splice_type: str
    left_breakpos: int
    right_breakpos: int
    junction_reads: int
    spanning_reads: int
    ffpm: float
    cds_left: str
    cds_right: str
    input_file: str
    
    class Config:
        extra = "ignore"


class StarFusionParser:
    """
    Parser for StarFusion/FusionInspector TSV output files.
    
    Expected columns from StarFusion:
    - #FusionName: GENE1--GENE2 format
    - JunctionRead: count of junction-spanning reads
    - SpanningRead: count of spanning reads
    - SpliceType: SPLICE_SITE_TYPE
    - LeftBreakpoint: genomic position
    - RightBreakpoint: genomic position
    - FFPM: Fusions Per Million (FPKM)
    - CDS_LEFT: coding sequence annotation left
    - CDS_RIGHT: coding sequence annotation right
    
    Exon columns may contain formats like:
    - E10^E11 (exon 10 to exon 11)
    - E10 (single exon)
    """
    
    # Standard StarFusion column names
    COLUMNS = [
        "#FusionName",
        "JunctionRead",
        "SpanningRead",
        "SpliceType",
        "LeftBreakpoint",
        "RightBreakpoint",
        "FFPM",
        "CDS_LEFT",
        "CDS_RIGHT",
    ]
    
    def __init__(self, sample_id: str = "unknown"):
        """
        Initialize parser.
        
        Args:
            sample_id: Default sample ID to use if not extractable from file
        """
        self.default_sample_id = sample_id
    
    def extract_exons(self, exon_info: str) -> Tuple[str, str]:
        """
        Extract exon structure from StarFusion format.
        
        Handles formats like:
        - E10^E11 -> returns ("E10", "E11")
        - E10 -> returns ("E10", "")
        - ENST:123:E10^E11 -> returns ("E10", "E11")
        
        Args:
            exon_info: Exon information string from StarFusion
            
        Returns:
            Tuple of (exon1_str, exon2_str)
        """
        if not exon_info:
            return ("", "")
        
        # Handle formats like "E10^E11" or "E10"
        if "^" in exon_info:
            parts = exon_info.split("^")
            return (parts[0] if len(parts) > 0 else "", 
                    parts[1] if len(parts) > 1 else "")
        
        # Handle formats like "ENST:123:E10^E11" - extract after last colon
        if ":" in exon_info:
            last_part = exon_info.split(":")[-1]
            if "^" in last_part:
                parts = last_part.split("^")
                return (parts[0], parts[1] if len(parts) > 1 else "")
        
        # Single exon or unknown format
        match = re.search(r'E(\d+)', exon_info)
        if match:
            return (f"E{match.group(1)}", "")
        
        return ("", "")
    
    def parse(self, file_path: Path) -> List[StarFusionEntry]:
        """
        Parse a StarFusion TSV file.
        
        Args:
            file_path: Path to the StarFusion TSV file
            
        Returns:
            List of StarFusionEntry objects
            
        Raises:
            FileNotFoundError: If file does not exist
            ValueError: If required columns are missing
        """
        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        entries = []
        
        with open(file_path, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            
            # Validate columns
            if reader.fieldnames is None:
                raise ValueError(f"Empty or invalid TSV file: {file_path}")
            
            required_cols = {"#FusionName", "JunctionRead", "SpanningRead"}
            missing = required_cols - set(reader.fieldnames)
            if missing:
                raise ValueError(
                    f"Missing required columns in {file_path}: {missing}. "
                    f"Found: {reader.fieldnames}"
                )
            
            for row_num, row in enumerate(reader, start=2):
                try:
                    entry = self._parse_row(row, file_path)
                    entries.append(entry)
                except Exception as e:
                    # Log warning but continue processing
                    import logging
                    logging.warning(
                        f"Error parsing row {row_num} in {file_path}: {e}"
                    )
                    continue
        
        return entries
    
    def _parse_row(self, row: dict, file_path: Path) -> StarFusionEntry:
        """
        Parse a single row into a StarFusionEntry.
        
        Args:
            row: Dictionary of column -> value
            file_path: Source file path
            
        Returns:
            StarFusionEntry object
        """
        fusion_name = row.get("#FusionName", "")
        
        # Extract gene1 and gene2 from fusion name (GENE1--GENE2)
        genes = self._parse_fusion_name(fusion_name)
        
        # Parse exon information from CDS columns
        cds_left = row.get("CDS_LEFT", "")
        cds_right = row.get("CDS_RIGHT", "")
        
        exon1, exon2 = self._parse_exons_from_cds(cds_left, cds_right)
        
        # Parse numeric fields
        junction_reads = self._safe_int(row.get("JunctionRead", "0"))
        spanning_reads = self._safe_int(row.get("SpanningRead", "0"))
        ffpm = self._safe_float(row.get("FFPM", "0.0"))
        
        left_breakpos = self._parse_breakpoint(row.get("LeftBreakpoint", ""))
        right_breakpos = self._parse_breakpoint(row.get("RightBreakpoint", ""))
        
        return StarFusionEntry(
            sample_id=self.default_sample_id,
            fusion_name=fusion_name,
            gene1=genes[0],
            gene2=genes[1],
            exon1=exon1,
            exon2=exon2,
            splice_type=row.get("SpliceType", "unknown"),
            left_breakpos=left_breakpos,
            right_breakpos=right_breakpos,
            junction_reads=junction_reads,
            spanning_reads=spanning_reads,
            ffpm=ffpm,
            cds_left=cds_left,
            cds_right=cds_right,
            input_file=str(file_path),
        )
    
    def _parse_fusion_name(self, fusion_name: str) -> Tuple[str, str]:
        """
        Parse fusion name into gene1 and gene2.
        
        Args:
            fusion_name: Fusion name like "GENE1--GENE2"
            
        Returns:
            Tuple of (gene1, gene2)
        """
        if "--" in fusion_name:
            parts = fusion_name.split("--")
            return (parts[0].strip(), parts[1].strip() if len(parts) > 1 else "")
        return (fusion_name.strip(), "")
    
    def _parse_exons_from_cds(self, cds_left: str, cds_right: str) -> Tuple[str, str]:
        """
        Parse exon information from CDS annotations.
        
        CDS_LEFT and CDS_RIGHT contain information like:
        - ENST00000123456.3:c.100+50_200-75:E10^E11
        
        Args:
            cds_left: CDS_LEFT column value
            cds_right: CDS_RIGHT column value
            
        Returns:
            Tuple of (exon1, exon2)
        """
        exon1 = ""
        exon2 = ""
        
        # Extract exon info from CDS_LEFT (gene1 side)
        if cds_left:
            exon1, _ = self.extract_exons(cds_left)
        
        # Extract exon info from CDS_RIGHT (gene2 side)
        if cds_right:
            _, exon2 = self.extract_exons(cds_right)
        
        return (exon1, exon2)
    
    def _parse_breakpoint(self, breakpos_str: str) -> int:
        """
        Parse genomic breakpoint position.
        
        Args:
            breakpos_str: Breakpoint string like "chr1:12345" or just "12345"
            
        Returns:
            Integer position (0 if unparseable)
        """
        if not breakpos_str or breakpos_str.strip() == "":
            return 0
        
        # Try to extract numeric position
        match = re.search(r':(\d+)', breakpos_str)
        if match:
            return int(match.group(1))
        
        # Try direct int conversion
        try:
            return int(breakpos_str.strip())
        except ValueError:
            return 0
    
    def _safe_int(self, value: str) -> int:
        """Safely convert string to int."""
        try:
            return int(value.strip())
        except (ValueError, TypeError):
            return 0
    
    def _safe_float(self, value: str) -> float:
        """Safely convert string to float."""
        try:
            return float(value.strip())
        except (ValueError, TypeError):
            return 0.0
