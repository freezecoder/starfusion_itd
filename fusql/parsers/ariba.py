"""Ariba fusion parser for FusionSQL."""

import re
from pathlib import Path
from typing import List, Tuple

import pandas as pd
from pydantic import BaseModel, Field

from fusql.utils.exceptions import ParserError
from fusql.utils.logging import get_logger

logger = get_logger(__name__)


class AribaFusion(BaseModel):
    """Pydantic model for Ariba fusion records.
    
    Maps to the ariba_fusions database table schema.
    """
    sample_id: str = Field(description="Sample identifier")
    fusion_id: str = Field(description="Unique fusion identifier")
    gene1: str = Field(description="Left gene symbol")
    gene2: str = Field(description="Right gene symbol")
    exon1: int = Field(description="Exon number in gene1")
    exon2: int = Field(description="Exon number in gene2")
    splice_site_1: str = Field(description="Type: acceptor/donor/unknown")
    splice_site_2: str = Field(description="Type: acceptor/donor/unknown")
    reads: int = Field(description="Supporting read count")
    coding: bool = Field(description="Is fusion in coding region")
    genes_strand: str = Field(description="Strand orientation")
    annotation: str = Field(description="Full annotation string")
    input_file: str = Field(description="Source file path")

    class Config:
        """Pydantic config."""
        frozen = False


class AribaParser:
    """Parser for Ariba fusion report TSV files.
    
    Ariba reports contain gene fusion calls with supporting read counts,
    exon information, and splice site annotations.
    """
    
    REQUIRED_COLUMNS = [
        "gene1",
        "gene2", 
        "gene1_strand",
        "gene2_strand",
        "assembled_seqs",
        "coding",
    ]
    
    def parse(self, file_path: Path) -> List[AribaFusion]:
        """Parse an Ariba fusion TSV file.
        
        Args:
            file_path: Path to the Ariba TSV file
            
        Returns:
            List of AribaFusion pydantic models
            
        Raises:
            ParserError: If file cannot be read or parsed
        """
        try:
            df = pd.read_csv(file_path, sep="\t")
        except Exception as e:
            raise ParserError(f"Failed to read Ariba file {file_path}: {e}")
        
        # Validate required columns
        missing = set(self.REQUIRED_COLUMNS) - set(df.columns)
        if missing:
            raise ParserError(f"Missing required columns in {file_path}: {missing}")
        
        fusions = []
        sample_id = self._extract_sample_id(file_path)
        
        for _, row in df.iterrows():
            annotation = str(row.get("assembled_seqs", ""))
            exon1, exon2 = self.extract_exon_number(annotation)
            splice_site_1, splice_site_2 = self.extract_splice_site(annotation)
            
            # Build fusion_id from gene pair
            fusion_id = f"{row['gene1']}--{row['gene2']}"
            
            coding_val = row.get("coding", "no")
            is_coding = str(coding_val).lower() in ("yes", "true", "1")
            
            genes_strand = f"{row.get('gene1_strand', '.')}/{row.get('gene2_strand', '.')}"
            
            fusion = AribaFusion(
                sample_id=sample_id,
                fusion_id=fusion_id,
                gene1=str(row["gene1"]),
                gene2=str(row["gene2"]),
                exon1=exon1,
                exon2=exon2,
                splice_site_1=splice_site_1,
                splice_site_2=splice_site_2,
                reads=int(row.get("reads", 0)),
                coding=is_coding,
                genes_strand=genes_strand,
                annotation=annotation,
                input_file=str(file_path),
            )
            fusions.append(fusion)
        
        logger.info(f"Parsed {len(fusions)} fusions from {file_path}")
        return fusions
    
    def extract_exon_number(self, annotation: str) -> Tuple[int, int]:
        """Extract exon numbers from Ariba annotation string.
        
        Ariba annotations contain exon information in formats like:
        - "E10-E11" (exon 10 to exon 11)
        - "E10-E11intron" (with intron notation)
        - "E10" (single exon)
        
        Args:
            annotation: The assembled_seqs annotation string
            
        Returns:
            Tuple of (exon1, exon2) - defaults to (0, 0) if not found
        """
        # Match patterns like E10, E11, etc.
        exon_pattern = r"E(\d+)"
        matches = re.findall(exon_pattern, annotation)
        
        if len(matches) >= 2:
            return int(matches[0]), int(matches[1])
        elif len(matches) == 1:
            return int(matches[0]), 0
        else:
            return 0, 0
    
    def extract_splice_site(self, annotation: str) -> Tuple[str, str]:
        """Identify splice site types from annotation string.
        
        Splice sites can be:
        - acceptor: 5' splice site (start of intron)
        - donor: 3' splice site (end of intron) 
        - unknown: Cannot determine
        
        Args:
            annotation: The assembled_seqs annotation string
            
        Returns:
            Tuple of (splice_site_1, splice_site_2) for gene1 and gene2
        """
        annotation_lower = annotation.lower()
        
        # Determine splice site types based on annotation content
        # Common patterns in Ariba output
        
        # Check for exon-exon junction notation
        if "e" in annotation_lower and re.search(r"e\d+[-_]e\d+", annotation_lower):
            # This is an exon-exon junction
            # First gene is typically at acceptor/donor depending on strand
            splice1 = "acceptor"
            splice2 = "donor"
        elif re.search(r"(intron|donor|acceptor)", annotation_lower):
            # Explicit splice site notation
            splice1 = "unknown"
            splice2 = "unknown"
            if "donor" in annotation_lower:
                splice2 = "donor"
            if "acceptor" in annotation_lower:
                splice1 = "acceptor"
        else:
            splice1 = "unknown"
            splice2 = "unknown"
        
        return splice1, splice2
    
    def _extract_sample_id(self, file_path: Path) -> str:
        """Extract sample ID from file path.
        
        Attempts to extract sample ID using common patterns.
        
        Args:
            file_path: Path to the Ariba file
            
        Returns:
            Sample ID string or filename stem if not found
        """
        path_str = str(file_path)
        
        # Common patterns: Sample_XXXX, run_XXXX, sample_XXXX
        patterns = [
            r"[/\\]Sample[_-]([^/\\]+)[/\\]",
            r"[/\\]sample[_-]([^/\\]+)[/\\]",
            r"[/\\]Run[_-]([^/\\]+)[/\\]",
            r"[/\\]run[_-]([^/\\]+)[/\\]",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, path_str, re.IGNORECASE)
            if match:
                return match.group(1)
        
        # Fall back to filename stem
        return file_path.stem
