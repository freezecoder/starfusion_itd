"""Ariba fusion parser for FusionSQL."""

import re
from pathlib import Path
from typing import List, Tuple, Optional

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
    exon1: int = Field(description="Exon number in gene1 (0 if unknown)")
    exon2: int = Field(description="Exon number in gene2 (0 if unknown)")
    splice_site_1: str = Field(description="Splice site type for gene1 (acceptor/donor/CDS/splice-site/unknown)")
    splice_site_2: str = Field(description="Splice site type for gene2 (acceptor/donor/CDS/splice-site/unknown)")
    reads: int = Field(description="Supporting read count (split_reads1 + split_reads2)")
    coding: bool = Field(description="Is fusion in coding region (CDS)")
    genes_strand: str = Field(description="Strand orientation from strand1(gene/fusion)")
    reading_frame: str = Field(description="Reading frame: in-frame, out-of-frame, or .")
    confidence: str = Field(description="Arriba confidence level: high, medium, low")
    fusion_type: str = Field(description="Fusion type: translocation, deletion, etc.")
    input_file: str = Field(description="Source file path")

    class Config:
        """Pydantic config."""
        frozen = False


class AribaParser:
    """Parser for Arriba fusion report TSV files.
    
    Arriba reports contain gene fusion calls with supporting read counts,
    exon information, and splice site annotations.
    
    Real Arriba columns:
    #gene1, gene2, strand1(gene/fusion), strand2(gene/fusion), breakpoint1, breakpoint2,
    site1, site2, type, split_reads1, split_reads2, discordant_mates, coverage1, coverage2,
    confidence, reading_frame, tags, retained_protein_domains, closest_genomic_breakpoint1,
    closest_genomic_breakpoint2, gene_id1, gene_id2, transcript_id1, transcript_id2,
    direction1, direction2, filters, fusion_transcript, peptide_sequence, read_identifiers
    """
    
    # Column name mappings (file columns -> model fields)
    COLUMN_MAPPING = {
        "#gene1": "gene1",
        "gene1": "gene1", 
        "gene2": "gene2",
        "strand1(gene/fusion)": "strand1",
        "strand2(gene/fusion)": "strand2", 
        "site1": "site1",
        "site2": "site2",
        "type": "fusion_type",
        "split_reads1": "split_reads1",
        "split_reads2": "split_reads2", 
        "confidence": "confidence",
        "reading_frame": "reading_frame",
        "fusion_transcript": "fusion_transcript",
        "filters": "filters",
        "tags": "tags",
    }
    
    def parse(self, file_path: Path) -> List[AribaFusion]:
        """Parse an Arriba fusion TSV file.
        
        Args:
            file_path: Path to the Arriba TSV file
            
        Returns:
            List of AribaFusion pydantic models
            
        Raises:
            ParserError: If file cannot be read or parsed
        """
        try:
            # Arriba files have header starting with # (e.g. "#gene1")
            # Don't use comment="#" — that would skip the header row entirely.
            # Instead read normally and strip the leading '#' from the first column name.
            df = pd.read_csv(file_path, sep="\t", header=0)
            # Strip '#' prefix from the gene1 column if present
            if df.columns[0] == "#gene1":
                df.columns = df.columns.str.lstrip("#")
        except Exception as e:
            raise ParserError(f"Failed to read Arriba file {file_path}: {e}")
        
        fusions = []
        sample_id = self._extract_sample_id(file_path)
        
        for _, row in df.iterrows():
            gene1 = str(row.get("gene1", ""))
            gene2 = str(row.get("gene2", ""))
            
            if not gene1 or not gene2 or gene1 == "nan" or gene2 == "nan":
                continue
            
            # Build fusion_id from gene pair
            fusion_id = f"{gene1}--{gene2}"
            
            # Extract exon numbers from fusion_transcript
            fusion_transcript = str(row.get("fusion_transcript", ""))
            exon1, exon2 = self.extract_exon_number(fusion_transcript)
            
            # Get splice sites from site1/site2 columns
            site1 = str(row.get("site1", "unknown"))
            site2 = str(row.get("site2", "unknown"))
            splice_site_1 = self._normalize_splice_site(site1)
            splice_site_2 = self._normalize_splice_site(site2)
            
            # Calculate total reads
            split_reads1 = int(row.get("split_reads1", 0) or 0)
            split_reads2 = int(row.get("split_reads2", 0) or 0)
            total_reads = split_reads1 + split_reads2
            
            # Check coding status from site1/site2 (CDS means in coding region)
            is_coding = "CDS" in site1 or "CDS" in site2
            
            # Strand orientation
            strand1 = str(row.get("strand1(gene/fusion)", "./."))
            strand2 = str(row.get("strand2(gene/fusion)", "./."))
            genes_strand = f"{strand1}/{strand2}"
            
            # Reading frame
            reading_frame = str(row.get("reading_frame", "."))
            
            # Confidence level  
            confidence = str(row.get("confidence", "unknown"))
            
            # Fusion type
            fusion_type = str(row.get("type", "unknown"))
            
            fusion = AribaFusion(
                sample_id=sample_id,
                fusion_id=fusion_id,
                gene1=gene1,
                gene2=gene2,
                exon1=exon1,
                exon2=exon2,
                splice_site_1=splice_site_1,
                splice_site_2=splice_site_2,
                reads=total_reads,
                coding=is_coding,
                genes_strand=genes_strand,
                reading_frame=reading_frame if reading_frame != "." else "unknown",
                confidence=confidence,
                fusion_type=fusion_type,
                input_file=str(file_path),
            )
            fusions.append(fusion)
        
        logger.info(f"Parsed {len(fusions)} fusions from {file_path}")
        return fusions
    
    def extract_exon_number(self, text: str) -> Tuple[int, int]:
        """Extract exon numbers from fusion transcript or annotation.
        
        Arriba uses E# notation for exons in fusion_transcript:
        - "GACCATC...|ATCGTACT..." (transcript with | as junction)
        
        Args:
            text: Text containing exon information
            
        Returns:
            Tuple of (exon1, exon2) - defaults to (0, 0) if not found
        """
        # Match patterns like E10, E11 at the junction
        exon_pattern = r'E(\d+)'
        matches = re.findall(exon_pattern, text)
        
        if len(matches) >= 2:
            return int(matches[0]), int(matches[1])
        elif len(matches) == 1:
            return int(matches[0]), 0
        else:
            return 0, 0
    
    def _normalize_splice_site(self, site: str) -> str:
        """Normalize splice site to standard categories.
        
        Args:
            site: Raw site string from Arriba (e.g., "CDS/splice-site", "CDS", "splice-site")
            
        Returns:
            Normalized site type
        """
        site_lower = site.lower()
        
        if "splice-site" in site_lower or "splice_site" in site_lower:
            return "splice-site"
        elif "cds" in site_lower:
            return "CDS"
        elif "donor" in site_lower:
            return "donor"
        elif "acceptor" in site_lower:
            return "acceptor"
        else:
            return "unknown"
    
    def _extract_sample_id(self, file_path: Path) -> str:
        """Extract sample ID from file path.
        
        Attempts to extract sample ID using common patterns.
        
        Args:
            file_path: Path to the Arriba file
            
        Returns:
            Sample ID string or filename stem if not found
        """
        path_str = str(file_path)
        
        # Common patterns: Sample_XXXX, Run133, sample_XXXX
        patterns = [
            r"[/\\]Run(\d+)[/\\]",
            r"[/\\]Sample[_-]([^/\\]+)[/\\]",
            r"[/\\]sample[_-]([^/\\]+)[/\\]",
            r"[/\\]run[_-]([^/\\]+)[/\\]",
        ]
        
        for pattern in patterns:
            match = re.search(pattern, path_str, re.IGNORECASE)
            if match:
                return match.group(1)
        
        # Fall back to filename stem
        return file_path.stem
