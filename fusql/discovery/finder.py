"""File finder functions for Ariba and StarFusion outputs."""

import re
from pathlib import Path
from typing import List, Optional


# Default file name patterns for each caller type
# Arriba: arriba_report.tsv
ARIBA_PATTERNS = [
    re.compile(r"^ariba_report\.tsv$", re.IGNORECASE),
    re.compile(r"ariba.*report\.tsv$", re.IGNORECASE),
]

# StarFusion: star-fusion.fusion_predictions.abridged.tsv
STARFUSION_PATTERNS = [
    re.compile(r"^star-fusion\.fusion_predictions.*\.tsv$", re.IGNORECASE),
    re.compile(r"^fusion_inspector.*\.tsv$", re.IGNORECASE),
]


def is_ariba_file(filename: str) -> bool:
    """Check if filename matches Ariba naming pattern."""
    return any(p.match(filename) for p in ARIBA_PATTERNS)


def is_starfusion_file(filename: str) -> bool:
    """Check if filename matches StarFusion naming pattern."""
    return any(p.match(filename) for p in STARFUSION_PATTERNS)


def find_ariba_files(
    root_dir: Path, 
    glob_pattern: str = "*.tsv",
    name_patterns: List[str] = None
) -> List[Path]:
    """
    Recursively find Ariba fusion report files.
    
    Args:
        root_dir: Root directory to search
        glob_pattern: Glob pattern for file matching (default: "*.tsv")
        name_patterns: List of regex patterns for matching filenames
        
    Returns:
        List of Path objects pointing to Ariba files
    """
    patterns = [re.compile(p) for p in (name_patterns or [])] or ARIBA_PATTERNS
    
    ariba_files = []
    for path in root_dir.rglob(glob_pattern):
        if path.is_file() and any(p.match(path.name) for p in patterns):
            ariba_files.append(path)
    return sorted(ariba_files)


def find_starfusion_files(
    root_dir: Path, 
    glob_pattern: str = "*.tsv",
    name_patterns: List[str] = None
) -> List[Path]:
    """
    Recursively find StarFusion/FusionInspector output files.
    
    Args:
        root_dir: Root directory to search
        glob_pattern: Glob pattern for file matching (default: "*.tsv")
        name_patterns: List of regex patterns for matching filenames
        
    Returns:
        List of Path objects pointing to StarFusion files
    """
    patterns = [re.compile(p) for p in (name_patterns or [])] or STARFUSION_PATTERNS
    
    starfusion_files = []
    for path in root_dir.rglob(glob_pattern):
        if path.is_file() and any(p.match(path.name) for p in patterns):
            starfusion_files.append(path)
    return sorted(starfusion_files)


def extract_sample_id(file_path: Path, pattern: Optional[str] = None) -> Optional[str]:
    """
    Extract sample ID from file path using regex pattern.
    
    Common patterns in file names:
    - Sample_([^/]+)
    - Run_([^/]+)
    - ([A-Za-z0-9_-]+)_(ariba|starfusion|report)
    
    Args:
        file_path: Path to extract sample ID from
        pattern: Optional regex pattern with a single capture group
        
    Returns:
        Sample ID string if found, None otherwise
    """
    if pattern is None:
        # Default pattern: capture sample-like identifiers
        default_pattern = r"([A-Za-z0-9][A-Za-z0-9_-]*[A-Za-z0-9])"
        pattern = default_pattern
    else:
        pattern = pattern

    try:
        regex = re.compile(pattern)
        # Search in the file path string
        match = regex.search(str(file_path))
        if match:
            return match.group(1)
    except re.error:
        pass
    
    return None
