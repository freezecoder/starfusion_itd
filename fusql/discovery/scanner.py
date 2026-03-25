"""Directory scanner with sample ID extraction and file grouping."""

from pathlib import Path
from typing import Dict, Optional

from .finder import find_ariba_files, find_starfusion_files, extract_sample_id


class FusionFileScanner:
    """
    Scans directories for Ariba and StarFusion fusion files,
    grouping them by sample ID.
    """
    
    def __init__(self, sample_id_pattern: Optional[str] = None):
        """
        Initialize scanner.
        
        Args:
            sample_id_pattern: Optional regex pattern with a single capture group
                               to extract sample ID from file paths.
        """
        self.sample_id_pattern = sample_id_pattern
    
    def identify_file_type(self, file_path: Path) -> Optional[str]:
        """
        Identify whether a file is an Ariba or StarFusion file.
        
        Args:
            file_path: Path to the file
            
        Returns:
            "ariba", "starfusion", or None if unknown type
        """
        from .finder import is_ariba_file, is_starfusion_file
        
        name = file_path.name
        if is_ariba_file(name):
            return "ariba"
        elif is_starfusion_file(name):
            return "starfusion"
        return None
    
    def scan_directory(self, root_dir: Path) -> Dict[str, Dict[str, Path]]:
        """
        Scan directory for all Ariba and StarFusion files.
        
        Groups files by sample ID. Returns structure:
        {
            sample_id: {
                "ariba": Path,        # or None if not found
                "starfusion": Path,   # or None if not found
            }
        }
        
        Args:
            root_dir: Root directory to scan
            
        Returns:
            Dictionary mapping sample IDs to their fusion files
        """
        ariba_files = find_ariba_files(root_dir)
        starfusion_files = find_starfusion_files(root_dir)
        
        # Build mapping: sample_id -> {ariba: Path, starfusion: Path}
        sample_map: Dict[str, Dict[str, Path]] = {}
        
        # Process Ariba files
        for fpath in ariba_files:
            sample_id = extract_sample_id(fpath, self.sample_id_pattern)
            if sample_id is None:
                # Use the file stem as fallback identifier
                sample_id = fpath.stem
            
            if sample_id not in sample_map:
                sample_map[sample_id] = {"ariba": None, "starfusion": None}
            sample_map[sample_id]["ariba"] = fpath
        
        # Process StarFusion files
        for fpath in starfusion_files:
            sample_id = extract_sample_id(fpath, self.sample_id_pattern)
            if sample_id is None:
                sample_id = fpath.stem
            
            if sample_id not in sample_map:
                sample_map[sample_id] = {"ariba": None, "starfusion": None}
            sample_map[sample_id]["starfusion"] = fpath
        
        return sample_map
    
    def find_samples(self, root_dir: Path) -> Dict[str, Dict[str, Path]]:
        """
        Alias for scan_directory for semantic clarity.
        
        Returns all discovered sample IDs with their associated files.
        """
        return self.scan_directory(root_dir)
