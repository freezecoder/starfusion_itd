"""Directory scanner with sample ID extraction and file grouping."""

import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .finder import find_ariba_files, find_starfusion_files


class FusionFileScanner:
    """
    Scans directories for Ariba and StarFusion fusion files,
    grouping them by run_id + sample_id.
    
    Expected directory structure:
        {run_id}/pipelineout/{sample_id}/{fusion_files}
    
    Example:
        Z3AT9/pipelineout/Z3AT9_IonCode_0125/star-fusion.fusion_predictions.abridged.tsv
    """
    
    def __init__(
        self,
        run_id_pattern: str = r"^([^/]+)",  # Top-level folder is run_id
        sample_id_pattern: str = r"([^/]+)$",  # Last folder is sample_id
        subdir_pattern: str = "pipelineout",  # Intermediate folder
    ):
        """
        Initialize scanner.
        
        Args:
            run_id_pattern: Regex with capture group for run_id (default: top-level folder)
            sample_id_pattern: Regex with capture group for sample_id (default: last folder)
            subdir_pattern: Expected subdirectory name between run_id and sample_id
        """
        self.run_id_pattern = run_id_pattern
        self.sample_id_pattern = sample_id_pattern
        self.subdir_pattern = subdir_pattern
    
    def extract_ids(self, file_path: Path) -> Tuple[Optional[str], Optional[str]]:
        """
        Extract run_id and sample_id from file path.
        
        Args:
            file_path: Path to a fusion file
            
        Returns:
            Tuple of (run_id, sample_id) or (None, None) if not matched
        """
        path_str = str(file_path)
        parts = path_str.replace("\\", "/").split("/")
        
        # Find the subdir (pipelineout) position
        subdir_idx = -1
        for i, part in enumerate(parts):
            if part == self.subdir_pattern:
                subdir_idx = i
                break
        
        if subdir_idx == -1 or subdir_idx >= len(parts) - 2:
            # No pipelineout found, or not enough parts
            return None, None
        
        run_id = parts[subdir_idx - 1] if subdir_idx > 0 else None
        sample_id = parts[subdir_idx + 1] if subdir_idx + 1 < len(parts) else None
        
        return run_id, sample_id
    
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
    
    def scan_directory(
        self, 
        root_dir: Path,
        ariba_patterns: Optional[List[str]] = None,
        starfusion_patterns: Optional[List[str]] = None,
    ) -> Dict[str, Dict[str, Path]]:
        """
        Scan directory for all Ariba and StarFusion files.
        
        Groups files by run_id + sample_id. Returns structure:
        {
            "run_id/sample_id": {
                "run_id": str,
                "sample_id": str,
                "ariba": Path or None,
                "starfusion": Path or None,
            }
        }
        
        Args:
            root_dir: Root directory to scan
            ariba_patterns: Custom regex patterns for Ariba files
            starfusion_patterns: Custom regex patterns for StarFusion files
            
        Returns:
            Dictionary mapping "run_id/sample_id" to their fusion files
        """
        ariba_files = find_ariba_files(root_dir, name_patterns=ariba_patterns)
        starfusion_files = find_starfusion_files(root_dir, name_patterns=starfusion_patterns)
        
        # Build mapping: "run_id/sample_id" -> {run_id, sample_id, ariba, starfusion}
        sample_map: Dict[str, Dict] = {}
        
        # Process Ariba files
        for fpath in ariba_files:
            run_id, sample_id = self.extract_ids(fpath)
            if run_id is None or sample_id is None:
                continue
            
            key = f"{run_id}/{sample_id}"
            if key not in sample_map:
                sample_map[key] = {
                    "run_id": run_id,
                    "sample_id": sample_id,
                    "ariba": None,
                    "starfusion": None,
                }
            sample_map[key]["ariba"] = fpath
        
        # Process StarFusion files
        for fpath in starfusion_files:
            run_id, sample_id = self.extract_ids(fpath)
            if run_id is None or sample_id is None:
                continue
            
            key = f"{run_id}/{sample_id}"
            if key not in sample_map:
                sample_map[key] = {
                    "run_id": run_id,
                    "sample_id": sample_id,
                    "ariba": None,
                    "starfusion": None,
                }
            sample_map[key]["starfusion"] = fpath
        
        return sample_map
    
    def find_samples(
        self, 
        root_dir: Path,
        ariba_patterns: Optional[List[str]] = None,
        starfusion_patterns: Optional[List[str]] = None,
    ) -> Dict[str, Dict]:
        """
        Alias for scan_directory for semantic clarity.
        
        Returns all discovered samples with their associated files.
        """
        return self.scan_directory(root_dir, ariba_patterns, starfusion_patterns)
