"""Tests for discovery module (finder and scanner)."""

import tempfile
from pathlib import Path

from fusql.discovery.finder import (
    find_ariba_files,
    find_starfusion_files,
    extract_sample_id,
    is_ariba_file,
    is_starfusion_file,
)
from fusql.discovery.scanner import FusionFileScanner


class TestFinder:
    """Tests for file finder functions."""

    def test_find_ariba_files(self, tmp_path):
        """find_ariba_files locates Ariba report files by pattern."""
        # Create a fake directory structure
        sample_dir = tmp_path / "Sample_123" / "pipelineout" / "sample_123"
        sample_dir.mkdir(parents=True)
        ariba_file = sample_dir / "ariba_report.tsv"
        ariba_file.write_text("#gene1\tgene2\nGENE1\tGENE2\n")

        results = find_ariba_files(tmp_path)
        assert len(results) == 1
        assert results[0].name == "ariba_report.tsv"

    def test_find_starfusion_files(self, tmp_path):
        """find_starfusion_files locates StarFusion output files by pattern."""
        sample_dir = tmp_path / "Sample_456" / "pipelineout" / "sample_456"
        sample_dir.mkdir(parents=True)
        sf_file = sample_dir / "star-fusion.fusion_predictions.abridged.tsv"
        sf_file.write_text("#FusionName\tJunctionReadCount\nGENE1--GENE2\t100\n")

        results = find_starfusion_files(tmp_path)
        assert len(results) == 1
        assert "star-fusion" in results[0].name

    def test_extract_sample_id(self):
        """extract_sample_id pulls sample ID from file path using default regex."""
        path = Path("/data/Sample_ABC/pipelineout/Sample_ABC/file.tsv")
        sample_id = extract_sample_id(path)
        assert sample_id == "Sample_ABC"

    def test_scanner_extracts_run_and_sample_ids(self, tmp_path):
        """FusionFileScanner.extract_ids returns correct run_id and sample_id."""
        scanner = FusionFileScanner()
        # Simulate path: {run_id}/pipelineout/{sample_id}/star-fusion.tsv
        file_path = Path("/Z3AT9/pipelineout/Z3AT9_IonCode_0125/star-fusion.fusion_predictions.abridged.tsv")
        run_id, sample_id = scanner.extract_ids(file_path)
        assert run_id == "Z3AT9"
        assert sample_id == "Z3AT9_IonCode_0125"
