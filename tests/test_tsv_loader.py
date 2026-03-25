"""Tests for TSV loader."""

import csv
from pathlib import Path

from fusql.loaders.tsv import TSVLoader, TABLE_COLUMNS


class TestTSVLoader:
    """Tests for TSVLoader.write()."""

    def test_write_creates_tsv_file(self, temp_output_dir):
        """write() creates a .tsv file in the output directory."""
        loader = TSVLoader(temp_output_dir)
        rows = [
            {
                "id": 1,
                "run_id": "RUN1",
                "sample_id": "SAMPLE1",
                "gene1": "BCR",
                "gene2": "ABL1",
            }
        ]
        output_path = loader.write("ariba_fusions", rows)
        assert output_path.exists()
        assert output_path.suffix == ".tsv"

    def test_tsv_has_correct_headers(self, temp_output_dir):
        """Written TSV has headers matching the table schema."""
        loader = TSVLoader(temp_output_dir)
        rows = [
            {
                "id": 1,
                "run_id": "RUN1",
                "sample_id": "SAMPLE1",
                "gene1": "BCR",
                "gene2": "ABL1",
            }
        ]
        loader.write("ariba_fusions", rows)
        output_path = loader.get_written_files()["ariba_fusions"]
        with open(output_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            headers = reader.fieldnames
        expected = TABLE_COLUMNS["ariba_fusions"]
        assert headers == expected

    def test_tsv_is_tab_separated(self, temp_output_dir):
        """Written file uses tab as the field delimiter."""
        loader = TSVLoader(temp_output_dir)
        rows = [{"id": 1, "run_id": "RUN1", "sample_id": "SAMPLE1", "gene1": "G1", "gene2": "G2"}]
        output_path = loader.write("ariba_fusions", rows)
        with open(output_path) as f:
            first_line = f.readline()
        # No commas in tab-separated output
        assert "\t" in first_line

    def test_tsv_has_run_id_sample_id(self, temp_output_dir):
        """Written TSV rows contain run_id and sample_id fields."""
        loader = TSVLoader(temp_output_dir)
        rows = [
            {
                "id": 1,
                "run_id": "RUN42",
                "sample_id": "SAMPLE99",
                "gene1": "GENE1",
                "gene2": "GENE2",
            }
        ]
        output_path = loader.write("ariba_fusions", rows)
        with open(output_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            written_rows = list(reader)
        assert len(written_rows) == 1
        assert written_rows[0]["run_id"] == "RUN42"
        assert written_rows[0]["sample_id"] == "SAMPLE99"
