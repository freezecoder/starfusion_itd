"""TSV loader for fusql test mode.

Mirrors the exact schema of MSSQL tables.
"""

import csv
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List

from fusql.loaders.base import BaseLoader


# TSV column order matching MSSQL schemas
TABLE_COLUMNS: Dict[str, List[str]] = {
    "ariba_fusions": [
        "id",
        "run_id",
        "sample_id",
        "fusion_id",
        "gene1",
        "gene2",
        "exon1",
        "exon2",
        "splice_site_1",
        "splice_site_2",
        "reads",
        "coding",
        "genes_strand",
        "reading_frame",
        "confidence",
        "fusion_type",
        "input_file",
        "loaded_at",
    ],
    "starfusion_fusions": [
        "id",
        "run_id",
        "sample_id",
        "fusion_name",
        "gene1",
        "gene2",
        "exon1",
        "exon2",
        "splice_type",
        "left_breakpos",
        "right_breakpos",
        "junction_reads",
        "spanning_reads",
        "FFPM",
        "cds_left",
        "cds_right",
        "prot_fusion_type",
        "large_anchor_support",
        "annots",
        "input_file",
        "loaded_at",
    ],
    "fusion_concordance": [
        "id",
        "run_id",
        "sample_id",
        "fusion_id",
        "gene1",
        "gene2",
        "ariba_found",
        "starfusion_found",
        "concordance_status",
        "splice_concordant",
        "ariba_exons",
        "starfusion_exons",
        "notes",
        "analyzed_at",
    ],
}


class TSVLoader(BaseLoader):
    """TSV loader for test mode.

    Writes rows to tab-separated files matching MSSQL schema exactly.
    """

    def __init__(self, output_dir: Path):
        """Initialize TSV loader.

        Args:
            output_dir: Directory to write TSV files to.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._written_files: Dict[str, Path] = {}

    def write(self, table_name: str, rows: List[Dict[str, Any]]) -> Path:
        """Write rows to TSV file.

        Args:
            table_name: Name of the table (becomes filename).
            rows: List of row dictionaries.

        Returns:
            Path to the written TSV file.
        """
        output_path = self.output_dir / f"{table_name}.tsv"
        columns = TABLE_COLUMNS.get(
            table_name, list(rows[0].keys()) if rows else []
        )

        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f, fieldnames=columns, extrasaction="ignore", delimiter="\t"
            )
            writer.writeheader()

            for i, row in enumerate(rows, start=1):
                # Ensure id is set for identity column simulation
                row_copy = dict(row)
                if "id" in columns and "id" not in row_copy:
                    row_copy["id"] = i

                # Convert datetime objects to ISO format
                for key, val in row_copy.items():
                    if isinstance(val, datetime):
                        row_copy[key] = val.isoformat()

                writer.writerow(row_copy)

        self._written_files[table_name] = output_path
        return output_path

    def table_exists(self, table: str) -> bool:
        """Check if TSV file exists for table.

        Args:
            table: Table name.

        Returns:
            True if TSV file exists, False otherwise.
        """
        return (self.output_dir / f"{table}.tsv").exists()

    def create_table(self, table: str, schema: Dict[str, str]) -> None:
        """Create empty TSV file with headers.

        Args:
            table: Table name.
            schema: Dictionary mapping column names to types (ignored).
        """
        columns = TABLE_COLUMNS.get(table, list(schema.keys()))
        output_path = self.output_dir / f"{table}.tsv"

        if not output_path.exists():
            with open(output_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(
                    f, fieldnames=columns, delimiter="\t"
                )
                writer.writeheader()

    def get_written_files(self) -> Dict[str, Path]:
        """Get mapping of table names to their written file paths.

        Returns:
            Dictionary mapping table names to file paths.
        """
        return dict(self._written_files)
