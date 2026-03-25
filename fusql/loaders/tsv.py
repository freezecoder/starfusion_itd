"""CSV loader for fusql test mode.

Mirrors the exact schema of MSSQL tables.
"""

from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from fusql.loaders.base import BaseLoader


# CSV column order matching MSSQL schemas
TABLE_COLUMNS: Dict[str, List[str]] = {
    "ariba_fusions": [
        "id",
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
        "annotation",
        "input_file",
        "loaded_at",
    ],
    "starfusion_fusions": [
        "id",
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
        "CDS_left",
        "CDS_right",
        "input_file",
        "loaded_at",
    ],
    "fusion_concordance": [
        "id",
        "sample_id",
        "fusion_id",
        "gene1",
        "gene2",
        "ariba_found",
        "starfusion_found",
        "concordance_status",
        "ariba_exons",
        "starfusion_exons",
        "splice_concordant",
        "notes",
        "analyzed_at",
    ],
}


class CSVLoader(BaseLoader):
    """CSV loader for test mode.

    Writes rows to CSV files matching MSSQL schema exactly.
    """

    def __init__(self, output_dir: Path):
        """Initialize CSV loader.

        Args:
            output_dir: Directory to write CSV files to.
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._written_files: Dict[str, Path] = {}

    def write(self, table_name: str, rows: List[Dict[str, Any]]) -> Path:
        """Write rows to CSV file.

        Args:
            table_name: Name of the table (becomes filename).
            rows: List of row dictionaries.

        Returns:
            Path to the written CSV file.
        """
        import csv

        output_path = self.output_dir / f"{table_name}.csv"
        columns = TABLE_COLUMNS.get(
            table_name, list(rows[0].keys()) if rows else []
        )

        with open(output_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=columns, extrasaction="ignore")
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
        """Check if CSV file exists for table.

        Args:
            table: Table name.

        Returns:
            True if CSV file exists, False otherwise.
        """
        return (self.output_dir / f"{table}.csv").exists()

    def create_table(self, table: str, schema: Dict[str, str]) -> None:
        """Create empty CSV file with headers (no-op for CSV mode).

        Since CSV is just a file format, this creates an empty file
        with headers matching the schema.

        Args:
            table: Table name.
            schema: Dictionary mapping column names to types (ignored).
        """
        columns = TABLE_COLUMNS.get(table, list(schema.keys()))
        output_path = self.output_dir / f"{table}.csv"

        if not output_path.exists():
            with open(output_path, "w", newline="", encoding="utf-8") as f:
                import csv

                writer = csv.DictWriter(f, fieldnames=columns)
                writer.writeheader()

    def get_written_files(self) -> Dict[str, Path]:
        """Get mapping of table names to their written file paths.

        Returns:
            Dictionary mapping table names to file paths.
        """
        return dict(self._written_files)
