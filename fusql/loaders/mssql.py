"""MSSQL loader for fusql using pyodbc."""

from datetime import datetime
from typing import Any, Dict, List, Optional

import pyodbc

from fusql.loaders.base import BaseLoader


# Schema definitions for each table
TABLE_SCHEMAS: Dict[str, Dict[str, str]] = {
    "ariba_fusions": {
        "id": "BIGINT IDENTITY PRIMARY KEY",
        "sample_id": "VARCHAR(100)",
        "fusion_id": "VARCHAR(200)",
        "gene1": "VARCHAR(50)",
        "gene2": "VARCHAR(50)",
        "exon1": "INT",
        "exon2": "INT",
        "splice_site_1": "VARCHAR(10)",
        "splice_site_2": "VARCHAR(10)",
        "reads": "INT",
        "coding": "BIT",
        "genes_strand": "VARCHAR(10)",
        "annotation": "VARCHAR(500)",
        "input_file": "VARCHAR(500)",
        "loaded_at": "DATETIME",
    },
    "starfusion_fusions": {
        "id": "BIGINT IDENTITY PRIMARY KEY",
        "sample_id": "VARCHAR(100)",
        "fusion_name": "VARCHAR(300)",
        "gene1": "VARCHAR(50)",
        "gene2": "VARCHAR(50)",
        "exon1": "VARCHAR(20)",
        "exon2": "VARCHAR(20)",
        "splice_type": "VARCHAR(50)",
        "left_breakpos": "BIGINT",
        "right_breakpos": "BIGINT",
        "junction_reads": "INT",
        "spanning_reads": "INT",
        "FFPM": "FLOAT",
        "CDS_left": "VARCHAR(50)",
        "CDS_right": "VARCHAR(50)",
        "input_file": "VARCHAR(500)",
        "loaded_at": "DATETIME",
    },
    "fusion_concordance": {
        "id": "BIGINT IDENTITY PRIMARY KEY",
        "sample_id": "VARCHAR(100)",
        "fusion_id": "VARCHAR(300)",
        "gene1": "VARCHAR(50)",
        "gene2": "VARCHAR(50)",
        "ariba_found": "BIT",
        "starfusion_found": "BIT",
        "concordance_status": "VARCHAR(20)",
        "ariba_exons": "VARCHAR(50)",
        "starfusion_exons": "VARCHAR(50)",
        "splice_concordant": "BIT",
        "notes": "VARCHAR(500)",
        "analyzed_at": "DATETIME",
    },
}


class MSSQLLoader(BaseLoader):
    """MSSQL bulk loader using pyodbc.

    Supports bulk insert, upsert, and automatic table creation.
    """

    def __init__(self, connection_string: str):
        """Initialize MSSQL loader.

        Args:
            connection_string: ODBC connection string for MSSQL.
        """
        self.connection_string = connection_string
        self._conn: Optional[pyodbc.Connection] = None

    def _get_connection(self) -> pyodbc.Connection:
        """Get or create database connection."""
        if self._conn is None:
            self._conn = pyodbc.connect(self.connection_string)
        return self._conn

    def _ensure_table(self, table: str) -> None:
        """Ensure table exists, creating if necessary."""
        if not self.table_exists(table):
            self.create_table(table, TABLE_SCHEMAS.get(table, {}))

    def write(self, table: str, rows: List[Dict[str, Any]]) -> int:
        """Write rows to table using bulk insert.

        Args:
            table: Table name.
            rows: List of row dictionaries.

        Returns:
            Number of rows inserted.
        """
        self._ensure_table(table)
        if not rows:
            return 0

        conn = self._get_connection()
        cursor = conn.cursor()

        columns = list(rows[0].keys())
        placeholders = ", ".join(["?" for _ in columns])
        col_names = ", ".join(columns)

        batch_size = 1000
        total_inserted = 0

        for batch_start in range(0, len(rows), batch_size):
            batch = rows[batch_start : batch_start + batch_size]
            for row in batch:
                values = []
                for col in columns:
                    val = row.get(col)
                    if isinstance(val, datetime):
                        values.append(val)
                    elif val is None:
                        values.append(None)
                    else:
                        values.append(val)
                cursor.execute(
                    f"INSERT INTO {table} ({col_names}) VALUES ({placeholders})",
                    values,
                )
            conn.commit()
            total_inserted += len(batch)

        cursor.close()
        return total_inserted

    def bulk_insert(self, table: str, rows: List[Dict[str, Any]]) -> int:
        """Bulk insert rows into table.

        Uses batch INSERT for compatibility. Falls back to BULK INSERT
        if available and appropriate.

        Args:
            table: Table name.
            rows: List of row dictionaries.

        Returns:
            Number of rows inserted.
        """
        return self.write(table, rows)

    def upsert(
        self, table: str, rows: List[Dict[str, Any]], key_columns: List[str]
    ) -> int:
        """Insert or update rows on key conflict.

        Args:
            table: Table name.
            rows: List of row dictionaries.
            key_columns: Columns to use as unique key.

        Returns:
            Number of rows affected.
        """
        self._ensure_table(table)
        if not rows:
            return 0

        conn = self._get_connection()
        cursor = conn.cursor()

        columns = [c for c in rows[0].keys() if c not in ("id", "loaded_at", "analyzed_at")]
        col_names = ", ".join(columns)

        total_affected = 0

        for row in rows:
            set_clause = ", ".join([f"{c} = ?" for c in columns])
            where_clause = " AND ".join([f"{c} = ?" for c in key_columns])

            # Check if exists
            check_sql = f"SELECT 1 FROM {table} WHERE {where_clause}"
            check_vals = [row.get(c) for c in key_columns]
            cursor.execute(check_sql, check_vals)
            exists = cursor.fetchone() is not None

            if exists:
                # Update
                update_vals = [row.get(c) for c in columns]
                update_vals.extend([row.get(c) for c in key_columns])
                update_sql = f"UPDATE {table} SET {set_clause} WHERE {where_clause}"
                cursor.execute(update_sql, update_vals)
            else:
                # Insert
                placeholders = ", ".join(["?" for _ in columns])
                insert_vals = [row.get(c) for c in columns]
                insert_sql = f"INSERT INTO {table} ({col_names}) VALUES ({placeholders})"
                cursor.execute(insert_sql, insert_vals)

            total_affected += cursor.rowcount

        conn.commit()
        cursor.close()
        return total_affected

    def table_exists(self, table: str) -> bool:
        """Check if table exists in database.

        Args:
            table: Table name.

        Returns:
            True if table exists, False otherwise.
        """
        try:
            conn = self._get_connection()
            cursor = conn.cursor()
            cursor.execute(
                "SELECT 1 FROM sys.tables WHERE name = ? AND schema_id = SCHEMA_ID('dbo')",
                table,
            )
            exists = cursor.fetchone() is not None
            cursor.close()
            return exists
        except pyodbc.Error:
            return False

    def create_table(self, table: str, schema: Dict[str, str]) -> None:
        """Create table with given schema.

        Args:
            table: Table name.
            schema: Dictionary mapping column names to SQL types.
        """
        if not schema:
            raise ValueError(f"No schema provided for table {table}")

        columns = []
        for col_name, col_type in schema.items():
            if col_name != "id":
                columns.append(f"{col_name} {col_type}")

        create_sql = f"CREATE TABLE {table} ({', '.join(columns)})"

        conn = self._get_connection()
        cursor = conn.cursor()
        cursor.execute(create_sql)
        conn.commit()
        cursor.close()

    def close(self) -> None:
        """Close the database connection."""
        if self._conn is not None:
            self._conn.close()
            self._conn = None

    def __enter__(self) -> "MSSQLLoader":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()
