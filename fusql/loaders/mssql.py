"""MSSQL loader for fusql using SQLAlchemy + pyodbc."""

from datetime import datetime
from typing import Any, Dict, List, Optional

from sqlalchemy import create_engine, text, inspect
from sqlalchemy.engine import Engine

from fusql.loaders.base import BaseLoader


class MSSQLLoader(BaseLoader):
    """MSSQL bulk loader using SQLAlchemy + pyodbc.

    Append-only mode: partitioned by run_id + sample_id.
    """

    def __init__(self, connection_string: str):
        """Initialize MSSQL loader.

        Args:
            connection_string: SQLAlchemy connection string for MSSQL.
                Format: mssql+pyodbc://user:pass@server/database?driver=ODBC+Driver+17+for+SQL+Server
        """
        self.connection_string = connection_string
        self._engine: Optional[Engine] = None

    @property
    def engine(self) -> Engine:
        """Get or create SQLAlchemy engine."""
        if self._engine is None:
            self._engine = create_engine(self.connection_string)
        return self._engine

    def write(self, table: str, rows: List[Dict[str, Any]]) -> int:
        """Write rows to table using bulk insert.

        Args:
            table: Table name.
            rows: List of row dictionaries.

        Returns:
            Number of rows inserted.
        """
        if not rows:
            return 0

        with self.engine.connect() as conn:
            for row in rows:
                row_copy = dict(row)
                
                # Handle datetime conversion
                for key, val in row_copy.items():
                    if isinstance(val, datetime):
                        row_copy[key] = val
                
                # Remove id if present (let DB generate it)
                row_copy.pop("id", None)
                
                # Build INSERT statement
                cols = list(row_copy.keys())
                placeholders = ", ".join([f":{c}" for c in cols])
                insert_sql = text(f"INSERT INTO {table} ({', '.join(cols)}) VALUES ({placeholders})")
                conn.execute(insert_sql, row_copy)
            
            conn.commit()
        
        return len(rows)

    def bulk_insert(self, table: str, rows: List[Dict[str, Any]]) -> int:
        """Bulk insert rows into table.

        Args:
            table: Table name.
            rows: List of row dictionaries.

        Returns:
            Number of rows inserted.
        """
        return self.write(table, rows)

    def table_exists(self, table: str) -> bool:
        """Check if table exists in database.

        Args:
            table: Table name.

        Returns:
            True if table exists, False otherwise.
        """
        inspector = inspect(self.engine)
        return table in inspector.get_table_names()

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
            if col_name == "id":
                columns.append(f"{col_name} {col_type}")
            else:
                columns.append(f"{col_name} {col_type}")

        create_sql = f"CREATE TABLE {table} ({', '.join(columns)})"

        with self.engine.connect() as conn:
            conn.execute(text(create_sql))
            conn.commit()

    def close(self) -> None:
        """Close the database engine."""
        if self._engine is not None:
            self._engine.dispose()
            self._engine = None

    def __enter__(self) -> "MSSQLLoader":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()
