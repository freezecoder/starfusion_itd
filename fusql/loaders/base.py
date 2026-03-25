"""Base loader class for fusql data loaders."""

from abc import ABC, abstractmethod
from typing import Any, Dict, List


class BaseLoader(ABC):
    """Abstract base class for data loaders.

    Subclasses must implement write(), table_exists(), and create_table().
    """

    @abstractmethod
    def write(self, table: str, rows: List[Dict[str, Any]]) -> Any:
        """Write rows to the specified table.

        Args:
            table: Table name to write to.
            rows: List of dictionaries representing rows.

        Returns:
            Implementation-specific result (e.g., path to CSV or row count).
        """
        ...

    @abstractmethod
    def table_exists(self, table: str) -> bool:
        """Check if the specified table exists.

        Args:
            table: Table name to check.

        Returns:
            True if table exists, False otherwise.
        """
        ...

    @abstractmethod
    def create_table(self, table: str, schema: Dict[str, str]) -> None:
        """Create a table with the given schema if it does not exist.

        Args:
            table: Table name to create.
            schema: Dictionary mapping column names to SQL type strings.
        """
        ...
