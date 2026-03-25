"""Configuration management for FusionSQL.

Supports loading from:
- Environment variables
- YAML config file (~/.fusql.yaml or ./fusql.yaml)
- Defaults

Settings:
- MSSQL connection string
- Default table names
- File patterns
- Log level
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, List
import os
import yaml


@dataclass
class MSSQLConfig:
    connection_string: Optional[str] = None
    driver: str = "ODBC Driver 17 for SQL Server"


@dataclass
class FusionSQLConfig:
    mssql: MSSQLConfig = field(default_factory=MSSQLConfig)
    table_ariba: str = "ariba_fusions"
    table_starfusion: str = "starfusion_fusions"
    table_concordance: str = "fusion_concordance"
    ariba_patterns: List[str] = field(
        default_factory=lambda: [r"^ariba_report\.tsv$", r"ariba.*report\.tsv$"]
    )
    starfusion_patterns: List[str] = field(
        default_factory=lambda: [
            r"^star-fusion\.fusion_predictions.*\.tsv$",
            r"^fusion_inspector.*\.tsv$",
        ]
    )
    log_level: str = "INFO"

    def to_dict(self) -> dict:
        """Convert config to a plain dict suitable for YAML serialization."""
        return {
            "mssql": asdict(self.mssql),
            "table_ariba": self.table_ariba,
            "table_starfusion": self.table_starfusion,
            "table_concordance": self.table_concordance,
            "ariba_patterns": self.ariba_patterns,
            "starfusion_patterns": self.starfusion_patterns,
            "log_level": self.log_level,
        }

    @classmethod
    def from_dict(cls, d: dict) -> "FusionSQLConfig":
        """Reconstruct FusionSQLConfig from a plain dict (e.g., from YAML)."""
        mssql = MSSQLConfig(**d.get("mssql", {}))
        return cls(
            mssql=mssql,
            table_ariba=d.get("table_ariba", "ariba_fusions"),
            table_starfusion=d.get("table_starfusion", "starfusion_fusions"),
            table_concordance=d.get("table_concordance", "fusion_concordance"),
            ariba_patterns=d.get(
                "ariba_patterns",
                [r"^ariba_report\.tsv$", r"ariba.*report\.tsv$"],
            ),
            starfusion_patterns=d.get(
                "starfusion_patterns",
                [
                    r"^star-fusion\.fusion_predictions.*\.tsv$",
                    r"^fusion_inspector.*\.tsv$",
                ],
            ),
            log_level=d.get("log_level", "INFO"),
        )


def _env_overrides(config: FusionSQLConfig) -> FusionSQLConfig:
    """Apply environment variable overrides to config.

    Environment variables (all prefixed with FUSQL_):
    - FUSQL_MSSQL: MSSQL connection string
    - FUSQL_DRIVER: ODBC driver name
    - FUSQL_TABLE_ARIBA: Ariba table name
    - FUSQL_TABLE_STARFUSION: STARFusion table name
    - FUSQL_TABLE_CONCORDANCE: Concordance table name
    - FUSQL_ARIBA_PATTERNS: JSON-encoded list of Ariba file patterns
    - FUSQL_STARFUSION_PATTERNS: JSON-encoded list of STARFusion file patterns
    - FUSQL_LOG_LEVEL: Logging level
    """
    import json

    if os.environ.get("FUSQL_MSSQL"):
        config.mssql.connection_string = os.environ["FUSQL_MSSQL"]
    if os.environ.get("FUSQL_DRIVER"):
        config.mssql.driver = os.environ["FUSQL_DRIVER"]
    if os.environ.get("FUSQL_TABLE_ARIBA"):
        config.table_ariba = os.environ["FUSQL_TABLE_ARIBA"]
    if os.environ.get("FUSQL_TABLE_STARFUSION"):
        config.table_starfusion = os.environ["FUSQL_TABLE_STARFUSION"]
    if os.environ.get("FUSQL_TABLE_CONCORDANCE"):
        config.table_concordance = os.environ["FUSQL_TABLE_CONCORDANCE"]
    if os.environ.get("FUSQL_ARIBA_PATTERNS"):
        config.ariba_patterns = json.loads(os.environ["FUSQL_ARIBA_PATTERNS"])
    if os.environ.get("FUSQL_STARFUSION_PATTERNS"):
        config.starfusion_patterns = json.loads(
            os.environ["FUSQL_STARFUSION_PATTERNS"]
        )
    if os.environ.get("FUSQL_LOG_LEVEL"):
        config.log_level = os.environ["FUSQL_LOG_LEVEL"]
    return config


def load_config(config_path: Optional[Path] = None) -> FusionSQLConfig:
    """Load config from env vars and config files.

    Priority (highest to lowest):
    1. Environment variables (FUSQL_*)
    2. Explicit config_path argument
    3. ./fusql.yaml (current directory)
    4. ~/.fusql.yaml (home directory)

    Args:
        config_path: Explicit path to a YAML config file. If provided, skips
                     the default search paths.

    Returns:
        FusionSQLConfig with all settings resolved.
    """
    config = FusionSQLConfig()

    # Search for config file if not explicitly provided
    search_paths: List[Path] = []
    if config_path is None:
        search_paths = [Path.cwd() / "fusql.yaml", Path.home() / ".fusql.yaml"]
    else:
        search_paths = [Path(config_path)]

    for path in search_paths:
        if path.exists():
            with open(path, "r") as f:
                data = yaml.safe_load(f)
            if data:
                config = FusionSQLConfig.from_dict(data)
            break

    # Environment variables always win
    config = _env_overrides(config)

    return config


def save_config(config: FusionSQLConfig, path: Path) -> None:
    """Save config to a YAML file.

    Args:
        config: FusionSQLConfig to serialize.
        path: Destination file path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(config.to_dict(), f, default_flow_style=False, sort_keys=False)
