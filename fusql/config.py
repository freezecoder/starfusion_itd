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
- Output templates (field selection)
"""

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional, List, Dict, Any
import os
import yaml


@dataclass
class MSSQLConfig:
    connection_string: Optional[str] = None
    driver: str = "ODBC Driver 17 for SQL Server"


@dataclass
class OutputTemplate:
    """Template for customizing output fields."""
    include: bool = True  # Whether to include this output
    fields: List[str] = None  # List of field names to include (None = all fields)
    
    def __post_init__(self):
        if self.fields is None:
            self.fields = []
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "include": self.include,
            "fields": self.fields if self.fields else [],
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> "OutputTemplate":
        if d is None:
            return cls()
        return cls(
            include=d.get("include", True),
            fields=d.get("fields", []),
        )


@dataclass
class OutputTemplates:
    """Collection of output templates for each output type."""
    ariba: OutputTemplate = field(default_factory=OutputTemplate)
    starfusion: OutputTemplate = field(default_factory=OutputTemplate)
    concordance: OutputTemplate = field(default_factory=OutputTemplate)
    
    def to_dict(self) -> Dict[str, Any]:
        return {
            "ariba": self.ariba.to_dict(),
            "starfusion": self.starfusion.to_dict(),
            "concordance": self.concordance.to_dict(),
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> "OutputTemplates":
        if d is None:
            return cls()
        return cls(
            ariba=OutputTemplate.from_dict(d.get("ariba", {})),
            starfusion=OutputTemplate.from_dict(d.get("starfusion", {})),
            concordance=OutputTemplate.from_dict(d.get("concordance", {})),
        )
    
    def filter_fields(self, output_type: str, all_fields: Dict[str, Any]) -> Dict[str, Any]:
        """Filter fields based on template.
        
        Args:
            output_type: 'ariba', 'starfusion', or 'concordance'
            all_fields: Dict of all field values
            
        Returns:
            Dict with only the selected fields (or all if no template specified)
        """
        template = getattr(self, output_type, None)
        if template is None or not template.fields:
            return all_fields
        
        return {k: v for k, v in all_fields.items() if k in template.fields}


# Default templates - show most useful fields
DEFAULT_ARIBA_FIELDS = [
    "run_id", "sample_id", "fusion_id", "gene1", "gene2",
    "exon1", "exon2", "splice_site_1", "splice_site_2",
    "reads", "coding", "reading_frame", "confidence"
]

DEFAULT_STARFUSION_FIELDS = [
    "run_id", "sample_id", "fusion_name", "gene1", "gene2",
    "exon1", "exon2", "splice_type", "junction_reads",
    "spanning_reads", "ffpm", "prot_fusion_type"
]

DEFAULT_CONCORDANCE_FIELDS = [
    "run_id", "sample_id", "fusion_id", "gene1", "gene2",
    "ariba_found", "starfusion_found", "concordance_status",
    "splice_concordant"
]


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
    output_templates: OutputTemplates = field(default_factory=OutputTemplates)

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
            "output_templates": self.output_templates.to_dict(),
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
            output_templates=OutputTemplates.from_dict(d.get("output_templates", {})),
        )


def _env_overrides(config: FusionSQLConfig) -> FusionSQLConfig:
    """Apply environment variable overrides to config."""
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
    """Load config from env vars and config files."""
    config = FusionSQLConfig()

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

    config = _env_overrides(config)
    return config


def save_config(config: FusionSQLConfig, path: Path) -> None:
    """Save config to a YAML file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        yaml.dump(config.to_dict(), f, default_flow_style=False, sort_keys=False)
