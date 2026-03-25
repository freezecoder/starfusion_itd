# FusionSQL — Gene Fusion Analysis to Microsoft SQL Server

**Project:** `starfusion_itd` monorepo  
**Purpose:** Parse gene fusion callers (Arriba, STARFusion+FusionInspector) → load into MSSQL, merge results, analyze concordance  
**Testing Mode:** TSV output mirroring MSSQL schemas  

## Input Sources

### 1. Arriba Fusion Caller
- **Tool:** Arriba (https://github.com/suhrig/arriba)
- **Output:** `arriba_fusions.tsv` — gene fusion calls with splice sites, reading frame, coverage

### 2. STARFusion + FusionInspector Pipeline
- **Tool:** STAR-Fusion (https://github.com/STAR-Fusion/STAR-Fusion) + FusionInspector (https://github.com/FusionInspector/FusionInspector)
- **Workflow:** STARFusion → candidate fusions → FusionInspector → validated fusions with deep annotation
- **Output:** `fusion_inspector.fusions.tsv` — STARFusion results validated by FusionInspector
- **Note:** STARFusion and FusionInspector are from the same tool suite; FusionInspector extends STARFusion analysis  

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                        FusionSQL Package                             │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────────────┐  │
│  │  File Finder │───▶│  Parsers     │───▶│  Data Transform      │  │
│  │  (Crawler)   │    │  - Ariba    │    │  - Schema Mapping   │  │
│  │              │    │  - StarFusion│    │  - Validation        │  │
│  └──────────────┘    └──────────────┘    └──────────────────────┘  │
│                                                      │               │
│                                                      ▼               │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │                    MSSQL Loader                              │   │
│  │  - Bulk insert (bcp/BULK INSERT)                           │   │
│  │  - Upsert logic for re-runs                                 │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                              │                                      │
│                              ▼                                      │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │                 Concordance Analyzer                          │   │
│  │  - Sample-matched Ariba + StarFusion merge                  │   │
│  │  - Shared/Unique/Discordant fusion classification           │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                              │                                      │
│                              ▼                                      │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │                   Workflow Orchestrator                       │   │
│  │  - Chain parsers → loader → concordance                      │   │
│  │  - CLI + Python API                                         │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
                               │
                               ▼
                    ┌─────────────────────┐
                    │   Output Formats    │
                    │  - MSSQL Tables    │
                    │  - TSV (test mode) │
                    │  - JSON (summary)   │
                    └─────────────────────┘
```

---

## Package Structure

```
starfusion_itd/
├── fusql/                      # Main Python package
│   ├── __init__.py
│   ├── pyproject.toml          # Package config
│   ├── cli.py                  # CLI entry point
│   ├── config.py               # Config management
│   │
│   ├── parsers/               # File parsers
│   │   ├── __init__.py
│   │   ├── base.py             # Base parser class
│   │   ├── ariba.py            # Ariba fusion parser
│   │   └── starfusion.py        # StarFusion/FusionInspector parser
│   │
│   ├── loaders/                # Data loaders
│   │   ├── __init__.py
│   │   ├── base.py             # Base loader class
│   │   ├── mssql.py            # MSSQL bulk loader
│   │   └── tsv.py              # TSV loader (test mode)
│   │
│   ├── concordance/             # Fusion concordance analysis
│   │   ├── __init__.py
│   │   ├── merger.py           # Merge Ariba + StarFusion
│   │   └── analyzer.py         # Concordance statistics
│   │
│   ├── discovery/               # File discovery
│   │   ├── __init__.py
│   │   ├── finder.py           # Recursive file finder
│   │   └── scanner.py          # Directory scanner with sample ID extraction
│   │
│   ├── workflow/               # Workflow orchestration
│   │   ├── __init__.py
│   │   ├── runner.py           # Main workflow runner
│   │   └── validators.py        # Input validation
│   │
│   └── utils/                  # Utilities
│       ├── __init__.py
│       ├── logging.py           # Logging setup
│       └── exceptions.py        # Custom exceptions
│
├── tests/                      # Test suite
│   ├── __init__.py
│   ├── test_ariba_parser.py
│   ├── test_starfusion_parser.py
│   ├── test_concordance.py
│   ├── test_discovery.py
│   └── fixtures/               # Test data
│
├── ITD/                        # Existing ITD detection code
│   └── ...
│
├── docs/                       # Documentation
│   ├── SPEC.md                 # This file
│   ├── ARCHITECTURE.md         # Detailed architecture
│   └── schema.md               # DB schema documentation
│
├── pyproject.toml              # Root project config
├── README.md
└── setup.py
```

---

## Database Schema

### Table: `ariba_fusions`
| Column | Type | Description |
|--------|------|-------------|
| `id` | BIGINT IDENTITY | Primary key |
| `run_id` | VARCHAR(100) | Sequencer run identifier |
| `sample_id` | VARCHAR(100) | Sample identifier |
| `fusion_id` | VARCHAR(200) | Unique fusion identifier |
| `gene1` | VARCHAR(50) | Left gene symbol |
| `gene2` | VARCHAR(50) | Right gene symbol |
| `exon1` | INT | Exon number in gene1 |
| `exon2` | INT | Exon number in gene2 |
| `splice_site_1` | VARCHAR(10) | Type: acceptor/donor/unknown |
| `splice_site_2` | VARCHAR(10) | Type: acceptor/donor/unknown |
| `reads` | INT | Supporting read count |
| `coding` | BIT | Is fusion in coding region |
| `genes_strand` | VARCHAR(10) | Strand orientation |
| `annotation` | VARCHAR(500) | Full annotation string |
| `input_file` | VARCHAR(500) | Source file path |
| `loaded_at` | DATETIME | Timestamp of load |

### Table: `starfusion_fusions`
| Column | Type | Description |
|--------|------|-------------|
| `id` | BIGINT IDENTITY | Primary key |
| `run_id` | VARCHAR(100) | Sequencer run identifier |
| `sample_id` | VARCHAR(100) | Sample identifier |
| `fusion_name` | VARCHAR(300) | Fusion display name |
| `gene1` | VARCHAR(50) | Left gene symbol |
| `gene2` | VARCHAR(50) | Right gene symbol |
| `exon1` | VARCHAR(20) | Exon structure (e.g., "E10") |
| `exon2` | VARCHAR(20) | Exon structure (e.g., "E11") |
| `splice_type` | VARCHAR(50) | SPLICE_SITE_TYPE from input |
| `left_breakpos` | BIGINT | Genomic position left |
| `right_breakpos` | BIGINT | Genomic position right |
| `junction_reads` | INT | Junction-spanning reads |
| `spanning_reads` | INT | Spanning reads |
| `FFPM` | FLOAT | Fusions per million FPMS |
| `CDS_left` | VARCHAR(50) | CDS annotation left |
| `CDS_right` | VARCHAR(50) | CDS annotation right |
| `input_file` | VARCHAR(500) | Source file path |
| `loaded_at` | DATETIME | Timestamp of load |

### Table: `fusion_concordance`
| Column | Type | Description |
|--------|------|-------------|
| `id` | BIGINT IDENTITY | Primary key |
| `run_id` | VARCHAR(100) | Sequencer run identifier |
| `sample_id` | VARCHAR(100) | Sample identifier |
| `fusion_id` | VARCHAR(300) | Normalized fusion ID |
| `gene1` | VARCHAR(50) | Left gene |
| `gene2` | VARCHAR(50) | Right gene |
| `ariba_found` | BIT | Found by Ariba |
| `starfusion_found` | BIT | Found by StarFusion |
| `concordance_status` | VARCHAR(20) | shared/unique_ariba/unique_starfusion |
| `ariba_exons` | VARCHAR(50) | Ariba exon info |
| `starfusion_exons` | VARCHAR(50) | StarFusion exon info |
| `splice_concordant` | BIT | Both at splice sites |
| `notes` | VARCHAR(500) | Additional notes |
| `analyzed_at` | DATETIME | Analysis timestamp |

---

## Function Specifications

### 1. File Discovery (`discovery/`)

#### `find_ariba_files(root_dir: Path, glob_pattern: str = "*.tsv", name_patterns: List[str] = None) -> List[Path]`
- Recursively find Ariba fusion report files
- **glob_pattern:** File extension filter (default: `*.tsv`)
- **name_patterns:** List of regex patterns for matching filenames
  - Default: `ariba.*report\.tsv$`, `.*ariba_report\.tsv$`
  - Custom: `["my_ariba_.*\\.tsv$", "run.*_ariba\\.tsv$"]`
- Return list of Path objects

#### `find_starfusion_files(root_dir: Path, glob_pattern: str = "*.tsv", name_patterns: List[str] = None) -> List[Path]`
- Recursively find StarFusion/FusionInspector output files
- **glob_pattern:** File extension filter (default: `*.tsv`)
- **name_patterns:** List of regex patterns for matching filenames
  - Default: `star.*fusion.*\.tsv$`, `fusion_inspector.*\.tsv$`
  - Custom: `["starfusion_.*\\.tsv$", "fi_.*\\.tsv$"]`
- Return list of Path objects

#### `extract_sample_id(file_path: Path, pattern: str) -> Optional[str]`
- Extract sample ID from file path using regex pattern
- Common patterns: `Sample_([^/]+)`, `Run_([^/]+)`, etc.
- Return sample ID string or None

#### `scan_directory(root_dir: Path, ariba_patterns: List[str] = None, starfusion_patterns: List[str] = None) -> Dict[str, Dict[str, Path]]`
- Scan directory for all Ariba and StarFusion files using custom patterns
- Group by sample ID
- Return structure: `{sample_id: {ariba: Path, starfusion: Path}}`

### 2. Parsers (`parsers/`)

#### Ariba Parser (`ariba.py`)

```python
class AribaParser:
    def parse(self, file_path: Path) -> List[AribaFusion]:
        """Parse Ariba fusion TSV file"""
        
    def extract_exon_number(self, annotation: str) -> Tuple[int, int]:
        """Extract exon numbers from annotation string"""
        
    def extract_splice_site(self, annotation: str) -> Tuple[str, str]:
        """Identify if fusion is at acceptor/donor site"""
```

**Ariba TSV Expected Columns:**
- `gene1`, `gene2` — gene symbols
- `gene1_strand`, `gene2_strand` — strand orientation
- `assembled_seqs` — alignment info with exon positions
- `coding` — yes/no if in coding region

#### StarFusion Parser (`starfusion.py`)

```python
class StarFusionParser:
    def parse(self, file_path: Path) -> List[StarFusionEntry]:
        """Parse StarFusion FusionInspector TSV"""
        
    def extract_exons(self, exon_info: str) -> Tuple[str, str]:
        """Extract exon structure from E10^E11 format"""
```

**StarFusion TSV Expected Columns:**
- `#FusionName` — fusion name (GENE1--GENE2)
- `JunctionRead`, `SpanningRead` — read counts
- `SpliceType` — SPLICE_SITE_TYPE
- `LeftBreakpoint`, `RightBreakpoint` — genomic positions
- `CDS_LEFT`, `CDS_RIGHT` — coding sequence annotation
- `FFPM` — fusions per million

### 3. Loaders (`loaders/`)

#### MSSQL Loader (`mssql.py`) — Production Mode

```python
class MSSQLLoader:
    def __init__(self, connection_string: str):
        """Initialize with MSSQL connection string"""
        
    def bulk_insert(self, table: str, rows: List[Dict]) -> int:
        """Append rows to table (partitioned by run_id + sample_id)"""
        
    def table_exists(self, table: str) -> bool:
        """Check if table exists"""
        
    def create_table(self, table: str, schema: Dict[str, str]) -> None:
        """Create table if not exists"""
```

**Design:**
- Append-only mode (no upsert)
- Partition key: `run_id` + `sample_id`
- Configurable table names per tool type (ariba_fusions, starfusion_fusions, fusion_concordance)

#### TSV Loader (`tsv.py`) — Test Mode

```python
class TSVLoader:
    def __init__(self, output_dir: Path):
        """Initialize with output directory"""
        
    def write(self, table_name: str, rows: List[Dict]) -> Path:
        """Write rows to TSV matching MSSQL schema"""
        # Output: output_dir/table_name.tsv
```

### 4. Concordance (`concordance/`)

#### Merger (`merger.py`)

```python
def merge_fusions(
    ariba_fusions: List[AribaFusion],
    starfusion_fusions: List[StarFusionEntry]
) -> List[FusionConcordance]:
    """
    Merge fusion calls from both callers.
    Normalize fusion IDs for matching.
    Classify as: shared, unique_ariba, unique_starfusion
    """
    
def normalize_fusion_id(gene1: str, gene2: str) -> str:
    """Normalize fusion ID: sort genes alphabetically, format as GENE1--GENE2"""
```

#### Analyzer (`analyzer.py`)

```python
def calculate_concordance_stats(merged: List[FusionConcordance]) -> Dict:
    """Calculate concordance statistics"""
    # Return: {total, shared, unique_ariba, unique_starfusion, 
    #          splice_concordance_rate, etc.}
    
def generate_summary_report(merged: List[FusionConcordance], stats: Dict) -> str:
    """Generate human-readable summary"""
```

### 5. Workflow (`workflow/`)

```python
class FusionWorkflow:
    def __init__(self, config: WorkflowConfig):
        self.loader = TSVLoader(...) if config.test_mode else MSSQLLoader(...)
        
    def run(
        self,
        input_dir: Path,
        sample_id: Optional[str] = None,
        output_dir: Optional[Path] = None
    ) -> WorkflowResult:
        """
        Main workflow:
        1. Discover files
        2. Parse Ariba (if found)
        3. Parse StarFusion (if found)
        4. Load to MSSQL/TSV
        5. Merge and analyze concordance
        6. Return results
        """
```

### 6. CLI (`cli.py`)

```bash
# Full workflow (DB mode) - auto-discover files with default patterns
fusql run /path/to/samples --run-id RUN001 --sample-id SAMPLE_001 \
  --mssql "mssql+pyodbc://user:pass@server/database?driver=ODBC+Driver+17+for+SQL+Server" \
  --table-ariba ariba_fusions \
  --table-starfusion starfusion_fusions

# With custom file patterns
fusql run /path/to/samples --run-id RUN001 --sample-id SAMPLE_001 \
  --ariba-patterns "my_ariba.*\\.tsv$" "ariba_.*report\\.tsv$" \
  --starfusion-patterns "starfusion.*\\.tsv$" "fi_.*\\.tsv$" \
  --mssql "Server=...;Database=..."

# Test mode (TSV output)
fusql run /path/to/samples --run-id RUN001 --sample-id SAMPLE_001 \
  --test-mode --output ./tsv_results

# Parse only
fusql parse-ariba /path/to/ariba.tsv --run-id RUN001 --sample-id SAMPLE_001 --output ariba_parsed.tsv
fusql parse-starfusion /path/to/starfusion.tsv --run-id RUN001 --sample-id SAMPLE_001 --output starfusion_parsed.tsv

# Concordance analysis
fusql merge --ariba ariba.tsv --starfusion starfusion.tsv --run-id RUN001 --sample-id SAMPLE_001 --output merged.tsv

# Discovery with patterns
fusql discover /path/to/samples \
  --ariba-patterns ".*ariba.*\\.tsv$" \
  --starfusion-patterns ".*fusion.*\\.tsv$"

# Discovery
fusql discover /path/to/samples --output discovered_files.json
```

---

## Configuration (`config.py`)

```python
@dataclass
class MSSQLConfig:
    server: str
    database: str
    username: str
    password: str
    driver: str = "ODBC Driver 17 for SQL Server"
    
@dataclass  
class WorkflowConfig:
    test_mode: bool = False
    output_dir: Optional[Path] = None
    mssql: Optional[MSSQLConfig] = None
    log_level: str = "INFO"
    sample_id_pattern: str = r"([A-Za-z0-9_-]+)"
    
# Load from config file or environment variables
# Config file: fusql.yaml or fusql.json
```

---

## Development Tasks (for sub-agents)

1. **Package Setup** — Create `fusql/` package structure, `pyproject.toml`, `__init__.py` files
2. **Ariba Parser** — Implement `parsers/ariba.py` with exon/splice site extraction
3. **StarFusion Parser** — Implement `parsers/starfusion.py` 
4. **MSSQL/TSV Loader** — Implement `loaders/mssql.py` and `loaders/tsv.py`
5. **File Discovery** — Implement `discovery/finder.py` and `scanner.py`
6. **Concordance Merger** — Implement `concordance/merger.py` and `analyzer.py`
7. **Workflow & CLI** — Implement `workflow/runner.py` and `cli.py`
8. **Tests** — Create test suite with fixtures

---

## Dependencies

```toml
[project]
dependencies = [
    "pandas>=2.0.0",
    "pyodbc>=5.0.0",          # MSSQL connection
    "psycopg2-binary>=2.9.0", # Backup DB option
    "click>=8.0.0",           # CLI framework
    "pyyaml>=6.0",            # Config files
    "pydantic>=2.0.0",        # Data validation
    "loguru>=0.7.0",          # Logging
    "rich>=13.0.0",           # CLI output formatting
]
```

---

## Next Steps

1. Create package structure
2. Implement core parsers first
3. Test with sample data
4. Add MSSQL loader
5. Implement concordance
6. CLI and workflow
