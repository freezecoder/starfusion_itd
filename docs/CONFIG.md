# FusionSQL Configuration Guide

FusionSQL uses a layered configuration system. Settings are resolved in this priority order (highest to lowest):

1. **CLI arguments** (e.g., `--config`, `--log-level`)
2. **Environment variables** (`FUSQL_*`)
3. **Config file** (`./fusql.yaml` or `~/.fusql.yaml`)
4. **Built-in defaults**

## Config File Format

Create a `fusql.yaml` file in your home directory or project root:

```yaml
mssql:
  connection_string: "mssql+pyodbc://user:pass@server/db?driver=ODBC+Driver+17+for+SQL+Server"
  driver: "ODBC Driver 17 for SQL Server"
table_ariba: ariba_fusions
table_starfusion: starfusion_fusions
table_concordance: fusion_concordance
ariba_patterns:
  - "^ariba_report\\.tsv$"
  - "ariba.*report\\.tsv$"
starfusion_patterns:
  - "^star-fusion\\.fusion_predictions.*\\.tsv$"
  - "^fusion_inspector.*\\.tsv$"
log_level: INFO
```

## Configuration Options

### `mssql.connection_string`

Full MSSQL connection string in SQLAlchemy URL format.

Example:
```
mssql+pyodbc://user:pass@hostname:1433/database?driver=ODBC+Driver+17+for+SQL+Server
```

### `mssql.driver`

ODBC driver name. Default: `ODBC Driver 17 for SQL Server`

### `table_ariba`

MSSQL table name for Ariba fusion data. Default: `ariba_fusions`

### `table_starfusion`

MSSQL table name for STARFusion/FusionInspector data. Default: `starfusion_fusions`

### `table_concordance`

MSSQL table name for concordance analysis results. Default: `fusion_concordance`

### `ariba_patterns`

List of regex patterns for identifying Ariba report files. Files matching any pattern are treated as Ariba output.

Default patterns:
- `^ariba_report\.tsv$`
- `ariba.*report\.tsv$`

### `starfusion_patterns`

List of regex patterns for identifying STARFusion/FusionInspector files.

Default patterns:
- `^star-fusion\.fusion_predictions.*\.tsv$`
- `^fusion_inspector.*\.tsv$`

### `log_level`

Logging verbosity. Options: `DEBUG`, `INFO`, `WARNING`, `ERROR`. Default: `INFO`

## Environment Variables

All configuration options can be set via environment variables prefixed with `FUSQL_`:

| Variable | Description | Example |
|----------|-------------|---------|
| `FUSQL_MSSQL` | MSSQL connection string | `mssql+pyodbc://...` |
| `FUSQL_DRIVER` | ODBC driver name | `ODBC Driver 17 for SQL Server` |
| `FUSQL_TABLE_ARIBA` | Ariba table name | `ariba_fusions` |
| `FUSQL_TABLE_STARFUSION` | STARFusion table name | `starfusion_fusions` |
| `FUSQL_TABLE_CONCORDANCE` | Concordance table name | `fusion_concordance` |
| `FUSQL_ARIBA_PATTERNS` | JSON list of Ariba patterns | `["^ariba.*\\.tsv$"]` |
| `FUSQL_STARFUSION_PATTERNS` | JSON list of STARFusion patterns | `["^starfusion.*\\.tsv$"]` |
| `FUSQL_LOG_LEVEL` | Logging level | `DEBUG` |

**Note:** For list-type variables (`*_PATTERNS`), pass a JSON-encoded list.

Example:
```bash
export FUSQL_MSSQL="mssql+pyodbc://user:pass@server/db?driver=ODBC+Driver+17+for+SQL+Server"
export FUSQL_LOG_LEVEL=DEBUG
export FUSQL_ARIBA_PATTERNS='["^my_ariba.*\\.tsv$", "^ariba_report\\.tsv$"]'
```

## CLI Usage

### Specify Config File

```bash
fusql --config /path/to/fusql.yaml run /data/samples
```

### Override Log Level

```bash
fusql --log-level DEBUG run /data/samples
```

### Show Current Configuration

```bash
fusql --show-config
```

This displays the fully-resolved config (with env vars and CLI overrides applied).

## Creating a Config File Programmatically

```python
from fusql.config import FusionSQLConfig, save_config
from pathlib import Path

config = FusionSQLConfig()
config.mssql.connection_string = "mssql+pyodbc://user:pass@server/db"
config.table_ariba = "my_ariba_fusions"
config.log_level = "DEBUG"

save_config(config, Path("~/.fusql.yaml"))
```

## Loading Config in Python

```python
from fusql.config import load_config

# Auto-detect config file (~/.fusql.yaml or ./fusql.yaml)
config = load_config()

# Explicit config file
config = load_config(Path("/path/to/fusql.yaml"))

# Access settings
print(config.mssql.connection_string)
print(config.table_ariba)
print(config.ariba_patterns)
```

## Database Behavior

### Append-Only Mode

**FusionSQL always operates in append-only mode for MSSQL.** This means:

- ✅ **INSERT only** — never UPDATE, DELETE, or MERGE
- ✅ **Idempotent runs** — re-processing same sample adds new rows (does not overwrite)
- ✅ **Data preserved** — existing data is never modified or deleted
- ✅ **Partitioned by run_id + sample_id** — each run's data is independently identifiable

This ensures safe operation in production environments where data integrity is critical.

### Example: Re-running a Sample

If you re-run the same sample:

```bash
# First run
fusql run /data --run-id RUN1 --sample-id SAMPLE1 --mssql "..."

# Second run (same sample)
fusql run /data --run-id RUN1 --sample-id SAMPLE1 --mssql "..."
```

Both runs insert independently. You will have duplicate rows with different `id` values. Query by `run_id` to distinguish them.

## Output Formats

### TSV Output (Test Mode)

Default format. Tab-separated values matching the MSSQL schema exactly.

```bash
fusql parse-ariba input.tsv --run-id RUN1 --sample-id S1 -o output.tsv
```

### Excel Output

For Excel (.xlsx) output, specify the file extension:

```bash
fusql parse-ariba input.tsv --run-id RUN1 --sample-id S1 -o output.xlsx
fusql parse-starfusion input.tsv --run-id RUN1 --sample-id S1 -o output.xlsx
```

Excel files include:
- Formatted headers (blue background, white bold text)
- Auto-adjusted column widths
- All data preserved in a spreadsheet-friendly format

## Output Templates

You can customize which fields appear in output files using `output_templates` in the config:

```yaml
output_templates:
  ariba:
    fields:
      - run_id
      - sample_id
      - gene1
      - gene2
      - exon1
      - exon2
      - reads
  starfusion:
    fields:
      - run_id
      - sample_id
      - gene1
      - gene2
      - junction_reads
      - spanning_reads
  concordance:
    fields:
      - run_id
      - sample_id
      - fusion_id
      - concordance_status
```

When `fields` is empty or omitted, all fields are included.

Use `--all-fields` flag to override and include all fields regardless of template.

## Loading Credentials from Python

For sensitive credentials, you can load the MSSQL connection string from a Python file:

```bash
# creds.py
def get_url():
    return "mssql+pyodbc://user:pass@server/db?driver=ODBC+Driver+17+for+SQL+Server"
```

```bash
fusql testdb --creds /path/to/creds.py
fusql run /data --run-id RUN1 --sample-id S1 --creds /path/to/creds.py
```

Supports:
- Function: `get_url()`, `get_connection_string()`
- Variable: `URL`, `url`, `CONNECTION_STRING`
