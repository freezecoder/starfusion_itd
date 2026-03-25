# FusionSQL CLI Usage Guide

Complete reference for the `fusql` command-line interface.

---

## Installation

```bash
# Create conda environment
conda create -n fusql python=3.13 && conda activate fusql

# Install dependencies
pip install pandas sqlalchemy pyodbc click pyyaml pydantic loguru rich

# Install fusql in editable mode
pip install -e /path/to/starfusion_itd/fusql
```

---

## Quick Start

### 1. Discover Files in a Directory

```bash
fusql discover /path/to/run_data
```

Output:
```json
{
  "Z3AT9/Z3AT9_IonCode_0125": {
    "run_id": "Z3AT9",
    "sample_id": "Z3AT9_IonCode_0125",
    "ariba": "/path/to/Z3AT9/pipelineout/Z3AT9_IonCode_0125/ariba_report.tsv",
    "starfusion": "/path/to/Z3AT9/pipelineout/Z3AT9_IonCode_0125/star-fusion.fusion_predictions.abridged.tsv"
  },
  "Z3AT9/Z3AT9_IonCode_0126": {
    "run_id": "Z3AT9",
    "sample_id": "Z3AT9_IonCode_0126",
    "ariba": "/path/to/Z3AT9/pipelineout/Z3AT9_IonCode_0126/ariba_report.tsv",
    "starfusion": "/path/to/Z3AT9/pipelineout/Z3AT9_IonCode_0126/star-fusion.fusion_predictions.abridged.tsv"
  }
}
```

### 2. Run Full Workflow (Test Mode → TSV Output)

```bash
fusql run /path/to/run_data \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --test-mode \
  --output ./results
```

Output files:
```
results/
├── ariba_fusions.tsv
├── starfusion_fusions.tsv
└── fusion_concordance.tsv
```

### 3. Run Full Workflow (Database Mode → MSSQL)

```bash
fusql run /path/to/run_data \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --mssql "mssql+pyodbc://user:pass@server/database?driver=ODBC+Driver+17+for+SQL+Server" \
  --table-ariba ariba_fusions \
  --table-starfusion starfusion_fusions
```

---

## Directory Structure

FusionSQL expects this directory structure:

```
{run_id}/
└── pipelineout/
    └── {sample_id}/
        ├── arriba_report.tsv          # Arriba output
        └── star-fusion.fusion_predictions.abridged.tsv  # StarFusion output
```

Example:
```
Z3AT9/
└── pipelineout/
    ├── Z3AT9_IonCode_0125/
    │   ├── arriba_report.tsv
    │   └── star-fusion.fusion_predictions.abridged.tsv
    └── Z3AT9_IonCode_0126/
        ├── arriba_report.tsv
        └── star-fusion.fusion_predictions.abridged.tsv
```

---

## Commands

### `fusql discover`

Scan a directory and discover all fusion files.

```bash
fusql discover [OPTIONS] ROOT_DIR
```

**Options:**
| Option | Description | Default |
|--------|-------------|---------|
| `--ariba-patterns` | Custom regex patterns for Ariba files (space-separated) | `ariba_report.tsv` |
| `--starfusion-patterns` | Custom regex patterns for StarFusion files | `star-fusion.fusion_predictions.*.tsv` |
| `--output` | Output file (JSON) | stdout |

**Examples:**

```bash
# Basic discovery
fusql discover /data/run_folder

# Save to file
fusql discover /data/run_folder --output discovered.json

# With custom patterns
fusql discover /data/run_folder \
  --ariba-patterns "custom_ariba.*\.tsv$" \
  --starfusion-patterns "my_fusion.*\.tsv$"
```

---

### `fusql run`

Run the full workflow: discover → parse → load → analyze concordance.

```bash
fusql run [OPTIONS] ROOT_DIR
```

**Required Options:**
| Option | Description |
|--------|-------------|
| `--run-id` | Sequencer run identifier (e.g., `Z3AT9`) |
| `--sample-id` | Sample identifier (e.g., `Z3AT9_IonCode_0125`) |

**Mode Options:**
| Option | Description | Default |
|--------|-------------|---------|
| `--test-mode` | Output TSV files instead of MSSQL | `False` |
| `--mssql` | MSSQL connection string (required if not `--test-mode`) | - |
| `--output` | Output directory for TSV files (required if `--test-mode`) | `./output` |

**Table Options:**
| Option | Description | Default |
|--------|-------------|---------|
| `--table-ariba` | Table name for Arriba fusions | `ariba_fusions` |
| `--table-starfusion` | Table name for StarFusion fusions | `starfusion_fusions` |
| `--table-concordance` | Table name for concordance results | `fusion_concordance` |

**File Pattern Options:**
| Option | Description | Default |
|--------|-------------|---------|
| `--ariba-patterns` | Custom regex patterns for Ariba files | see patterns |
| `--starfusion-patterns` | Custom regex patterns for StarFusion files | see patterns |

**Examples:**

```bash
# Test mode with TSV output
fusql run /data/Z3AT9 \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --test-mode \
  --output ./results

# Database mode
fusql run /data/Z3AT9 \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --mssql "mssql+pyodbc://user:pass@mssql.example.com/mydb?driver=ODBC+Driver+17+for+SQL+Server" \
  --table-ariba ariba_fusions \
  --table-starfusion starfusion_fusions

# With custom file patterns
fusql run /data/Z3AT9 \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --ariba-patterns "my_ariba.*\.tsv$" \
  --starfusion-patterns "starfusion.*\.tsv$" \
  --test-mode \
  --output ./results

# Custom MSSQL tables
fusql run /data/Z3AT9 \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --mssql "mssql+pyodbc://user:pass@server/db" \
  --table-ariba fusion_calls_ariba \
  --table-starfusion fusion_calls_starfusion \
  --table-concordance fusion_comparison
```

---

### `fusql parse-ariba`

Parse an Arriba TSV file without loading to DB.

```bash
fusql parse-ariba [OPTIONS] INPUT_FILE
```

**Options:**
| Option | Description | Required |
|--------|-------------|----------|
| `--run-id` | Run identifier | Yes |
| `--sample-id` | Sample identifier | Yes |
| `--output` | Output TSV file | Yes |

**Example:**
```bash
fusql parse-ariba /data/Z3AT9/pipelineout/Z3AT9_IonCode_0125/ariba_report.tsv \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --output ./ariba_parsed.tsv
```

---

### `fusql parse-starfusion`

Parse a StarFusion/FusionInspector TSV file without loading to DB.

```bash
fusql parse-starfusion [OPTIONS] INPUT_FILE
```

**Options:**
| Option | Description | Required |
|--------|-------------|----------|
| `--run-id` | Run identifier | Yes |
| `--sample-id` | Sample identifier | Yes |
| `--output` | Output TSV file | Yes |

**Example:**
```bash
fusql parse-starfusion /data/Z3AT9/pipelineout/Z3AT9_IonCode_0125/star-fusion.fusion_predictions.abridged.tsv \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --output ./starfusion_parsed.tsv
```

---

### `fusql merge`

Analyze concordance between Arriba and StarFusion results.

```bash
fusql merge [OPTIONS]
```

**Options:**
| Option | Description | Required |
|--------|-------------|----------|
| `--ariba` | Arriba parsed TSV file | Yes |
| `--starfusion` | StarFusion parsed TSV file | Yes |
| `--run-id` | Run identifier | Yes |
| `--sample-id` | Sample identifier | Yes |
| `--output` | Output TSV file | Yes |

**Example:**
```bash
fusql merge \
  --ariba ./ariba_parsed.tsv \
  --starfusion ./starfusion_parsed.tsv \
  --run-id Z3AT9 \
  --sample-id Z3AT9_IonCode_0125 \
  --output ./concordance.tsv
```

---

## Output Formats

### TSV Output (Test Mode)

Files are tab-separated with headers matching the database schema:

**`ariba_fusions.tsv`**
```tsv
id	run_id	sample_id	fusion_id	gene1	gene2	exon1	exon2	splice_site_1	splice_site_2	reads	coding	genes_strand	reading_frame	confidence	fusion_type	input_file	loaded_at
1	Z3AT9	Z3AT9_IonCode_0125	BCR--ABL1	BCR	ABL1	1	2	donor	acceptor	45	true	+/+	in-frame	high	translocation	/data/.../ariba_report.tsv	2026-03-25T10:00:00
```

**`starfusion_fusions.tsv`**
```tsv
id	run_id	sample_id	fusion_name	gene1	gene2	exon1	exon2	splice_type	left_breakpos	right_breakpos	junction_reads	spanning_reads	ffpm	cds_left	cds_right	prot_fusion_type	large_anchor_support	annots	input_file	loaded_at
```

**`fusion_concordance.tsv`**
```tsv
id	run_id	sample_id	fusion_id	gene1	gene2	ariba_found	starfusion_found	concordance_status	splice_concordant	ariba_exons	starfusion_exons	notes	analyzed_at
```

---

## Concordance Status Values

| Status | Meaning |
|--------|---------|
| `shared` | Found by both Arriba AND StarFusion |
| `unique_ariba` | Only Arriba detected this fusion |
| `unique_starfusion` | Only StarFusion detected this fusion |

---

## MSSQL Connection String Format

```
mssql+pyodbc://user:password@server:port/database?driver=ODBC+Driver+17+for+SQL+Server
```

**Example:**
```
mssql+pyodbc://sa:MyPassword@localhost:1433/FusionDB?driver=ODBC+Driver+17+for+SQL+Server
```

**Connection String Components:**
| Component | Description | Example |
|-----------|-------------|---------|
| `user` | Username | `sa` |
| `password` | Password | `MyPassword` |
| `server` | Server hostname/IP | `localhost` or `192.168.1.100` |
| `port` | Port (optional, default 1433) | `1433` |
| `database` | Database name | `FusionDB` |
| `driver` | ODBC driver name | `ODBC+Driver+17+for+SQL+Server` |

---

## Batch Processing Multiple Samples

To process multiple samples from the same run:

```bash
# Discover all samples first
fusql discover /data/Z3AT9 --output samples.json

# Process each sample (loop)
for sample in $(cat samples.json | jq -r 'keys[]'); do
  fusql run /data/Z3AT9 \
    --run-id Z3AT9 \
    --sample-id "$sample" \
    --mssql "mssql+pyodbc://user:pass@server/db?driver=ODBC+Driver+17+for+SQL+Server" \
    --table-ariba ariba_fusions \
    --table-starfusion starfusion_fusions
done
```

Or use the discovered JSON directly:

```bash
# Get sample list
SAMPLES=$(cat samples.json | jq -r 'keys[]')

# Process in parallel (if multiple samples)
echo "$SAMPLES" | xargs -P 4 -I {} fusql run /data/Z3AT9 \
  --run-id Z3AT9 \
  --sample-id {} \
  --mssql "mssql+pyodbc://user:pass@server/db?driver=ODBC+Driver+17+for+SQL+Server"
```

---

## `fusql sync`

Incremental sync — process only new samples not already in DB.

```bash
fusql sync [OPTIONS] WATCH_DIR
```

**Description:**
1. Query database for existing `(run_id, sample_id)` pairs
2. Scan watch directory for new samples
3. Process only samples not already loaded
4. Load to MSSQL

**Options:**
| Option | Description | Required |
|--------|-------------|----------|
| `--watch-table` | DB table to check for existing samples | Yes |
| `--mssql` | MSSQL connection string | Yes (if not using config) |
| `--schedule` | Schedule: `daily`, `hourly`, or custom cron | No |
| `--cron` | Custom cron expression (e.g., `0 2 * * *`) | No |

**Examples:**

```bash
# Single sync run
fusql sync /data/runs \
  --watch-table ariba_fusions \
  --mssql "mssql+pyodbc://user:pass@server/db?driver=ODBC+Driver+17+for+SQL+Server"

# With config file
fusql sync /data/runs --watch-table ariba_fusions

# Daily schedule (adds to cron)
fusql sync /data/runs --watch-table ariba_fusions --schedule daily

# Custom cron (every day at 2 AM)
fusql sync /data/runs --watch-table ariba_fusions --cron "0 2 * * *"

# Hourly sync
fusql sync /data/runs --watch-table ariba_fusions --schedule hourly
```

**Output:**
```
[INFO] Found 10 samples in database
[INFO] Discovered 15 samples in /data/runs
[INFO] 5 new samples to process
[INFO] Processing: Z3AT9/Z3AT9_IonCode_0125
[INFO] Processing: Z3AT9/Z3AT9_IonCode_0126
...
[INFO] Loaded 5 new samples in 45.2s
```

---

## Help

For full help:

```bash
fusql --help
fusql <command> --help
```

---

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Error (file not found, parse error, DB error) |
| 2 | Invalid arguments |
