# Image Reference

This document serves as a reference for all visual assets in the FusionSQL project.

---

## System Architecture Diagram

**File:** `docs/fusql_architecture.png`  
**Size:** ~500 KB  
**Resolution:** 2K  
**Generated:** 2026-03-25

### Description

High-level system architecture showing the complete FusionSQL data pipeline:

```
┌─────────────────────────────────────────────────────────────────────┐
│                        FusionSQL Package                             │
├─────────────────────────────────────────────────────────────────────┤
│                                                                     │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────────────┐  │
│  │  File Finder │───▶│  Parsers     │───▶│  Data Transform      │  │
│  │  (Crawler)   │    │  - Arriba  │    │  - Schema Mapping   │  │
│  │              │    │  - STARFusion│   │  - Validation        │  │
│  └──────────────┘    └──────────────┘    └──────────────────────┘  │
│                                                      │               │
│                                                      ▼               │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │                    MSSQL Loader                              │   │
│  │  - Bulk insert (bcp/BULK INSERT)                           │   │
│  │  - Append mode (partitioned by run_id + sample_id)           │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                              │                                      │
│                              ▼                                      │
│  ┌──────────────────────────────────────────────────────────────┐   │
│  │                 Concordance Analyzer                          │   │
│  │  - Sample-matched Arriba + STARFusion merge                  │   │
│  │  - Shared/Unique/Discordant fusion classification            │   │
│  └──────────────────────────────────────────────────────────────┘   │
│                                                                     │
└─────────────────────────────────────────────────────────────────────┘
```

### Usage

```markdown
![System Architecture](docs/fusql_architecture.png)
```

---

## Database Schema Diagram

**File:** `docs/fusql_schema.png`  
**Size:** ~640 KB  
**Resolution:** 2K  
**Generated:** 2026-03-25

### Description

Microsoft SQL Server database schema showing the three main tables and their relationships:

### Tables

#### 1. `ariba_fusions`
Stores Arriba fusion caller results.

| Column | Type | Description |
|--------|------|-------------|
| id | BIGINT (PK) | Primary key |
| run_id | VARCHAR(100) | Sequencer run identifier |
| sample_id | VARCHAR(100) | Sample identifier |
| fusion_id | VARCHAR(200) | Unique fusion identifier |
| gene1 | VARCHAR(50) | Left gene symbol |
| gene2 | VARCHAR(50) | Right gene symbol |
| exon1 | INT | Exon number in gene1 |
| exon2 | INT | Exon number in gene2 |
| splice_site_1 | VARCHAR(10) | Type: acceptor/donor/splice-site |
| splice_site_2 | VARCHAR(10) | Type: acceptor/donor/splice-site |
| reads | INT | Supporting read count |
| coding | BIT | Is fusion in coding region |
| genes_strand | VARCHAR(10) | Strand orientation |
| reading_frame | VARCHAR(20) | in-frame, out-of-frame, unknown |
| confidence | VARCHAR(20) | high, medium, low |
| fusion_type | VARCHAR(50) | translocation, deletion, etc. |
| input_file | VARCHAR(500) | Source file path |
| loaded_at | DATETIME | Timestamp of load |

#### 2. `starfusion_fusions`
Stores STARFusion + FusionInspector results.

| Column | Type | Description |
|--------|------|-------------|
| id | BIGINT (PK) | Primary key |
| run_id | VARCHAR(100) | Sequencer run identifier |
| sample_id | VARCHAR(100) | Sample identifier |
| fusion_name | VARCHAR(300) | Fusion display name |
| gene1 | VARCHAR(50) | Left gene symbol |
| gene2 | VARCHAR(50) | Right gene symbol |
| exon1 | VARCHAR(20) | Exon info (e.g., E1-923) |
| exon2 | VARCHAR(20) | Exon info (e.g., E1698-3270) |
| splice_type | VARCHAR(50) | INCL_NON_REF_SPLICE, ONLY_REF_SPLICE |
| left_breakpos | BIGINT | Left genomic breakpoint |
| right_breakpos | BIGINT | Right genomic breakpoint |
| junction_reads | INT | Junction-spanning reads |
| spanning_reads | INT | Spanning fragments |
| ffpm | FLOAT | Fusions per million FPMS |
| cds_left | VARCHAR(50) | CDS annotation for left gene |
| cds_right | VARCHAR(50) | CDS annotation for right gene |
| prot_fusion_type | VARCHAR(50) | INFRAME, FRAMESHIFT |
| large_anchor_support | VARCHAR(10) | YES/NO |
| annots | VARCHAR(1000) | Annotation categories |
| input_file | VARCHAR(500) | Source file path |
| loaded_at | DATETIME | Timestamp of load |

#### 3. `fusion_concordance`
Stores merged concordance analysis results.

| Column | Type | Description |
|--------|------|-------------|
| id | BIGINT (PK) | Primary key |
| run_id | VARCHAR(100) | Sequencer run identifier |
| sample_id | VARCHAR(100) | Sample identifier |
| fusion_id | VARCHAR(300) | Normalized fusion ID |
| gene1 | VARCHAR(50) | Left gene |
| gene2 | VARCHAR(50) | Right gene |
| ariba_found | BIT | Found by Arriba (0/1) |
| starfusion_found | BIT | Found by STARFusion (0/1) |
| concordance_status | VARCHAR(30) | shared, unique_ariba, unique_starfusion |
| splice_concordant | BIT | Both at splice sites (0/1) |
| ariba_exons | VARCHAR(50) | Arriba exon info |
| starfusion_exons | VARCHAR(50) | STARFusion exon info |
| notes | VARCHAR(500) | Additional notes |
| analyzed_at | DATETIME | Analysis timestamp |

### Relationships

```
ariba_fusions ─────┐
                   ├──▶ fusion_concordance ◀──┤
starfusion_fusions ┘
```

Both fusion tables feed into the concordance table for merged analysis.

### Usage

```markdown
![Database Schema](docs/fusql_schema.png)
```

---

## Generating New Images

Both images were generated using the Gemini image generation skill:

```bash
cd /data/.openclaw/workspace/skills/gemini-image
uv run scripts/generate_image.py \
  --prompt "your prompt here" \
  --filename output.png \
  --resolution 2K
```

### Image Specifications

- **Format:** PNG
- **Resolution:** 2K (2048px width)
- **Aspect Ratio:** Auto (prompt-guided)
- **Style:** Professional, clean diagrams, white background

---

## File Locations

All images are stored in the `docs/` directory:

```
starfusion_itd/
└── docs/
    ├── fusql_architecture.png  # System architecture
    ├── fusql_schema.png        # Database schema
    ├── PRESENTATION.md         # Slide deck
    ├── IMAGES.md              # This reference
    └── SPEC.md                # Full specification
```
