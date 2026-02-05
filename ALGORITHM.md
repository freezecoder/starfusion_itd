# Intra-Gene Fusion Detection Algorithm Documentation

## Overview

This algorithm identifies intra-gene fusions (within the same gene) between non-continuous exons from chimeric junction data. It processes chimeric read alignments to detect fusion breakpoints and determines which exons are involved in each fusion event.

## Input Data Format

### Chimeric Junction File (`Chimeric.out.junction.tsv`)

The algorithm expects a tab-separated file with the following columns:

- `chrom1`, `coord1`, `strand1`: Chromosome, coordinate, and strand of the first breakpoint
- `chrom2`, `coord2`, `strand2`: Chromosome, coordinate, and strand of the second breakpoint
- `breakpoint1_pos`, `CIGAR1`: Precise breakpoint position and CIGAR string for breakpoint 1
- `breakpoint2_pos`, `CIGAR2`: Precise breakpoint position and CIGAR string for breakpoint 2
- `read_name`: Identifier for the supporting read
- `junction_type`: Type of junction (0, 1, or 2)

### Exon Annotation File (`EGFR_exons_hg19.tsv`)

A tab-separated file containing exon coordinates:
- `chromosome`: Chromosome name
- `start`, `end`: Exon boundaries
- `gene`: Gene name
- `exon`: Exon number
- `strand`: Strand orientation (+ or -)

## Algorithm Flow

### Step 1: Data Loading

1. Load exon annotations from the exon file
2. Load chimeric junction data from the junction file
3. Extract gene name from exon annotations (assumes all exons belong to the same gene)

### Step 2: Junction Processing

For each chimeric junction in the input file:

#### 2.1 Extract Breakpoint Information

The algorithm extracts breakpoint positions from two sources, in order of preference:

1. **Primary source**: `breakpoint1_pos` and `breakpoint2_pos` columns
   - These represent the **precise breakpoint positions** determined by the aligner
   - These positions are typically derived from **soft-clipped alignments** where:
     - Soft-clips (S) in the CIGAR string indicate unaligned sequence at read ends
     - The breakpoint position corresponds to where the soft-clip occurs in the reference
     - For example, a CIGAR like `68M75S` indicates 68 matched bases followed by 75 soft-clipped bases, with the breakpoint at position 68

2. **Fallback source**: `coord1` and `coord2` columns
   - Used only if `breakpoint1_pos` or `breakpoint2_pos` are missing/NaN
   - These are the general junction coordinates

#### 2.2 Extract CIGAR Information

- `CIGAR1` and `CIGAR2` are extracted from the input
- CIGAR strings encode the alignment structure:
  - **M**: Match/Mismatch
  - **S**: Soft clip (unaligned sequence at read ends - indicates breakpoint)
  - **D**: Deletion
  - **I**: Insertion
  - **N**: Skipped region (intron)

**Key Point**: The CIGAR strings contain soft-clip information (S) that indicates where reads are split, which helps identify the precise breakpoint location. However, the algorithm **does not parse the CIGAR strings** to extract breakpoint positions - it relies on the pre-computed `breakpoint1_pos` and `breakpoint2_pos` values provided in the input file.

#### 2.3 Strandedness Enforcement

**YES, the algorithm enforces strandedness** in the following ways:

1. **Exon Matching**: When determining which exon a breakpoint falls into, the algorithm filters exons by:
   - Matching chromosome (`chromosome == chrom1` or `chrom2`)
   - **Matching strand** (`strand == strand1` or `strand2`)
   
   This means a breakpoint on the `+` strand will only match exons annotated on the `+` strand, and vice versa.

2. **Strand Consistency**: The algorithm uses `strand1` and `strand2` from each junction to:
   - Match breakpoints to exons with the same strand orientation
   - Report strand information in the output (`left_strand`, `right_strand`)

**Important**: The algorithm does NOT require both breakpoints to be on the same strand. Each breakpoint is matched independently to exons based on its own strand orientation.

#### 2.4 Filtering Criteria

The algorithm applies several filters:

1. **Intra-chromosomal**: Both breakpoints must be on the same chromosome (`chrom1 == chrom2`)
2. **Exon Assignment**: Both breakpoints must map to known exons (within or near exon boundaries)
3. **Different Exons**: The two breakpoints must map to different exons (same-exon events are excluded)
4. **Same Gene**: Both exons must belong to the same gene (enforced by using a single gene's exon file)

### Step 3: Exon Assignment

For each breakpoint position, the algorithm determines which exon it belongs to using `find_exon_for_position()`:

#### 3.1 Matching Strategy (in order of priority):

1. **Within Exon Boundaries**:
   - If `exon_start <= position <= exon_end`: Assign to that exon
   - Distance to exon = 0

2. **Near Exon Boundaries** (within 50bp window):
   - If position is within 50bp of exon start or end: Assign to that exon
   - Distance = absolute difference from boundary

3. **Within Gene Span** (with 10kb buffer):
   - If position is within the gene's overall span (first exon start - 10kb to last exon end + 10kb):
   - Find the closest exon
   - Assign if distance ≤ 50kb from exon boundary

4. **No Match**: Return `None` if no exon is found

**Strand Matching**: Throughout this process, only exons matching the breakpoint's strand orientation are considered.

### Step 4: Fusion Classification

For each valid junction:

1. **Calculate Exon Distance**: `abs(exon1 - exon2)`
2. **Determine Continuity**:
   - **Continuous**: Exons are adjacent (`abs(exon1 - exon2) == 1`)
     - Example: exon2-exon3
     - These are **penalized** (`penalty_score = 1`)
   - **Non-continuous**: Exons are not adjacent (`abs(exon1 - exon2) > 1`)
     - Example: exon21-exon18
     - These are **preferred** (`penalty_score = 0`)

### Step 5: Aggregation and Summary

Junctions are grouped by fusion type (exon1-exon2 combination), and for each group:

#### 5.1 Aggregated Breakpoint Position

- **Method**: Most common position (mode) among all supporting reads
- **Format**: `chr:pos` (e.g., `chr7:55268981`)
- Calculated separately for left and right breakpoints

#### 5.2 Position-CIGAR Pairs

- Each breakpoint position is paired with its CIGAR string
- Format: `position:CIGAR,position:CIGAR,...`
- Example: `55268981:668M75S,55268971:JEM75S`
- Duplicates are removed while preserving order

#### 5.3 Read IDs

- All read IDs supporting each breakpoint are collected
- Format: Comma-separated list

#### 5.4 Output Columns

The final output includes:

1. `gene`: Gene name
2. `supporting_reads`: Count of reads supporting this fusion
3. `left_exon_anno`, `right_exon_anno`: Exon numbers involved
4. `left_breakpoint_aggregated`, `right_breakpoint_aggregated`: Consensus breakpoint positions (chr:pos)
5. `left_strand`, `right_strand`: Strand orientations
6. `breakpoint1_positions`, `breakpoint2_positions`: All positions with CIGARs
7. `breakpoint1_read_ids`, `breakpoint2_read_ids`: Supporting read IDs

## Key Algorithm Characteristics

### Breakpoint Determination

- **Source**: Breakpoint positions come from the input file's `breakpoint1_pos` and `breakpoint2_pos` columns
- **Soft-clip Information**: While CIGAR strings contain soft-clip (S) information that indicates split reads, the algorithm does not parse CIGAR strings to extract breakpoints
- **Aggregation**: The final aggregated breakpoint uses the most common position (mode) across all supporting reads

### Strandedness

- **Enforced**: Yes, strictly enforced
- **Matching**: Each breakpoint is matched to exons based on its strand orientation
- **Independence**: Left and right breakpoints can have different strand orientations
- **Output**: Strand information is reported for both breakpoints

### Continuous vs Non-continuous Exons

- **Non-continuous** (preferred): Exon distance > 1
  - Example: exon21-exon18 (distance = 3)
  - `penalty_score = 0`
  
- **Continuous** (penalized): Adjacent exons (distance = 1)
  - Example: exon2-exon3
  - `penalty_score = 1`
  - Still reported but marked as penalized

### Sorting

Results are sorted by:
1. `penalty_score` (ascending): Non-continuous first
2. `exon_distance` (descending): Larger exon distances first

## Limitations and Notes

1. **CIGAR Parsing**: The algorithm does not parse CIGAR strings to extract breakpoint positions. It relies on pre-computed breakpoint positions in the input file.

2. **Single Gene**: The algorithm is designed for a single gene (all exons must belong to the same gene).

3. **Distance Thresholds**:
   - Exon boundary window: 50bp
   - Gene span buffer: 10kb
   - Maximum intron distance: 50kb

4. **Breakpoint Precision**: The precision of breakpoint positions depends on the quality of the input chimeric junction file and how breakpoints were originally called by the aligner.

5. **Soft-clip Interpretation**: While CIGAR strings are stored and reported, the algorithm does not use them to refine breakpoint positions. The soft-clip information (S in CIGAR) is preserved in the output for manual inspection.

## Example

For a junction with:
- `breakpoint1_pos = 55268981`, `CIGAR1 = "668M75S"` (68 matched, 75 soft-clipped)
- `breakpoint2_pos = 55241614`, `CIGAR2 = "68S7SM"` (68 soft-clipped, 7 matched)
- `strand1 = "+"`, `strand2 = "+"`

The algorithm:
1. Uses positions 55268981 and 55241614 directly
2. Matches to exons on the `+` strand only
3. Assigns to exon21 and exon18 respectively
4. Classifies as non-continuous (distance = 3)
5. Aggregates with other reads supporting the same fusion
6. Outputs: `chr7:55268981` and `chr7:55241614` as aggregated breakpoints
