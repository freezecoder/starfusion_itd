# Intra-Gene Fusion Detection Tool

A Python tool for identifying intra-gene fusions (ITDs - Internal Tandem Duplications) from chimeric junction data. This tool processes STAR chimeric junction output to detect fusion events between non-continuous exons within the same gene.

## Features

- **Intra-gene fusion detection**: Identifies fusions between different exons within the same gene
- **Non-continuous exon preference**: Prioritizes fusions between non-adjacent exons
- **Breakpoint span calculation**: Calculates the full span of breakpoint positions from supporting reads
- **Multiple output formats**: TSV, CSV, and BEDPE formats for different use cases
- **Genome browser visualization**: BEDPE output for viewing fusions in IGV or UCSC Genome Browser
- **Read filtering**: Configurable minimum supporting reads threshold
- **Multi-gene support**: Processes multiple genes from a single exon annotation file
- **Robust input handling**: Handles files with or without headers, comment lines, and various chromosome naming conventions

## Requirements

- Python 3.6 or higher
- pandas (tested with pandas >= 1.0.0)

## Installation

### Option 1: Direct Download

1. Download the script:
   ```bash
   wget https://github.com/yourusername/egfr_itd_dev/identify_intra_gene_fusions.py
   ```

2. Make it executable:
   ```bash
   chmod +x identify_intra_gene_fusions.py
   ```

### Option 2: Clone Repository

```bash
git clone https://github.com/yourusername/egfr_itd_dev.git
cd egfr_itd_dev
```

### Install Dependencies

```bash
pip install pandas
```

Or using conda:

```bash
conda install pandas
```

## Usage

### Basic Usage

```bash
python identify_intra_gene_fusions.py <exons_file> <junctions_file>
```

### Required Arguments

- `exons_file`: Path to exon annotation file (TSV format)
- `junctions_file`: Path to STAR Chimeric.out.junction.tsv file

### Optional Arguments

- `-o, --output OUTPUT`: Base name for output files (default: `fusions_itd`)
  - Creates: `{output_base}.summary.{ext}`, `{output_base}.detail.{ext}`, `{output_base}.bedpe`
  
- `--format {tsv,csv}`: Output file format (default: `tsv`)
  - `tsv`: Tab-separated values
  - `csv`: Comma-separated values
  
- `--min-reads MIN_READS`: Minimum number of supporting reads required (default: `10`)
  - Only fusion events with at least this many supporting reads will be reported

### Examples

#### Basic usage with default settings:
```bash
python identify_intra_gene_fusions.py EGFR_exons_hg19.tsv Chimeric.out.junction.tsv
```

This creates:
- `fusions_itd.summary.tsv`
- `fusions_itd.detail.tsv`
- `fusions_itd.bedpe`

#### Custom output name and CSV format:
```bash
python identify_intra_gene_fusions.py \
    EGFR_exons_hg19.tsv \
    Chimeric.out.junction.tsv \
    -o my_results \
    --format csv
```

This creates:
- `my_results.summary.csv`
- `my_results.detail.csv`
- `my_results.bedpe`

#### Filter with minimum 5 supporting reads:
```bash
python identify_intra_gene_fusions.py \
    EGFR_exons_hg19.tsv \
    Chimeric.out.junction.tsv \
    --min-reads 5
```

#### All options combined:
```bash
python identify_intra_gene_fusions.py \
    EGFR_exons_hg19.tsv \
    Chimeric.out.junction.tsv \
    -o egfr_itd_results \
    --format csv \
    --min-reads 15
```

## Input File Formats

### Exon Annotation File

A tab-separated file containing exon coordinates. Can have headers or not.

**Required columns** (in order if no headers):
- `chromosome`: Chromosome name (e.g., `chr7`, `7`)
- `start`: Exon start position (1-based)
- `end`: Exon end position (1-based)
- `gene`: Gene name
- `exon`: Exon number
- `strand`: Strand orientation (`+` or `-`)

**Example** (with headers):
```tsv
chromosome	start	end	gene	exon	strand
chr7	55242464	55242521	EGFR	1	+
chr7	55248914	55249023	EGFR	2	+
```

**Example** (without headers):
```tsv
chr7	55242464	55242521	EGFR	1	+
chr7	55248914	55249023	EGFR	2	+
```

### Chimeric Junction File

STAR's `Chimeric.out.junction.tsv` file format. The tool handles:
- Files with or without headers
- Comment lines starting with `#` (automatically skipped)
- Various chromosome naming conventions (`chr7`, `7`, `Chr7`)

**Expected columns** (if headers present):
- `chrom1`, `coord1`, `strand1`: First breakpoint chromosome, coordinate, strand
- `chrom2`, `coord2`, `strand2`: Second breakpoint chromosome, coordinate, strand
- `breakpoint1_pos`, `CIGAR1`: Precise breakpoint position and CIGAR for breakpoint 1
- `breakpoint2_pos`, `CIGAR2`: Precise breakpoint position and CIGAR for breakpoint 2
- `read_name`: Read identifier
- `junction_type`: Junction type (0, 1, or 2)

**Note**: If the file has no headers, the tool assumes standard STAR column order.

## Output Files

### Summary File (`{output_base}.summary.{ext}`)

Contains fusion events without read-level details. Columns:

- `gene`: Gene name
- `supporting_reads`: Number of reads supporting this fusion
- `left_location_anno`: Left breakpoint location (e.g., "exon 21" or "exon 23 (150bp)")
- `right_location_anno`: Right breakpoint location
- `left_breakpoint_aggregated`: Consensus left breakpoint position (`chr:pos`)
- `right_breakpoint_aggregated`: Consensus right breakpoint position (`chr:pos`)
- `left_fusion_pos`: Left breakpoint span (`chr:min-max`)
- `right_fusion_pos`: Right breakpoint span (`chr:min-max`)
- `left_strand`: Left breakpoint strand (`+` or `-`)
- `right_strand`: Right breakpoint strand (`+` or `-`)

### Detail File (`{output_base}.detail.{ext}`)

Contains all information including read-level details. Includes all summary columns plus:

- `breakpoint1_positions`: All breakpoint1 positions with CIGARs (`pos:CIGAR,pos:CIGAR,...`)
- `breakpoint2_positions`: All breakpoint2 positions with CIGARs
- `breakpoint1_read_ids`: Comma-separated list of read IDs supporting breakpoint1
- `breakpoint2_read_ids`: Comma-separated list of read IDs supporting breakpoint2

### BEDPE File (`{output_base}.bedpe`)

Standard BEDPE format for genome browser visualization. Format:

```
chr1  start1  end1  chr2  start2  end2  name  score  strand1  strand2
```

- Uses 0-based start coordinates (BEDPE standard)
- Fusion names: `gene_left_location_right_location`
- Score: Number of supporting reads
- Breakpoint spans represent min-max range from supporting reads

**Loading in IGV**:
1. File → Load from File → Select `.bedpe` file
2. Fusions will appear as arcs connecting breakpoint regions

**Loading in UCSC Genome Browser**:
1. Add Custom Tracks → Upload file
2. Select "BEDPE" format
3. Fusions will appear as paired-end intervals

## Algorithm Overview

1. **Data Loading**: Loads exon annotations and chimeric junctions
2. **Junction Processing**: For each junction:
   - Extracts breakpoint positions (prefers `breakpoint1_pos`/`breakpoint2_pos`, falls back to `coord1`/`coord2`)
   - Matches breakpoints to exons based on chromosome, position, and strand
   - Filters for intra-chromosomal, intra-gene fusions
3. **Exon Assignment**: Determines which exon each breakpoint belongs to (within exon, near exon, or in intron)
4. **Fusion Classification**: Classifies as continuous (adjacent exons) or non-continuous (preferred)
5. **Aggregation**: Groups junctions by fusion type and calculates:
   - Most common breakpoint positions (mode)
   - Breakpoint position spans (min-max range)
   - Supporting read counts
6. **Filtering**: Filters by minimum supporting reads threshold
7. **Output**: Generates summary, detail, and BEDPE files

For detailed algorithm documentation, see [ALGORITHM.md](ALGORITHM.md).

## Filtering and Statistics

The tool provides detailed filtering statistics:

```
Filtering by minimum supporting reads (min_reads=10):
  Total fusion events before filtering: 25
  Events filtered out (< 10 reads): 5
  Events remaining after filtering: 20
```

This helps you understand how many events were filtered and adjust the `--min-reads` threshold if needed.

## Troubleshooting

### "No intra-gene fusions found"

Possible causes:
1. **No junctions match gene chromosomes**: Check that your chimeric file contains junctions on the same chromosome as your gene
2. **Strand mismatch**: Ensure breakpoint strands match exon strand annotations
3. **Minimum reads threshold too high**: Try lowering `--min-reads`
4. **Column misalignment**: The tool auto-detects and corrects column misalignment, but check diagnostic messages

### "KeyError: 'chrom1'"

This indicates column detection issues. The tool should auto-correct, but if it persists:
- Check that your chimeric file follows STAR format
- Ensure the file is tab-separated
- Verify comment lines start with `#`

### Performance Issues

For large files (>100K junctions):
- The tool uses optimized pandas operations
- Processing time scales with number of junctions and genes
- Consider pre-filtering junctions by chromosome if processing multiple genes

## Citation

If you use this tool in your research, please cite:

```
Intra-Gene Fusion Detection Tool
[Your citation information]
```

## License

[Specify your license]

## Contributing

[Contributing guidelines]

## Support

For issues, questions, or feature requests, please open an issue on GitHub.

## Version History

- **v1.0**: Initial release
  - Basic intra-gene fusion detection
  - TSV/CSV output formats
  - BEDPE output for genome browser visualization
  - Minimum reads filtering
  - Breakpoint span calculation
