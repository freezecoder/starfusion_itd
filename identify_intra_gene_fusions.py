#!/usr/bin/env python3
"""
Identify intra-gene fusions between non-continuous exons from chimeric junction data.

This script identifies fusions within the same gene (EGFR) that involve
non-continuous exon breakpoints. Continuous exon breakpoints (e.g., exon2-exon3)
are penalized but still reported.
"""

import pandas as pd
import sys
from collections import defaultdict


def load_exons(exon_file):
    """Load exon information from TSV file."""
    exons_df = pd.read_csv(exon_file, sep='\t')
    return exons_df


def load_chimeric_junctions(junction_file):
    """
    Load chimeric junction data from TSV file.
    Handles both files with and without headers.
    
    Expected columns (in order):
    chrom1, coord1, strand1, chrom2, coord2, strand2, junction_type, 
    repeat_left, repeat_right, read_name, breakpoint1_pos, CIGAR1, 
    breakpoint2_pos, CIGAR2
    """
    # Expected column names (used when file has no headers)
    expected_columns = [
        'chrom1', 'coord1', 'strand1', 'chrom2', 'coord2', 'strand2',
        'junction_type', 'repeat_left', 'repeat_right', 'read_name',
        'breakpoint1_pos', 'CIGAR1', 'breakpoint2_pos', 'CIGAR2'
    ]
    
    # Read first line to check if it looks like headers
    with open(junction_file, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            raise ValueError(f"File {junction_file} appears to be empty")
        first_values = [val.strip() for val in first_line.split('\t')]
        
        # Check if first line looks like headers
        # Headers typically contain words like 'chrom', 'coord', 'strand', 'read_name'
        header_keywords = ['chrom', 'coord', 'strand', 'read', 'cigar', 'junction', 'repeat']
        first_line_lower = first_line.lower()
        looks_like_headers = any(keyword in first_line_lower for keyword in header_keywords)
        
        # Check if first line looks like data (starts with 'chr' followed by chromosome identifier)
        looks_like_data = (
            len(first_values) > 0 and 
            (first_values[0].startswith('chr') or first_values[0].startswith('Chr'))
        )
    
    # Read file based on header detection
    if looks_like_headers and not looks_like_data:
        # File has headers - read with headers
        junctions_df = pd.read_csv(junction_file, sep='\t')
        # Strip whitespace from column names
        junctions_df.columns = junctions_df.columns.str.strip()
    else:
        # File doesn't have headers - read without headers and assign expected names
        junctions_df = pd.read_csv(junction_file, sep='\t', header=None, names=expected_columns)
    
    return junctions_df


def find_exon_for_position(position, exons_df, chromosome, strand, window=50, max_intron_dist=50000):
    """
    Find which exon a genomic position falls into or near.
    
    Args:
        position: Genomic coordinate
        exons_df: DataFrame with exon information
        chromosome: Chromosome name
        strand: Strand orientation
        window: Window size for boundary matching (default: 50bp)
        max_intron_dist: Maximum distance from exon for intronic breakpoints (default: 50kb)
    
    Returns:
        Tuple of (exon number, distance_to_exon) if found, (None, None) otherwise
    """
    # Filter for matching chromosome and strand
    matching_exons = exons_df[
        (exons_df['chromosome'] == chromosome) &
        (exons_df['strand'] == strand)
    ].copy()
    
    if len(matching_exons) == 0:
        return None, None
    
    # Check if position is within exon boundaries
    within_exon = matching_exons[
        (matching_exons['start'] <= position) &
        (matching_exons['end'] >= position)
    ]
    
    if len(within_exon) > 0:
        return int(within_exon.iloc[0]['exon']), 0
    
    # Check if position is near exon boundaries (within window)
    # Near start boundary
    near_start = matching_exons[
        (matching_exons['start'] - window <= position) &
        (matching_exons['start'] + window >= position)
    ]
    
    if len(near_start) > 0:
        dist = abs(position - near_start.iloc[0]['start'])
        return int(near_start.iloc[0]['exon']), dist
    
    # Near end boundary
    near_end = matching_exons[
        (matching_exons['end'] - window <= position) &
        (matching_exons['end'] + window >= position)
    ]
    
    if len(near_end) > 0:
        dist = abs(position - near_end.iloc[0]['end'])
        return int(near_end.iloc[0]['exon']), dist
    
    # Check if position is within or near gene span (between first and last exon)
    gene_start = matching_exons['start'].min()
    gene_end = matching_exons['end'].max()
    gene_span_buffer = 10000  # 10kb buffer around gene span
    
    # Check if position is within gene span or within buffer
    if (gene_start - gene_span_buffer) <= position <= (gene_end + gene_span_buffer):
        # Position is within or near gene span, find closest exon
        distances = []
        for _, exon_row in matching_exons.iterrows():
            # Distance to exon start
            dist_to_start = abs(position - exon_row['start'])
            # Distance to exon end
            dist_to_end = abs(position - exon_row['end'])
            # Distance to exon region (0 if within, else min distance to boundary)
            if exon_row['start'] <= position <= exon_row['end']:
                dist_to_exon = 0
            else:
                dist_to_exon = min(dist_to_start, dist_to_end)
            distances.append((dist_to_exon, int(exon_row['exon'])))
        
        if distances:
            min_dist, closest_exon = min(distances, key=lambda x: x[0])
            if min_dist <= max_intron_dist:
                return closest_exon, min_dist
    
    return None, None


def is_continuous_exon_pair(exon1, exon2):
    """
    Check if two exons are continuous (adjacent).
    
    Args:
        exon1: First exon number
        exon2: Second exon number
    
    Returns:
        True if exons are continuous (e.g., exon2-exon3), False otherwise
    """
    if exon1 is None or exon2 is None:
        return False
    
    # Check if exons are adjacent (difference of 1)
    return abs(exon1 - exon2) == 1


def create_summary_output(results_df, gene_name):
    """
    Create aggregated summary output grouped by fusion type.
    
    Args:
        results_df: DataFrame with individual read results
        gene_name: Name of the gene
    
    Returns:
        DataFrame with aggregated summary by fusion type
    """
    summary_rows = []
    
    # Group by fusion_type
    for fusion_type, group in results_df.groupby('fusion_type'):
        # Extract exon numbers from fusion_type (e.g., "exon21-exon18")
        exon1 = group['breakpoint1_exon'].iloc[0]
        exon2 = group['breakpoint2_exon'].iloc[0]
        
        # Get strand information
        strand1 = group['breakpoint1_strand'].iloc[0]
        strand2 = group['breakpoint2_strand'].iloc[0]
        
        # Get chromosome (should be same for both breakpoints in intra-gene fusions)
        chromosome = group['chromosome'].iloc[0]
        
        # Calculate aggregated breakpoint positions (most common position)
        bp1_positions_list = group['breakpoint1_coord'].tolist()
        bp2_positions_list = group['breakpoint2_coord'].tolist()
        
        # Find most common position (mode) for each breakpoint
        from collections import Counter
        bp1_pos_counter = Counter(bp1_positions_list)
        bp2_pos_counter = Counter(bp2_positions_list)
        
        # Get most common position, or median if tie
        bp1_agg_pos = bp1_pos_counter.most_common(1)[0][0]
        bp2_agg_pos = bp2_pos_counter.most_common(1)[0][0]
        
        # Format as chr:pos
        bp1_aggregated = f"{chromosome}:{bp1_agg_pos}"
        bp2_aggregated = f"{chromosome}:{bp2_agg_pos}"
        
        # Count supporting reads
        supporting_reads = len(group)
        
        # Collect positions with CIGARs for breakpoint1
        # Format: position:CIGAR,position:CIGAR,...
        bp1_pos_cigar_pairs = []
        for _, row in group.iterrows():
            pos = row['breakpoint1_coord']
            cigar = row.get('breakpoint1_cigar', '')
            bp1_pos_cigar_pairs.append(f"{pos}:{cigar}")
        # Remove duplicates while preserving order
        seen = set()
        bp1_pos_cigar_unique = []
        for item in bp1_pos_cigar_pairs:
            if item not in seen:
                seen.add(item)
                bp1_pos_cigar_unique.append(item)
        bp1_positions_cigars_str = ','.join(bp1_pos_cigar_unique)
        
        # Collect positions with CIGARs for breakpoint2
        bp2_pos_cigar_pairs = []
        for _, row in group.iterrows():
            pos = row['breakpoint2_coord']
            cigar = row.get('breakpoint2_cigar', '')
            bp2_pos_cigar_pairs.append(f"{pos}:{cigar}")
        # Remove duplicates while preserving order
        seen = set()
        bp2_pos_cigar_unique = []
        for item in bp2_pos_cigar_pairs:
            if item not in seen:
                seen.add(item)
                bp2_pos_cigar_unique.append(item)
        bp2_positions_cigars_str = ','.join(bp2_pos_cigar_unique)
        
        # Collect all read IDs for each breakpoint
        bp1_read_ids = sorted(group['read_name'].unique())
        bp2_read_ids = sorted(group['read_name'].unique())
        
        # Format as comma-separated strings
        bp1_read_ids_str = ','.join(bp1_read_ids)
        bp2_read_ids_str = ','.join(bp2_read_ids)
        
        # Get penalty score and exon distance from first row
        penalty_score = group['penalty_score'].iloc[0]
        exon_distance = group['exon_distance'].iloc[0]
        
        summary_rows.append({
            'gene': gene_name,
            'supporting_reads': supporting_reads,
            'left_exon_anno': exon1,
            'right_exon_anno': exon2,
            'left_breakpoint_aggregated': bp1_aggregated,
            'right_breakpoint_aggregated': bp2_aggregated,
            'left_strand': strand1,
            'right_strand': strand2,
            'breakpoint1_positions': bp1_positions_cigars_str,
            'breakpoint2_positions': bp2_positions_cigars_str,
            'breakpoint1_read_ids': bp1_read_ids_str,
            'breakpoint2_read_ids': bp2_read_ids_str,
            'penalty_score': penalty_score,
            'exon_distance': exon_distance,
        })
    
    summary_df = pd.DataFrame(summary_rows)
    
    # Sort by penalty score (non-continuous first) and exon distance
    summary_df = summary_df.sort_values(
        by=['penalty_score', 'exon_distance'],
        ascending=[True, False]
    )
    
    return summary_df


def identify_intra_gene_fusions(exons_file, junctions_file, output_file=None):
    """
    Identify intra-gene fusions with non-continuous exon breakpoints.
    
    Args:
        exons_file: Path to EGFR exons TSV file
        junctions_file: Path to chimeric junctions TSV file
        output_file: Optional path to output file (if None, prints to stdout)
    """
    # Load data
    print("Loading exon information...", file=sys.stderr)
    exons_df = load_exons(exons_file)
    
    print("Loading chimeric junctions...", file=sys.stderr)
    junctions_df = load_chimeric_junctions(junctions_file)
    
    # Get gene name from exons (assuming all exons are from same gene)
    gene_name = exons_df['gene'].iloc[0]
    print(f"Analyzing fusions for gene: {gene_name}", file=sys.stderr)
    
    # Store results
    results = []
    
    # Process each junction
    for idx, row in junctions_df.iterrows():
        chrom1 = row['chrom1']
        coord1 = row['coord1']
        strand1 = row['strand1']
        chrom2 = row['chrom2']
        coord2 = row['coord2']
        strand2 = row['strand2']
        read_name = row['read_name']
        
        # Use breakpoint positions if available, otherwise use coord1/coord2
        bp1_pos = row.get('breakpoint1_pos', coord1)
        bp2_pos = row.get('breakpoint2_pos', coord2)
        
        # Get CIGAR strings
        cigar1 = row.get('CIGAR1', '')
        cigar2 = row.get('CIGAR2', '')
        
        # Handle NaN values
        if pd.isna(bp1_pos):
            bp1_pos = coord1
        if pd.isna(bp2_pos):
            bp2_pos = coord2
        if pd.isna(cigar1):
            cigar1 = ''
        if pd.isna(cigar2):
            cigar2 = ''
        
        # Convert to int
        try:
            bp1_pos = int(bp1_pos)
            bp2_pos = int(bp2_pos)
        except (ValueError, TypeError):
            continue
        
        # Check if both breakpoints are on the same chromosome (intra-chromosomal)
        if chrom1 != chrom2:
            continue
        
        # Find exons for both breakpoints
        exon1, dist1 = find_exon_for_position(bp1_pos, exons_df, chrom1, strand1)
        exon2, dist2 = find_exon_for_position(bp2_pos, exons_df, chrom2, strand2)
        
        # Check if both breakpoints are in exons
        if exon1 is None or exon2 is None:
            continue
        
        # Check if it's the same exon (not a fusion)
        if exon1 == exon2:
            continue
        
        # Determine if it's continuous or non-continuous
        is_continuous = is_continuous_exon_pair(exon1, exon2)
        
        # Calculate exon distance
        exon_distance = abs(exon1 - exon2)
        
        # Create result entry
        result = {
            'read_name': read_name,
            'chromosome': chrom1,
            'breakpoint1_coord': bp1_pos,
            'breakpoint1_exon': exon1,
            'breakpoint1_dist_to_exon': dist1 if dist1 is not None else 0,
            'breakpoint1_strand': strand1,
            'breakpoint1_cigar': str(cigar1),
            'breakpoint2_coord': bp2_pos,
            'breakpoint2_exon': exon2,
            'breakpoint2_dist_to_exon': dist2 if dist2 is not None else 0,
            'breakpoint2_strand': strand2,
            'breakpoint2_cigar': str(cigar2),
            'exon_distance': exon_distance,
            'is_continuous': is_continuous,
            'fusion_type': f'exon{exon1}-exon{exon2}',
            'penalty_score': 0 if not is_continuous else 1,  # Penalty for continuous
            'junction_type': row['junction_type'],
        }
        
        results.append(result)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    if len(results_df) == 0:
        print("No intra-gene fusions found.", file=sys.stderr)
        return
    
    # Create aggregated summary output
    summary_df = create_summary_output(results_df, gene_name)
    
    # Print summary to stderr
    print(f"\nFound {len(summary_df)} unique fusion types with {len(results_df)} total supporting reads:", file=sys.stderr)
    print(f"  - Non-continuous (exon distance > 1): {(summary_df['penalty_score'] == 0).sum()}", file=sys.stderr)
    print(f"  - Continuous (adjacent exons, penalized): {(summary_df['penalty_score'] == 1).sum()}", file=sys.stderr)
    
    # Show unique fusion types
    if len(summary_df) > 0:
        print(f"\nFusion types found:", file=sys.stderr)
        for _, row in summary_df.iterrows():
            penalty_note = " (PENALIZED - continuous exons)" if row['penalty_score'] == 1 else ""
            print(f"  - exon{int(row['left_exon_anno'])}-exon{int(row['right_exon_anno'])}: {int(row['supporting_reads'])} reads{penalty_note}", file=sys.stderr)
    
    # Select only the required columns for output (excluding penalty_score and exon_distance)
    output_columns = [
        'gene',
        'supporting_reads',
        'left_exon_anno',
        'right_exon_anno',
        'left_breakpoint_aggregated',
        'right_breakpoint_aggregated',
        'left_strand',
        'right_strand',
        'breakpoint1_positions',
        'breakpoint2_positions',
        'breakpoint1_read_ids',
        'breakpoint2_read_ids'
    ]
    output_df = summary_df[output_columns].copy()
    
    # Output results
    if output_file:
        output_df.to_csv(output_file, sep='\t', index=False)
        print(f"\nResults written to: {output_file}", file=sys.stderr)
    else:
        # Print to stdout
        print("\n" + "="*100)
        print("INTRAGENE FUSION RESULTS")
        print("="*100)
        print(output_df.to_string(index=False))
    
    return output_df


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Identify intra-gene fusions with non-continuous exon breakpoints'
    )
    parser.add_argument(
        'exons_file',
        help='Path to EGFR exons TSV file (EGFR_exons_hg19.tsv)'
    )
    parser.add_argument(
        'junctions_file',
        help='Path to chimeric junctions TSV file (Chimeric.out.junction.tsv)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Output file path (default: print to stdout)',
        default=None
    )
    
    args = parser.parse_args()
    
    identify_intra_gene_fusions(
        args.exons_file,
        args.junctions_file,
        args.output
    )
