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
    """
    Load exon information from TSV file.
    Handles both files with and without headers.
    
    Expected columns (in order):
    chromosome, start, end, gene, exon, strand
    """
    # Expected column names (used when file has no headers)
    expected_columns = ['chromosome', 'start', 'end', 'gene', 'exon', 'strand']
    
    # Read first line to check if it looks like headers
    with open(exon_file, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            raise ValueError(f"File {exon_file} appears to be empty")
        first_values = [val.strip() for val in first_line.split('\t')]
        
        # Check if first value exactly matches expected header column names
        # This is more reliable than checking for keywords
        first_val_lower = first_values[0].lower() if len(first_values) > 0 else ''
        looks_like_headers = first_val_lower in [col.lower() for col in expected_columns]
        
        # Check if first line looks like data
        # Data would have: chromosome identifier (chr7, chr1, etc.) as first value
        # and numeric values for start/end positions
        looks_like_data = False
        if len(first_values) >= 3:
            # Check if first value is a chromosome identifier (chr followed by number or X/Y)
            first_val = first_values[0].lower()
            is_chr_id = (
                first_val.startswith('chr') and 
                len(first_val) > 3 and
                (first_val[3:].isdigit() or first_val[3:] in ['x', 'y', 'm', 'mt'])
            )
            # Check if second and third values look like numeric coordinates
            try:
                int(first_values[1])
                int(first_values[2])
                looks_like_data = is_chr_id
            except ValueError:
                looks_like_data = False
    
    # Read file based on header detection
    if looks_like_headers and not looks_like_data:
        # File has headers - read with headers
        exons_df = pd.read_csv(exon_file, sep='\t')
        # Strip whitespace from column names
        exons_df.columns = exons_df.columns.str.strip()
    else:
        # File doesn't have headers - read without headers and assign expected names
        exons_df = pd.read_csv(exon_file, sep='\t', header=None, names=expected_columns)
    
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
    Find which exon or intron a genomic position falls into or near.
    
    Args:
        position: Genomic coordinate
        exons_df: DataFrame with exon information
        chromosome: Chromosome name
        strand: Strand orientation
        window: Window size for boundary matching (default: 50bp)
        max_intron_dist: Maximum distance from exon for intronic breakpoints (default: 50kb)
    
    Returns:
        Tuple of (location_string, distance, is_within) where:
        - location_string: "exon X" or "intron X" format
        - distance: Distance to exon (0 if within, >0 if near)
        - is_within: True if within exon, False if near
        Returns (None, None, None) if not found
    """
    # Filter for matching chromosome and strand
    matching_exons = exons_df[
        (exons_df['chromosome'] == chromosome) &
        (exons_df['strand'] == strand)
    ].copy()
    
    if len(matching_exons) == 0:
        return None, None, None
    
    # Sort exons by exon number
    matching_exons = matching_exons.sort_values('exon')
    
    # Check if position is within exon boundaries
    within_exon = matching_exons[
        (matching_exons['start'] <= position) &
        (matching_exons['end'] >= position)
    ]
    
    if len(within_exon) > 0:
        exon_num = int(within_exon.iloc[0]['exon'])
        return f"exon {exon_num}", 0, True
    
    # Check if position is near exon boundaries (within window)
    # Near start boundary
    near_start = matching_exons[
        (matching_exons['start'] - window <= position) &
        (matching_exons['start'] + window >= position)
    ]
    
    if len(near_start) > 0:
        exon_num = int(near_start.iloc[0]['exon'])
        dist = abs(position - near_start.iloc[0]['start'])
        return f"exon {exon_num}", dist, False
    
    # Near end boundary
    near_end = matching_exons[
        (matching_exons['end'] - window <= position) &
        (matching_exons['end'] + window >= position)
    ]
    
    if len(near_end) > 0:
        exon_num = int(near_end.iloc[0]['exon'])
        dist = abs(position - near_end.iloc[0]['end'])
        return f"exon {exon_num}", dist, False
    
    # Check if position is within or near gene span (between first and last exon)
    gene_start = matching_exons['start'].min()
    gene_end = matching_exons['end'].max()
    gene_span_buffer = 10000  # 10kb buffer around gene span
    
    # Check if position is within gene span or within buffer
    if (gene_start - gene_span_buffer) <= position <= (gene_end + gene_span_buffer):
        # Find closest exon and determine if in intron
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
            distances.append((dist_to_exon, int(exon_row['exon']), exon_row['start'], exon_row['end']))
        
        if distances:
            min_dist, closest_exon, closest_start, closest_end = min(distances, key=lambda x: x[0])
            if min_dist <= max_intron_dist:
                # Check if position is in an intron (between two exons)
                # Intron N is between exon N and exon N+1
                for i in range(len(matching_exons) - 1):
                    exon_i = matching_exons.iloc[i]
                    exon_i_plus_1 = matching_exons.iloc[i + 1]
                    
                    # Check if position is between exon i and exon i+1
                    if exon_i['end'] < position < exon_i_plus_1['start']:
                        # Position is in intron i (between exon i and exon i+1)
                        intron_num = int(exon_i['exon'])
                        # Calculate distance to nearest exon boundary
                        dist_to_exon_i_end = position - exon_i['end']
                        dist_to_exon_i_plus_1_start = exon_i_plus_1['start'] - position
                        min_intron_dist = min(dist_to_exon_i_end, dist_to_exon_i_plus_1_start)
                        return f"intron {intron_num}", min_intron_dist, False
                
                # Not in an intron, but near an exon
                return f"exon {closest_exon}", min_dist, False
    
    return None, None, None


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
        # Get location information - use most common location or first if all same
        bp1_locations = group['breakpoint1_location'].tolist()
        bp2_locations = group['breakpoint2_location'].tolist()
        
        # Get most common location for each breakpoint
        from collections import Counter
        bp1_loc_counter = Counter(bp1_locations)
        bp2_loc_counter = Counter(bp2_locations)
        
        bp1_most_common_loc = bp1_loc_counter.most_common(1)[0][0]
        bp2_most_common_loc = bp2_loc_counter.most_common(1)[0][0]
        
        # Get distance and within status for the most common location
        bp1_rows = group[group['breakpoint1_location'] == bp1_most_common_loc]
        bp2_rows = group[group['breakpoint2_location'] == bp2_most_common_loc]
        
        # Use average distance if multiple reads have same location but different distances
        bp1_avg_dist = bp1_rows['breakpoint1_dist_to_exon'].mean()
        bp2_avg_dist = bp2_rows['breakpoint2_dist_to_exon'].mean()
        
        bp1_is_within = bp1_rows['breakpoint1_is_within'].iloc[0] if len(bp1_rows) > 0 else False
        bp2_is_within = bp2_rows['breakpoint2_is_within'].iloc[0] if len(bp2_rows) > 0 else False
        
        # Format location with distance indicator
        def format_location(location_str, distance, is_within):
            """Format location string with distance indicator."""
            if is_within:
                return location_str
            else:
                return f"{location_str} ({int(distance)}bp)"
        
        left_location_anno = format_location(bp1_most_common_loc, bp1_avg_dist, bp1_is_within)
        right_location_anno = format_location(bp2_most_common_loc, bp2_avg_dist, bp2_is_within)
        
        # Get exon numbers for sorting (extract from location string)
        def extract_exon_num(loc_str):
            if loc_str.startswith('exon '):
                return int(loc_str.split()[1])
            elif loc_str.startswith('intron '):
                return int(loc_str.split()[1])
            return 0
        
        exon1_num = extract_exon_num(bp1_most_common_loc)
        exon2_num = extract_exon_num(bp2_most_common_loc)
        
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
            'left_location_anno': left_location_anno,
            'right_location_anno': right_location_anno,
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
            'left_exon_num': exon1_num,  # For sorting
            'right_exon_num': exon2_num,  # For sorting
        })
    
    summary_df = pd.DataFrame(summary_rows)
    
    # Sort by penalty score (non-continuous first), then by exon numbers
    if 'left_exon_num' in summary_df.columns:
        summary_df = summary_df.sort_values(
            by=['penalty_score', 'left_exon_num', 'right_exon_num'],
            ascending=[True, True, True]
        )
    else:
        summary_df = summary_df.sort_values(
            by=['penalty_score', 'left_location_anno', 'right_location_anno'],
            ascending=[True, True, True]
        )
    
    return summary_df


def process_single_gene(exons_df, junctions_df, gene_name):
    """
    Process fusions for a single gene.
    
    Args:
        exons_df: DataFrame with exon information (filtered for one gene)
        junctions_df: DataFrame with chimeric junction data
        gene_name: Name of the gene being processed
    
    Returns:
        Tuple of (summary_df, detail_df) or (None, None) if no fusions found
    """
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
        
        # Find location for both breakpoints
        loc1, dist1, is_within1 = find_exon_for_position(bp1_pos, exons_df, chrom1, strand1)
        loc2, dist2, is_within2 = find_exon_for_position(bp2_pos, exons_df, chrom2, strand2)
        
        # Check if both breakpoints are in exons or introns
        if loc1 is None or loc2 is None:
            continue
        
        # Extract exon numbers for comparison (handle both "exon X" and "intron X" formats)
        def extract_exon_number(location_str):
            """Extract exon number from location string for comparison."""
            if location_str.startswith('exon '):
                return int(location_str.split()[1])
            elif location_str.startswith('intron '):
                # For intron N, it's between exon N and exon N+1, so we use N for comparison
                return int(location_str.split()[1])
            return None
        
        exon1_num = extract_exon_number(loc1)
        exon2_num = extract_exon_number(loc2)
        
        # Check if it's the same location (not a fusion)
        if loc1 == loc2:
            continue
        
        # Determine if it's continuous or non-continuous (only for exons)
        is_continuous = False
        if exon1_num is not None and exon2_num is not None:
            # Only check continuity if both are exons (not introns)
            if 'exon' in loc1 and 'exon' in loc2:
                is_continuous = is_continuous_exon_pair(exon1_num, exon2_num)
        
        # Calculate distance between locations
        if exon1_num is not None and exon2_num is not None:
            exon_distance = abs(exon1_num - exon2_num)
        else:
            exon_distance = 0
        
        # Format location with distance indicator
        def format_location(location_str, distance, is_within):
            """Format location string with distance indicator."""
            if is_within:
                return location_str
            else:
                return f"{location_str} ({distance}bp)"
        
        loc1_formatted = format_location(loc1, dist1, is_within1)
        loc2_formatted = format_location(loc2, dist2, is_within2)
        
        # Create fusion type string
        fusion_type = f"{loc1}-{loc2}"
        
        # Create result entry
        result = {
            'read_name': read_name,
            'chromosome': chrom1,
            'breakpoint1_coord': bp1_pos,
            'breakpoint1_location': loc1,
            'breakpoint1_location_formatted': loc1_formatted,
            'breakpoint1_dist_to_exon': dist1 if dist1 is not None else 0,
            'breakpoint1_is_within': is_within1,
            'breakpoint1_exon_num': exon1_num,  # For sorting/comparison
            'breakpoint1_strand': strand1,
            'breakpoint1_cigar': str(cigar1),
            'breakpoint2_coord': bp2_pos,
            'breakpoint2_location': loc2,
            'breakpoint2_location_formatted': loc2_formatted,
            'breakpoint2_dist_to_exon': dist2 if dist2 is not None else 0,
            'breakpoint2_is_within': is_within2,
            'breakpoint2_exon_num': exon2_num,  # For sorting/comparison
            'breakpoint2_strand': strand2,
            'breakpoint2_cigar': str(cigar2),
            'exon_distance': exon_distance,
            'is_continuous': is_continuous,
            'fusion_type': fusion_type,
            'penalty_score': 0 if not is_continuous else 1,  # Penalty for continuous
            'junction_type': row['junction_type'],
        }
        
        results.append(result)
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    if len(results_df) == 0:
        print(f"  No intra-gene fusions found for {gene_name}.", file=sys.stderr)
        return None, None
    
    # Create aggregated summary output
    summary_df = create_summary_output(results_df, gene_name)
    
    # Print summary to stderr
    print(f"  Found {len(summary_df)} unique fusion types with {len(results_df)} total supporting reads:", file=sys.stderr)
    print(f"    - Non-continuous (exon distance > 1): {(summary_df['penalty_score'] == 0).sum()}", file=sys.stderr)
    print(f"    - Continuous (adjacent exons, penalized): {(summary_df['penalty_score'] == 1).sum()}", file=sys.stderr)
    
    # Define columns for detail output (all columns)
    detail_columns = [
        'gene',
        'supporting_reads',
        'left_location_anno',
        'right_location_anno',
        'left_breakpoint_aggregated',
        'right_breakpoint_aggregated',
        'left_strand',
        'right_strand',
        'breakpoint1_positions',
        'breakpoint2_positions',
        'breakpoint1_read_ids',
        'breakpoint2_read_ids'
    ]
    
    # Define columns for summary output (without read alignment details and read IDs)
    summary_columns = [
        'gene',
        'supporting_reads',
        'left_location_anno',
        'right_location_anno',
        'left_breakpoint_aggregated',
        'right_breakpoint_aggregated',
        'left_strand',
        'right_strand'
    ]
    
    # Sort before selecting columns (so we can use exon_num for sorting)
    if 'left_exon_num' in summary_df.columns:
        summary_df = summary_df.sort_values(
            by=['penalty_score', 'left_exon_num', 'right_exon_num'],
            ascending=[True, True, True]
        )
    
    detail_df = summary_df[detail_columns].copy()
    summary_output_df = summary_df[summary_columns].copy()
    
    return summary_output_df, detail_df


def identify_intra_gene_fusions(exons_file, junctions_file, output_base='fusions_itd'):
    """
    Identify intra-gene fusions with non-continuous exon breakpoints.
    Processes multiple genes if present in the exon file.
    
    Args:
        exons_file: Path to exons TSV file (can contain multiple genes)
        junctions_file: Path to chimeric junctions TSV file
        output_base: Base name for output files (default: 'fusions_itd')
                    Creates: {output_base}.summary.tsv and {output_base}.detail.tsv
    """
    # Load data
    print("Loading exon information...", file=sys.stderr)
    exons_df = load_exons(exons_file)
    
    print("Loading chimeric junctions...", file=sys.stderr)
    junctions_df = load_chimeric_junctions(junctions_file)
    
    # Get unique genes
    unique_genes = exons_df['gene'].unique()
    print(f"\nFound {len(unique_genes)} gene(s) to process: {', '.join(unique_genes)}", file=sys.stderr)
    
    # Store results for all genes
    all_summary_dfs = []
    all_detail_dfs = []
    
    # Process each gene
    for gene_name in unique_genes:
        # Filter exons for this gene
        gene_exons_df = exons_df[exons_df['gene'] == gene_name].copy()
        
        # Process this gene
        summary_df, detail_df = process_single_gene(gene_exons_df, junctions_df, gene_name)
        
        if summary_df is not None and detail_df is not None:
            all_summary_dfs.append(summary_df)
            all_detail_dfs.append(detail_df)
    
    # Combine results from all genes
    if len(all_summary_dfs) == 0:
        print("\nNo intra-gene fusions found for any gene.", file=sys.stderr)
        # Create empty files with headers
        summary_columns = [
            'gene', 'supporting_reads', 'left_location_anno', 'right_location_anno',
            'left_breakpoint_aggregated', 'right_breakpoint_aggregated',
            'left_strand', 'right_strand'
        ]
        detail_columns = [
            'gene', 'supporting_reads', 'left_location_anno', 'right_location_anno',
            'left_breakpoint_aggregated', 'right_breakpoint_aggregated',
            'left_strand', 'right_strand', 'breakpoint1_positions', 'breakpoint2_positions',
            'breakpoint1_read_ids', 'breakpoint2_read_ids'
        ]
        empty_summary = pd.DataFrame(columns=summary_columns)
        empty_detail = pd.DataFrame(columns=detail_columns)
        summary_file = f"{output_base}.summary.tsv"
        detail_file = f"{output_base}.detail.tsv"
        empty_summary.to_csv(summary_file, sep='\t', index=False)
        empty_detail.to_csv(detail_file, sep='\t', index=False)
        print(f"\nEmpty output files created: {summary_file}, {detail_file}", file=sys.stderr)
        return
    
    # Concatenate all results
    combined_summary = pd.concat(all_summary_dfs, ignore_index=True)
    combined_detail = pd.concat(all_detail_dfs, ignore_index=True)
    
    # Sort by gene, then by exon numbers (if available)
    # Check if sorting columns exist
    sort_cols = ['gene']
    if 'left_exon_num' in combined_summary.columns:
        sort_cols.extend(['left_exon_num', 'right_exon_num'])
    else:
        sort_cols.extend(['left_location_anno', 'right_location_anno'])
    
    combined_summary = combined_summary.sort_values(
        by=sort_cols,
        ascending=[True] * len(sort_cols)
    )
    combined_detail = combined_detail.sort_values(
        by=sort_cols,
        ascending=[True] * len(sort_cols)
    )
    
    # Write output files
    summary_file = f"{output_base}.summary.tsv"
    detail_file = f"{output_base}.detail.tsv"
    
    combined_summary.to_csv(summary_file, sep='\t', index=False)
    combined_detail.to_csv(detail_file, sep='\t', index=False)
    
    print(f"\nResults written to:", file=sys.stderr)
    print(f"  Summary (without read alignment details): {summary_file}", file=sys.stderr)
    print(f"  Detail (with all columns): {detail_file}", file=sys.stderr)
    print(f"\nTotal fusions found: {len(combined_summary)} unique fusion types", file=sys.stderr)
    
    return combined_summary, combined_detail


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Identify intra-gene fusions with non-continuous exon breakpoints'
    )
    parser.add_argument(
        'exons_file',
        help='Path to exons TSV file (can contain multiple genes)'
    )
    parser.add_argument(
        'junctions_file',
        help='Path to chimeric junctions TSV file (Chimeric.out.junction.tsv)'
    )
    parser.add_argument(
        '-o', '--output',
        help='Base name for output files (default: fusions_itd). Creates {base}.summary.tsv and {base}.detail.tsv',
        default='fusions_itd'
    )
    
    args = parser.parse_args()
    
    identify_intra_gene_fusions(
        args.exons_file,
        args.junctions_file,
        args.output
    )
