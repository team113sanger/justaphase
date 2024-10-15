import argparse
import logging
import collections
import gzip
import pysam
# from pathlib import Path


VCFLine = collections.namedtuple('VCFLine', ['chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'normal', 'tumour'])
FormatFields = collections.namedtuple('FormatFields', ['gt', 'forward_a_count', 'forward_c_count', 'forward_g_count', 'forward_t_count', 'reverse_a_count', 'reverse_c_count', 'reverse_g_count', 'reverse_t_count', 'proportion_mut', 'phase_set'])

def setup_logging():
    '''Setup logging configuration'''
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')


def get_args():
    '''Get command-line arguments'''
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--variants_file', nargs='+',help='One or more VCF files to be phased. Currently only supports Caveman VCFs compressed with bgzip.')
    parser.add_argument('-p', '--phased_variants', nargs='+',help='A file containing a set of phased variants.')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory.')
    
    return parser.parse_args()

def check_mnv_phasing(phased_variants_file) -> list:
    '''Use phasing information from phased VCF to check if adjacent SNVs occur on same chromosome'''
    phased_variants_file: str
    with open(phased_variants_file, 'r') as phased_variants_file:
        sample_mnvs = [line.strip() for line in phased_variants_file.readlines()]
        mnvs_to_keep = []
        # Loop through VCF lines and check if they have the same phase set
        mnv_phase_sets = []
        
        for line in sample_mnvs:
            vcf_line = VCFLine(*line.split('\t'))
            if len(vcf_line.tumour.split(':')) == 11: # Ensure the given line has a phase set
                tumour_format_fields = FormatFields(*vcf_line.tumour.split(':'))
                # Add the phase set from each variant to the set
                mnv_phase_sets.append(tumour_format_fields.phase_set)
            else:
                continue # Skip this line as it has no phase set
        
        if len(mnv_phase_sets) > 1: # Exclude MNVs where only one variant is in phase
            if len(set(mnv_phase_sets)) == 1: # All SNVs have the same phase set
                mnvs_to_keep.append([VCFLine(*line.split('\t')) for line in mnv_vcf_lines])
        
        if len(set(mnv_phase_sets)) == 1: # All SNVs have the same phase set
            mnvs_to_keep.append([VCFLine(*line.split('\t')) for line in mnv_vcf_lines])

    return mnvs_to_keep

def modify_vcf_lines(phased_vcf_file, mnvs_to_keep, output_vcf_file) -> None:
    '''Modify adjacent SNV lines in the phased VCF and output to a new VCF file'''
    phased_vcf_file: str
    mnvs_to_keep: list
    output_vcf_file: str
    
    logging.info('Converting adjacent phased SNVs into MNVs and outputting to VCF...')
    
    if not mnvs_to_keep:
        logging.info('No MNVs found. The output VCF will be identical to the input VCF.')
        # Simply copy the input VCF to the output file
        with gzip.open(phased_vcf_file, 'rb') as infile, gzip.open(output_vcf_file, 'wb') as outfile:
            outfile.write(infile.read())
        return
    
    # Create a dictionary where the keys are chromosome+position and the values are the MNV objects
    mnv_positions = {f"{mnv.chrom}:{mnv.pos}": mnv for mnv_lines in mnvs_to_keep for mnv in mnv_lines}
    
    # Create a dictionary where the keys are chromosome+position of the first MNV in each group and the values are the entire MNV group (list of VCFLine objects)
    mnv_groups = {f"{mnv[0].chrom}:{mnv[0].pos}": mnv for mnv in mnvs_to_keep}
    
    with gzip.open(phased_vcf_file, 'rt') as infile, pysam.BGZFile(output_vcf_file, 'w') as outfile:
        # Read and write the header lines
        header_lines = [line for line in infile if line.startswith('#')]
        new_header_lines = create_new_header_lines(header_lines, mnvs_to_keep)
        
        for line in header_lines:
            if not line.startswith('#CHROM'):
                outfile.write(line.encode('utf-8'))
            else:
                for new_header in new_header_lines:
                    outfile.write(new_header.encode('utf-8'))
                outfile.write(line.encode('utf-8'))
        
        # Variables to track the current MNV group and lines belonging to the group
        current_group = None
        group_lines = []
        
        # Re-open the file to process the non-header lines
        with gzip.open(phased_vcf_file, 'rt') as infile:
            for line in infile:
                if line.startswith('#'):
                    continue
                
                # Split the line into its component fields
                vcf_line = VCFLine(*line.split('\t'))
                chrom_pos = f"{vcf_line.chrom}:{vcf_line.pos}"
                
                if chrom_pos in mnv_positions:
                    # If the current line's position matches an MNV position
                    if current_group is None:
                        # If no current group is being processed, start a new group
                        current_group = mnv_groups[chrom_pos]
                    
                    if chrom_pos in [f"{vcf.chrom}:{vcf.pos}" for vcf in current_group]:
                        # If the current position is within the current MNV group, add it to the group_lines
                        group_lines.append(vcf_line)
                    
                    if len(group_lines) == len(current_group):
                        # If all lines in the current MNV group have been collected, merge them
                        merged_line = merge_vcf_lines(group_lines)
                        # Write the merged line to the output file
                        outfile.write(merged_line.encode('utf-8'))
                        # Reset the current group and group_lines
                        current_group = None
                        group_lines = []
                else:
                    # Write the non-MNV lines to the output file
                    outfile.write(line.encode('utf-8'))
    
    logging.info(f'Successfully wrote VCF containing MNVs to {output_vcf_file}')

def create_new_header_lines(header_lines, mnvs_to_keep) -> list[str]:
    '''Create the appropriate number of new header lines based on the maximum MNV length observed'''
    header_lines: list[str]
    mnvs_to_keep: list
    
    if not mnvs_to_keep:
        # If no MNVs are found, return an empty list (no new headers needed)
        return []

    mnv_lengths = [len(mnv) for mnv in mnvs_to_keep]
    max_mnv_length = max(mnv_lengths)
    
    new_header_lines = []
    
    for line in header_lines:
        if line.startswith('##INFO') or line.startswith('##FORMAT'):
            line_parts = line.split(',')
            for i in (range(1, max_mnv_length + 1)):
                new_start = line_parts[0] + f'_{i}'
                middle = ','.join(line_parts[1:-1])
                new_end = line_parts[-1].strip().replace('">', '') + f' (MNV allele {i} in series)' + '">\n'
                new_line = ','.join([new_start, middle, new_end])
                new_header_lines.append(new_line)
    
    return new_header_lines


def merge_vcf_lines(lines) -> str:
    '''Merge multiple VCF lines into a single line'''
    lines: list[str]
    
    if not lines:
        return ''
    
    # Convert lines into VCFLine objects
    vcf_lines = [VCFLine(*line) for line in lines]

    # Create merged MNV line using SNV variant lines
    chrom = vcf_lines[0].chrom # All variants have the same chromosome
    pos = vcf_lines[0].pos # Take the position of the first SNV in the MNV
    ids = [] # List to store IDs that will be merged later
    modified_ref = ''
    modified_alt = ''
    qual = vcf_lines[0].qual # We don't have QUAL values. Take the "." from the first line
    filters = set() # Set to store all FILTER fields observed in MNV
    
    for line in vcf_lines:
        # Concatenate SNV REF and ALTs to create new REF and ALT for the MNV
        modified_ref += line.ref
        modified_alt += line.alt
        ids.append(line.id)
        if ';' in line.filter:
            for filter in line.filter.split(';'):
                filters.add(filter)
        else:
            filters.add(line.filter)
    
    # Construct new ID field
    modified_id = ';'.join(ids)
    
    # Construct new FILTER field
    if len(filters) == 1:
        modified_filter = filters.pop()
    else:
        modified_filter = ';'.join(filters)
    
    # Construct new INFO field
    info_fields = [line.info for line in vcf_lines]
    modified_vcf_info_fields = []
    
    for index, info in enumerate(info_fields, start=1):
        # Split the INFO fields into key-value pairs
        key_value_pairs = info.split(';')
        # Modify each key by appending the suffix "_{index}"
        modified_pairs = [f"{pair.split('=')[0]}_{index}={pair.split('=')[1]}" for pair in key_value_pairs if '=' in pair]
        # Join the modified key-value pairs back into a string
        modified_info = ';'.join(modified_pairs)
        # Append the modified INFO field to the new list
        modified_vcf_info_fields.append(modified_info)
    
    modified_info = ';'.join(modified_vcf_info_fields) 
    
    # Extract FORMAT field keys and values from both NORMAL and TUMOUR
    format_keys = vcf_lines[0].format.split(':')
    normal_format_values = [vcf_line.normal for vcf_line in vcf_lines]
    tumour_format_values = [vcf_line.tumour.strip() for vcf_line in vcf_lines]
    
    # Determine the number of SNVs in the MNV
    num_snvs = len(vcf_lines)
    
    # Generate the modified FORMAT keys with the correct number of suffixes
    modified_keys = [f"{key}_{i+1}" for i in range(num_snvs) for key in format_keys]
    modified_format = ':'.join(modified_keys)
    
    # Iterate over NORMAL and TUMOUR format fields and populate modified format keys
    modified_normal_format_values = ':'.join(normal_format_values)
    modified_tumour_format_values = ':'.join(tumour_format_values)
    
    # Construct merged line using modified fields
    merged_line = '\t'.join([
        chrom,
        pos,
        modified_id,
        modified_ref,
        modified_alt,
        qual,
        modified_filter,
        modified_info,
        modified_format,
        modified_normal_format_values,
        modified_tumour_format_values
        
    ])
    
    merged_line = merged_line + '\n'
    
    return merged_line

def main():

    args = get_args()
    logging.info(f'Received command line arguments: {args}')
    mnvs_to_keep = check_mnv_phasing(args.phased_variants)
    output_file = str(args.phased_variants).replace('.vcf.gz', '.mnvs.vcf.gz')
    modify_vcf_lines(args.variants_file, mnvs_to_keep, output_file)


if __name__ == '__main__':
    main()
