#version to fall back on

#%%
# from functools import partial
# from random import choices, randint, randrange, random, sample
# from typing import List, Optional, Callable, Tuple
import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
# from rich_argparse import ArgumentDefaultsRichHelpFormatter
import pandas as pd
# from collections import Counter
# from tqdm import tqdm
# from Bio.SeqUtils import MeltingTemp
# from Bio import SeqIO
from plotly import graph_objects as go
# import json
# from imp import reload
from toast_amplicon import primer_selection
# reload(primer_selection)
# import Amplicon_no
# reload(Amplicon_no)
# import argparse
# from functools import reduce
import os
# import plotting1
# reload(plotting1)

#%%
def value_counts_list(lst):
    """
    Computes the frequency count of unique elements in a list and returns a dictionary, sorted by frequency count in
    descending order.

    Args:
    - lst (list): List of elements

    Returns:
    - dict: Dictionary with unique elements as keys and their frequency count as values, sorted by frequency count
    in descending order
    """
    value_counts = {}
    for item in lst:
        if item in value_counts:
            value_counts[item] += 1
        else:
            value_counts[item] = 1
    sorted_value_counts = dict(sorted(value_counts.items(), key=lambda x: x[1], reverse=True))
    return sorted_value_counts

def print_full(x):
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', 2000)
    pd.set_option('display.float_format', '{:20,.2f}'.format)
    pd.set_option('display.max_colwidth', None)
    print(x)
    pd.reset_option('display.max_rows')
    pd.reset_option('display.max_columns')
    pd.reset_option('display.width')
    pd.reset_option('display.float_format')
    pd.reset_option('display.max_colwidth')


tb_drug_resistance_genes = {
    'gyrB': ['Levofloxacin'],
    'gyrA': ['Levofloxacin'],
    'mshA': ['Isoniazid'],
    'rpoB': ['Rifampicin'],
    'rpoC': ['Rifampicin'],
    'rpsL': ['Streptomycin'],
    'embR': ['Ethambutol'],
    'rrs': ['Kanamycin', 'Capreomycin', 'Amikacin', 'Streptomycin'],
    'fabG1': ['Isoniazid'],
    'inhA': ['Isoniazid'],
    'rpsA': ['Pyrazinamide'],
    'tlyA': ['Capreomycin'],
    'ndh': ['Isoniazid'],
    'katG': ['Isoniazid'],
    'pncA': ['Pyrazinamide'],
    'kasA': ['Isoniazid'],
    'eis': ['Kanamycin', 'Amikacin'],
    'ahpC': ['Isoniazid'],
    'rpoA': ['Rifampicin'],
    'panD': ['Pyrazinamide'],
    'embC': ['Ethambutol'],
    'embA': ['Ethambutol'],
    'embB': ['Ethambutol'],
    'ubiA': ['Ethambutol'],
    'gid': ['Streptomycin']
}


# full_data.loc[full_data['gene'].isin(tb_drug_resistance_genes.keys()), 'weight'] += 0.5
# full_data.to_csv('snp_priority.csv', index=False)

def rolling_sum(df, weight, window_size, genomic_pos):
    """
    Calculates the rolling sum of a list with a given window size.

    Parameters:
        lst (list): The list to calculate the rolling sum for.
        window_size (int): The size of the rolling window.

    Returns:
        list: A list containing the rolling sum values.
    """
    # Calculate the rolling sum using a list comprehension
    rolling_sum = []
    pos = np.unique(genomic_pos).tolist()
    for x in pos:
        start = x
        in_range = [i for i in pos if i <= start+window_size]
        end = min(in_range, key=lambda x:abs(x-(start+window_size)))
        # end = x + window_size
        freq_sum = df[(df['genome_pos']>=start) & (df['genome_pos']<=end)][str(weight)].sum()
        rolling_sum.append(freq_sum)
    return rolling_sum

def genome_size(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length

def extraction_prep(x, ref_size, ref_genome, padding=150):
    # padding = int(padding/2)
    low_b = x[0]
    high_b = x[1]
    if low_b <= padding:
        low_b = padding+1
    elif high_b >= ref_size-padding:
        high_b = ref_size-padding
    
    #     low_b-= 150
    elif (high_b - low_b) < padding*2+50:
        # print('======')
        high_b+= padding*3
    else:
        high_b+= padding
        low_b-= padding
    # seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=padding)
    seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=0, fasta_file=ref_genome) #set to 0 to avoid padding again
    # seq_template = primer_selection.extract_sequence_from_fasta(low_b, high_b, padding=int(padding/2), fasta_file=ref_genome) #set to 0 to avoid padding again
    # print(low_b, high_b)
    # print(seq_template)
    return seq_template, low_b, high_b

def nucleotide_to_iupac(nucleotides):
    iupac_codes = {
        'R': {'A', 'G'},
        'Y': {'C', 'T'},
        'S': {'G', 'C'},
        'W': {'A', 'T'},
        'K': {'G', 'T'},
        'M': {'A', 'C'},
        'B': {'C', 'G', 'T'},
        'D': {'A', 'G', 'T'},
        'H': {'A', 'C', 'T'},
        'V': {'A', 'C', 'G'},
        'N': {'A', 'C', 'G', 'T'}
    }

    # Find the IUPAC code that matches the set of nucleotides
    nucleotide_set = set(nucleotides.upper())
    for code, bases in iupac_codes.items():
        if nucleotide_set == bases:
            return code
    return None
# Example usage
# print(nucleotide_to_iupac("AG"))  # Should return 'R'
# print(nucleotide_to_iupac("CT"))  # Should return 'Y'
#%%
# ideal_range = []
#place_amplicon function
def place_amplicon(full_data, read_number, read_size, primer_pool, accepted_primers, no_primer_, ref_genome, global_args, graphic_output=False, padding=150, output_path = '.', check_snp = True, start_end_r=None):
    # print( read_number, read_size, primer_pool, accepted_primers, no_primer_, ref_genome, graphic_output, padding, output_path)
    window_size = read_size
    run = 0
    full_data_cp = full_data.copy()
    # priorities = []
    # print('Calculating weight sum...')
    weight_window_sum = rolling_sum(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())
    pos = full_data_cp['genome_pos'].unique()
    covered_positions = {}
    covered_ranges = []
    designed_ranges = []
    reduce_amplicon = 0
    
    # primer design storage
    # accepted_primers = pd.DataFrame(columns=['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC', 'pLeft_Sequences', 'pLeft_EndStability','pRight_ID', 'pRight_coord', 'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences', 'pRight_EndStability', 'Penalty', 'Product_size'])
    # primer_pool = []
    # no_primer_ = no_primer_
    print('Placing Amplicons...')
    # for run in tqdm(range(0,read_number)):
    while run < read_number:
        print(f'**Amplicon #{run+1}')

        if graphic_output == True:
            op = f'{output_path}/Running_graphs'
            os.makedirs(op, exist_ok=True) #output path
            # print(np.argmax(weight_window_sum))

            trace = go.Scatter(
            x=[item*window_size for item in list(range(1, len(weight_window_sum) + 1))],
            y=weight_window_sum,
            mode='lines',
            line=dict(color='blue'),
            fill='tozeroy',
            fillcolor='rgba(0, 0, 255, 0.3)')
            # Create the layout§
            layout = go.Layout(
                title=f'Weight sum calculation by {window_size}bps sliding windows - #{run+1}',
                xaxis=dict(title='Sliding Window Genomic Position(bps)'),
                yaxis=dict(title='Window Weight Sum'),
                shapes=[
                # Add a vertical line at x=8
                dict(
                    type='line',
                    x0=np.argmax(weight_window_sum)*window_size,
                    x1=np.argmax(weight_window_sum)*window_size+window_size,
                    y0=0,
                    y1=max(weight_window_sum),
                    line=dict(color='rgba(255, 0, 0, 0.5)', width=20)
                )
            ]
            )
            # Create the figure
            fig = go.Figure(data=[trace], layout=layout)
            # Display the plot
            # fig.show()
            fig.write_image(f"{op}/{window_size}bps_Amplicon-#{run+1}_df_size{full_data.shape[0]}.png")
            print(f'**Graphic output saved to: {op}')
        if start_end_r is None:
            start_r = pos[np.argmax(weight_window_sum)] # find the index of the max value in the rolling sum
            # in_range = [i for i in pos if i <= start+window_size] # find all values in the window
            # end = min(in_range, key=lambda x:abs(x-(start+window_size))) # find the closest value to the end of the window
            # print(start_r)
            end_r = start_r+window_size
            start_r = start_r
            end_r = end_r
                        
            
            _snps = full_data_cp[(full_data_cp['genome_pos'] >= start_r) & (full_data_cp['genome_pos'] <= end_r)]
            _snp_pos = _snps[_snps['weight'] > 0]['genome_pos'].tolist()
            
            if (len(_snp_pos)>0) and ((end_r - max(_snp_pos))/(min(_snp_pos) - start_r+10) > 10): # roughly to the the important snps to move to the middle of the amplicon instead of just being at hte start of the amplicon
                # print('-----------adjusted')
                mid = int(min(_snp_pos) + (max(_snp_pos) - min(_snp_pos))/2)
                _len = int((end_r - start_r)/2)
                start_r = mid - _len
                end_r = mid + _len
                        
        else:
            start_r, end_r = start_end_r[0], start_end_r[1]
                
        print(f'Designing primers...for {read_size}bps Amplicons...with {padding}bps padding for genomic region {start_r}-{end_r}')

        #print()
        if end_r > genome_size(ref_genome):
            end_r = genome_size(ref_genome)-padding
        designed_ranges.append([start_r, end_r])

        # ideal_range.append([start_r, end_r])

        seq_template, low_b, high_b = extraction_prep([start_r, end_r], ref_size = genome_size(ref_genome), ref_genome=ref_genome, padding=padding)
        try:
            primer_pool, accepted_primers, no_primer = primer_selection.result_extraction(primer_pool, accepted_primers, seq_template, run+1, padding, ref_genome, high_b, low_b, read_size, full_data_cp, check_snp, global_args, freq_cutoff=50000)
        except TypeError as e:
            # print("!!! No primer designed")
            print(f"Error details: {e}")
            print("!!! Try increasing the padding size")
            import sys
            sys.exit(1)
        try:
            no_primer = [no_primer[-1]]
        except TypeError as e:
            print(f"Error details: {e}")
            print("!!! Try increasing the padding size")
            import sys
            sys.exit(1)
            
        no_primer_.extend(no_primer)
        
        if accepted_primers.shape[0] != 0:
            start_p, end_p = accepted_primers.iloc[accepted_primers.shape[0]-1][['pLeft_coord','pRight_coord']].values
        else:
            run= max(0,run-1)
            print('! No suitable primers found')
            break
        # print('---WindowMin', min(weight_window_sum))
        # print('---WindowMean', np.mean(weight_window_sum))
        # print('---Min', full_data_cp['weight'].min())
        # print('---Mean', full_data_cp['weight'].mean())
        # c = full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)].shape[0]
        # c = full_data_cp.shape[0]
        # full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = full_data_cp['weight'].min()/10/c  # set the weight of the covered positions smaller
        full_data_cp['weight'] = full_data_cp['weight'].astype(float)
        full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'] = full_data_cp['weight'].min()/2 # set the weight of the covered positions smaller
        full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = -0.05 # set the weight of the covered positions smaller
        # full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = 0 # set the weight of the covered positions smaller
        # print(start_p, end_p)
        # full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight']/100 # set the weight of the covered positions smaller
        # print(full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)]['weight'])
        # print('min',full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)]['genome_pos'].min())
        # # print(full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'].sum())
        # print(f'Cover ranges {run+1}/{read_number}: {[start_p, end_p]}==========')
# full_data_cp[(full_data_cp['genome_pos']>=1673440) & (full_data_cp['genome_pos']<=1674487)]
# full_data_cp[(full_data_cp['genome_pos']>=1673373) & (full_data_cp['genome_pos']<=1674373)]
        if [start_p, end_p] in covered_ranges:
            # if pos[np.argmax(weight_window_sum)] == start_r-50:
                # print('!this happens')
                # c = full_data_cp[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r)].shape[0]
                # c = full_data_cp.shape[0]
            print( full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'])
                # full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'] = full_data_cp['weight'].min()/10/c # set the weight of the covered positions smaller
            full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'] -= 0.05 # set the weight of the covered positions smaller
            full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] -= 0.05
            # else:
                # print('!this happens111')
                # full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = 0 # set the weight of the covered positions smaller
            # this is when problem comes, there is a difference in range coverage according to the design by weighted sum, however the actual range obtained from designed primers are dont cover the same range, hence the sae primers are repeatedly designed 
            print('Already covered, consider reducing amplicon number...Finding alternative sequences...')
            print( full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'])
            
            reduce_amplicon += 1
            accepted_primers = accepted_primers.iloc[:-1]
            # print('***')
            designed_ranges = designed_ranges[:-1]
            no_primer_ = no_primer_[:-1]
        else:
            run += 1
            covered_ranges.append([start_p, end_p])

            covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start_p, 'End': end_p}, 'Markers':full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)][['genome_pos','gene','sublin','drtype','drugs','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
            print('***')
            
            # print(full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)]['weight'].values)
            # print('==============')
        # print(weight_window_sum)
        weight_window_sum = rolling_sum(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())
    if reduce_amplicon > 0:
        pass
        # print(f'Consider reducing number of amplicons by: {reduce_amplicon}')

    # print('====================')

    # return covered_positions, covered_ranges, full_data_cp, primer_pool, accepted_primers, no_primer_
    return covered_positions, designed_ranges, full_data_cp, primer_pool, accepted_primers, no_primer_

def modify_primer_name(primer, amplicon_type, L_R):
    if 'User' in primer:
        return primer
    else:
        prefix = 'sp' if 'Gene_specific' in amplicon_type else 'ns' if 'Non_specific' in amplicon_type else 'spol' if 'Spoligotype' in amplicon_type else ''
        parts = primer.split('-')
        if len(parts) > 1:
            parts[-1] = L_R + str(int(parts[-1][1:]) + 1)  # Increment the number
        parts.insert(1, prefix)
        return '-'.join(parts)