#%%
# from functools import partial
# from random import choices, randint, randrange, random, sample
# from typing import List, Optional, Callable, Tuple
import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
import pandas as pd
# from collections import Counter
# from tqdm import tqdm
# import time
# from Bio.SeqUtils import MeltingTemp
# from Bio import SeqIO
from plotly import graph_objects as go
# import json
from imp import reload
# import argparse
from functools import reduce
from toast_amplicon import working_algo_gene_2in1 as w
from toast_amplicon import primer_selection as p_s
import os

#%%
def genome_size(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length

def rolling_sum_search(df, weight, window_size, genomic_pos):
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
    # print(pos)
    for x in pos:
        start = x
        in_range = [i for i in pos if i <= start+window_size]
        end = min(in_range, key=lambda x:abs(x-(start+window_size)))
        freq_sum = df[(df['genome_pos']>=start) & (df['genome_pos']<=end)][f'{weight}'].sum()
        rolling_sum.append(freq_sum)
    return rolling_sum
#%%
def place_amplicon_search(full_data, target_coverage, read_size, ref_size, output_path, graphic_output=False):
    full_data_cp = full_data.copy()
    read_size = read_size
    window_size = read_size
    read_number = 0
    # priorities = []
    graph_output = graphic_output # set to True to see the graph output
    weight_window_sum = rolling_sum_search(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())
    pos = full_data_cp['genome_pos'].unique()
    # covered_positions = pd.DataFrame(columns = ['sample_id','genome_pos','gene','change','freq', 'type','sublin','drtype','drugs','weight'])

    dtypes = {}
    # Iterate over each column in the DataFrame
    for column in full_data.columns:
        # Get the data type of the column
        dtype = str(full_data[column].dtype)
        # Add the column and its data type to the dictionary
        dtypes[column] = dtype
    covered_positions = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

    coverage_trace = [0]
    coverage = coverage_trace[-1]
    amplicon_number = 0
    amplicon_number_list = []
    coverage_list = []
    coverage_range = []
    print('Placing Amplicons...')
    while coverage < target_coverage*1:
    # while coverage < target_coverage*0.995:
        amplicon_number += 1
        amplicon_number_list.append(amplicon_number)
        start = pos[np.argmax(weight_window_sum)] # find the index of the max value in the rolling sum
        # in_range = [i for i in pos if i <= start+window_size] # find all values in the window
        # end = min(in_range, key=lambda x:abs(x-(start+window_size))) # find the closest value to the end of the window
        end = start+window_size
        if end >  ref_size:
            end =  ref_size

        coverage_range.append([start, end])
        covered_positions = covered_positions.reset_index(drop=True)
        
        full_data_cp.loc[(full_data_cp['genome_pos']>=start) & (full_data_cp['genome_pos']<=end), 'weight']= full_data_cp.loc[(full_data_cp['genome_pos']>=start) & (full_data_cp['genome_pos']<=end), 'weight'] = 0 # set the weight of the covered positions to 0
        
        df2 = full_data_cp[(full_data_cp['genome_pos']>= start) & (full_data_cp['genome_pos']<=end)].reset_index(drop=True)
        
        # print(df2)
        # covered_positions = covered_positions.dropna(axis=1, how='all')
        # df2 = df2.dropna(axis=1, how='all')
        covered_positions = pd.concat([covered_positions, df2])
        # covered_positions = covered_positions.drop_duplicates()
        # print(covered_positions['genome_pos'].unique().shape[0], full_data_cp['genome_pos'].unique().shape[0])
        # print(covered_positions['genome_pos'].unique().shape[0]/full_data_cp['genome_pos'].unique().shape[0])
        coverage_trace.append(covered_positions.shape[0]/full_data.shape[0])
        # coverage = coverage_trace[-1]
        coverage = covered_positions['genome_pos'].unique().shape[0]/full_data_cp['genome_pos'].unique().shape[0]
        
                # covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}}  # concise version of output
        # Generating sample data
        # Displaying the plot
        print(f'Amplicon#{amplicon_number}: SNP-coverage: {round(coverage,3)*100}%')

        # full_data_cp.loc[(full_data_cp['genome_pos']>=start) & (full_data_cp['genome_pos']<=end), 'weight'] = 0 # set the weight of the covered positions to 0
        coverage_list.append(coverage*100)
        weight_window_sum = rolling_sum_search(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist()) # recalculate the rolling sum
    
    print(f'**{amplicon_number} amplicons needed to cover {round(coverage,4)*100}% SNPs')

    output = pd.DataFrame({'Amplicon_number': amplicon_number_list, 'SNP_coverage(%)': coverage_list})
    
    os.makedirs(f'{output_path}/Amplicon_num', exist_ok=True)
    output.to_csv(f'{output_path}/Amplicon_num/amplicon_no_estimation_-{read_size}bps.csv', index=False)
    print(f'**Output saved to: {output_path}/Amplicon_num/amplicon_no_estimation_-{read_size}bps.csv')
    if graphic_output:
        x = np.linspace(0, len(coverage_trace), len(coverage_trace))
        y = [x*100 for x in coverage_trace]

        # Creating a figure
        fig = go.Figure(data=go.Scatter(x=x, y=y, mode='lines'))

        # Setting title and labels
        fig.update_layout(title=f'Number of {read_size}bps Amplicon needed for {target_coverage*100}% SNP coverage',
                        xaxis_title='No. of Amplicones',
                        yaxis_title='SNP Coverage (%)')
        # fig.add_vline(x=amplicon_number, line_width=3, line_dash="dash", line_color="green")

        # fig.show() 
        if output_path:
            os.makedirs(f'{output_path}/Amplicon_num', exist_ok=True)
            fig.write_image(f'{output_path}/Amplicon_num/coverage_trace_{target_coverage}-{read_size}.pdf')
            print(f'**Graphic output saved to: {output_path}/Amplicon_num/coverage_trace_{target_coverage}-{read_size}bps.pdf')
        else:
            print('No output path specified, graph not saved')
    # gene_coverage = covered_positions['gene'].value_counts()/full_data_cp['gene'].value_counts()
    
    return coverage_range

#%%
# def place_amplicon(full_data, read_number, read_size, primer_pool, accepted_primers, no_primer_, ref_genome, graphic_output=False, padding=150, output_path = '.'):
# def place_amplicon_spol(full_data, target_coverage, read_size, ref_genome, primer_pool, accepted_primers, no_primer_, padding = 150, graphic_output=False, check_snp=False, redesign_flag=False):
#     full_data_cp = full_data.copy()
#     read_size = read_size
#     window_size = read_size
#     read_number = 0
#     ref_size = genome_size(ref_genome)
#     # priorities = []
#     graph_output = graphic_output # set to True to see the graph output
#     weight_window_sum = rolling_sum_search(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())    
#     pos = full_data_cp['genome_pos'].unique()
#     covered_positions = {}
#     covered_ranges = []
    
#     dtypes = {
#     'genome_pos': 'float64',
#     'weight': 'float64'
#     }
#     covered_positions = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)
#     # covered_positions = pd.DataFrame(columns = ['genome_pos','weight'])
#     run = 0
#     coverage_trace = [0]
#     coverage = coverage_trace[-1]
#     amplicon_number = 0
#     print('Placing Amplicones for Spoligotypes...')
#     while coverage < target_coverage*100:
#         print(f'**Amplicon #{amplicon_number+1}')
#         print(f'Designing primers...for {read_size}bps Amplicons...with {padding}bps padding')
#         print(f'Current coverage: {round(coverage,1)}%')
#         if redesign_flag:
#             response = input("!!!There has been previous redesigning going on. Do you want to proceed with further amplicon placement? [Y/n]: ").strip().lower()
#             # Setting default response to 'yes' if the user enters nothing
#             if response == '' or response == 'y':
#                 pass
#             else:
#                 break
#             print('***')
#         amplicon_number += 1
        
#         # print(amplicon_number, coverage)
#         start_r = pos[np.argmax(weight_window_sum)] # find the index of the max value in the rolling sum
#         # in_range = [i for i in pos if i <= start+window_size] # find all values in the window
#         # end = min(in_range, key=lambda x:abs(x-(start+window_size))) # find the closest value to the end of the window
#         end_r = start_r +window_size
#         if end_r >  ref_size:
#             end_r =  ref_size-200
        
#         # covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start, 'End': end}, 'Markers':full_data_cp[['genome_pos','gene','sublin','drtype','drugs','weight']][start:end].sort_values(by=['weight']).to_dict('records')} 
#         seq_template, low_b, high_b = w.extraction_prep([start_r, end_r], ref_size = ref_size, ref_genome=ref_genome, padding=padding)
#         primer_pool, accepted_primers, no_primer = p_s.result_extraction(primer_pool, accepted_primers, seq_template, run+1, padding, ref_genome, high_b, low_b, read_size, full_data_cp, check_snp, freq_cutoff=50000)
#         # primer_pool, accepted_primers, no_primer = p_s.result_extraction(primer_pool, accepted_primers, seq_template, run+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, read_size = read_size, priority = full_data_cp, freq_cutoff=50000, check_snp=check_snp)
#         no_primer_.extend(no_primer)
#         # print(no_primer_)

#         # if no_primer[-1] == 'Redesigned':
#         #     redesign_flag = True
        
#         if accepted_primers.shape[0] != 0:
#             start_p, end_p = accepted_primers.iloc[accepted_primers.shape[0]-1][['pLeft_coord','pRight_coord']].values
#         else:
#             run= max(0,run-1)
#             print('No suitable primers found')
#             break
        
#         # print('start_p', 'end_p')
#         # print(start_p, end_p)
#         # c = full_data_cp.shape[0]
#         # full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p), 'weight'] = full_data_cp['weight'].min()/10/c  # set the weight of the covered positions smaller
#         full_data_cp.loc[(full_data_cp['genome_pos']>=int(start_p)) & (full_data_cp['genome_pos']<=int(end_p)), 'weight'] = 0 # set the weight of the covered positions smaller

#         if [start_p, end_p] in covered_ranges:
#             if pos[np.argmax(weight_window_sum)] == start_r:
#                 # c = full_data_cp[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r)].shape[0]
#                 # c = full_data_cp.shape[0]
                
#                 # full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'] = full_data_cp['weight'].min()/10/c # set the weight of the covered positions smaller
#                 full_data_cp.loc[(full_data_cp['genome_pos']>=start_r) & (full_data_cp['genome_pos']<=end_r), 'weight'] = 0 # set the weight of the covered positions smaller
#             # this is when problem comes, there is a difference in range coverage according to the design by weighted sum, however the actual range obtained from designed primers are dont cover the same range, hence the sae primers are repeatedly designed 
#                 # covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start_r, 'End': end_r}, 'Markers':full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)][['genome_pos','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
#                 covered_ranges.append([start_r, end_r])

#             accepted_primers = accepted_primers.iloc[:-1]
#             # print('***')
            
#         else:
#             run += 1
#             covered_ranges.append([start_p, end_p])

#             # covered_positions[f'Amplicon_{run+1}'] = {'Range':{'Start': start_p, 'End': end_p}, 'Markers':full_data_cp[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)][['genome_pos','weight']].sort_values(by=['weight']).to_dict('records')}# verbose version of output
#             print('***')

#             # print(full_data_cp.loc[(full_data_cp['genome_pos']>=start_p) & (full_data_cp['genome_pos']<=end_p)]['weight'].values)
#             # print('==============')
#         # print(weight_window_sum)

#         covered_positions = covered_positions.reset_index(drop=True)
#         df2 = full_data_cp[(full_data_cp['genome_pos']>= start_p) & (full_data_cp['genome_pos']<=end_p)].reset_index(drop=True)
#         covered_positions = pd.concat([covered_positions, df2])
#         # covered_positions = covered_positions.drop_duplicates()
#         # print('--------')
#         # print(covered_positions)
#         # print(full_data_cp)
#         # print(len(covered_positions['genome_pos'].unique()), full_data_cp.shape[0])
#         # print(len(covered_positions['genome_pos'].unique())/full_data_cp.shape[0])
#         coverage_trace.append(len(covered_positions['genome_pos'].unique())/full_data_cp.shape[0]*100)
#         coverage = coverage_trace[-1]
        
#         weight_window_sum = rolling_sum_search(full_data_cp, 'weight', window_size, full_data_cp['genome_pos'].tolist())    
    
#     print(f'{amplicon_number} amplicons are needed cover {round(coverage,1)}% SNPs')
#     # print(coverage_trace)
#     return covered_positions, covered_ranges, full_data_cp, primer_pool, accepted_primers, no_primer_
    