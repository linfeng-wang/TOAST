#here I would try to gather al paddings together and make new genome
# the versions to fall back on
#%%
# from functools import partial
# from random import choices, randint, randrange, random, sample
# from typing import List, Optional, Callable, Tuple
# from importlib import reload

import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
import pandas as pd
# from collections import Counter
# from tqdm import tqdm
# import time
# from Bio.SeqUtils import MeltingTemp
# from Bio.Blast import NCBIWWW, NCBIXML
# import primer3
from Bio.Seq import Seq
# from Bio.SeqUtils import MeltingTemp
from primer3 import calc_heterodimer
from primer3 import bindings
from Bio import SeqIO
from toast_amplicon import Thermo
import json
from functools import lru_cache

# reload(Thermo)
# %%
def calculate_gc_content(sequence):
    """
    Calculate the percentage of G and C nucleotides in a DNA sequence.

    Args:
        sequence (str): DNA sequence string.

    Returns:
        float: Percentage of G and C nucleotides in the sequence.
    """
    gc_count = 0
    total_count = 0

    for nucleotide in sequence:
        if nucleotide.upper() in ['G', 'C']:
            gc_count += 1
        if nucleotide.upper() in ['A', 'T', 'G', 'C']:
            total_count += 1

    gc_percentage = (gc_count / total_count) * 100
    return gc_percentage

def extract_sequence_from_fasta(start_pos, end_pos, padding = 150, fasta_file= '../db/MTB-h37rv_asm19595v2-eg18.fa', sequence_id='Chromosome'):
    """
    Extracts a subsequence from a FASTA file based on the given sequence ID, start position, and end position.
    """
    padding = int(padding)
    # Iterate over the sequences in the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Check if the current sequence ID matches the desired sequence ID
        
        if record.id == sequence_id:
            # Extract the subsequence based on the start and end positions
            subsequence = record.seq[start_pos-padding:end_pos+padding]
            return str(subsequence)  # Return the subsequence as a string
    # If the sequence ID is not found, return None
    return None

def genome_size(fasta_file):
    total_length = 0
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                total_length += len(line.strip())
    return total_length
#%%
#original
# def complement_sequence(dna_sequence):
#     trans = str.maketrans('ATCGRYSWKMNBVDH-.', 'TAGCYRWSMKNVBHD-.')
#     return dna_sequence.upper().translate(trans)

# faster
@lru_cache(maxsize=None)
def complement_sequence(seq):
    return seq.translate(str.maketrans('ATCGRYSWKMNBVDH-.', 'TAGCYRWSMKNVBHD-.'))

# print(complement_sequence("ATGCGTA"))
# Example usage:
# dna_sequence = "ATGCGTA"
# complement_sequence = complement_sequence(dna_sequence)
# print(complement_sequence)  # Output: TACGCAT
#%%
#original
# def reverse_complement_sequence(seq):
#     complement = {"A": "T", "T": "A", "C": "G", "G": "C", 'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'N': 'N', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', '-':'-','.': '.'}
#     reverse_seq = seq[::-1]
#     return "".join(complement[base] for base in reverse_seq)

# faster
@lru_cache(maxsize=None)
def reverse_complement_sequence(seq):
    complement = str.maketrans('ATCGRYSWKMNBVDH-.', 'TAGCYRWSMKNVBHD-.')
    return seq.translate(complement)[::-1]


def updown_stream_primer_range(start_pos, end_pos, dis_range=0):
    up_stream = complement_sequence(extract_sequence_from_fasta(start_pos-dis_range, start_pos))
    down_stream = reverse_complement_sequence(extract_sequence_from_fasta(end_pos, end_pos+dis_range))
    return up_stream, down_stream

# %%
# def check_heterodimer(primer1, primer2):
#     # Calculate melting temperature (Tm) for each primer
#     # tm1 = MeltingTemp.Tm_NN(primer1)
#     # tm2 = MeltingTemp.Tm_NN(primer2)

#     # Check for heterodimer formation between the two primers
#     heterodimer = calc_heterodimer(primer1, primer2)
#     # Print the results
#     # print("Primer 1 Tm:", tm1)
#     # print("Primer 2 Tm:", tm2)
#     # print("Heterodimer:", heterodimer.structure_found)
#     return heterodimer.structure_found
# Example usage
# primer1 = "AGTCATCGATCGATCGATCG"
# primer2 = "CGATCGATCGATCGATCGAT"
# check_heterodimer(primer1, primer2)

def check_heterodimer(seq1, seq2, global_arg = None, length=5):
    """Checks if the last 'length' bases of seq1 are found anywhere in seq2."""
    # Get the last 'length' bases of seq1
    last_bases_seq1 = seq1[-length:]
    last_bases_seq1_com = complement_sequence(last_bases_seq1)
    last_bases_seq1_rev_com = last_bases_seq1_com[::-1]
    # Check if these bases are found anywhere in seq2
    is_binding_com = last_bases_seq1_com in seq2
    is_binding_rev_com = last_bases_seq1_rev_com in seq2
    is_binding = is_binding_com or is_binding_rev_com
    
    if is_binding == False:
        last_bases_seq1 = seq2[-length:]
        last_bases_seq1_com = complement_sequence(last_bases_seq1)
        last_bases_seq1_rev_com = last_bases_seq1_com[::-1]
        # Check if these bases are found anywhere in seq2
        is_binding_com = last_bases_seq1_com in seq2
        is_binding_rev_com = last_bases_seq1_rev_com in seq2
        is_binding = is_binding_com or is_binding_rev_com
    
    if global_arg == None:
        primer3_hetero = calc_heterodimer(seq1, seq2,  mv_conc=50, dv_conc=1.5, dntp_conc=0.8, dna_conc=50, temp_c=37, max_loop=30, output_structure=False)
    else:
        primer3_hetero = calc_heterodimer(seq1, seq2,  mv_conc=global_arg['PRIMER_SALT_MONOVALENT'], dv_conc=global_arg["PRIMER_SALT_DIVALENT"], dntp_conc=global_arg["PRIMER_DNTP_CONC"], dna_conc=global_arg['PRIMER_DNA_CONC'], temp_c=global_arg["PRIMER_ANNEALING_TEMP"])
        
    primer3_hetero_result = primer3_hetero.tm >= 50 and primer3_hetero.dg <= -9000 
    
    is_binding = is_binding or primer3_hetero_result
    return is_binding
#%%
# test = extract_sequence_from_fasta(1000, 2300)

def find_sequence_location(query_seq, fasta_file):
    ref_genome = SeqIO.read(fasta_file, "fasta")
    position = ref_genome.seq.find(query_seq)
    if position != -1:
        return [position, position+len(query_seq)]
    else:
        print('!!!Primer not found in the reference genome')
        return 0

# pLeft_coord = []
# pLeft_coord.append(find_sequence_location('CATCGCACGTCGTCTTTCCG', ref_genome)[0])

# find_sequence_location('CATCGCACGTCGTCTTTCCG')[0]
# find_sequence_location(complement_sequence('CATCGCACGTCGTCTTTCCG'))
#%%
def simplified_tm(seq):
    # Simplified melting temperature calculation based on base pair count
    return (seq.count('A') + seq.count('T')) * 2 + (seq.count('C') + seq.count('G')) * 4

# def simplified_dg(seq):
#     # Simplified Gibbs free energy calculation based on the nearest-neighbor model (not accurate)
#     nn_params = {
#         'AA': -1.0, 'TT': -1.0,
#         'AT': -0.9, 'TA': -0.9,
#         'CA': -1.7, 'TG': -1.7,
#         'GT': -1.5, 'AC': -1.5,
#         'CT': -1.6, 'AG': -1.6,
#         'GA': -1.5, 'TC': -1.5,
#         'CG': -2.8, 'GC': -2.3,
#         'GG': -1.9, 'CC': -1.9
#     }
#     dg = 0
#     for i in range(len(seq) - 1):
#         dinucleotide = seq[i:i+2]
#         dg += nn_params.get(dinucleotide, 0)
#     return dg
#original
# def calculate_similarity(seq1, seq2):
#     matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
#     return matches / len(seq1) * 100

## faster
def calculate_similarity(seq1, seq2):
    a = np.frombuffer(seq1.encode(), dtype='S1')
    b = np.frombuffer(seq2.encode(), dtype='S1')
    return (a == b).sum() / len(a) * 100
# def check_last_n_base_pairs(seq1, seq2, n):
#     """
#     Check if the last n base pairs of two sequences are the same.

#     Parameters:
#     - seq1: The first sequence (string)
#     - seq2: The second sequence (string)
#     - n: The number of base pairs to compare from the end (integer)

#     Returns:
#     - bool: True if the last n base pairs are the same, otherwise False
#     """
    
#     # Check if either sequence is shorter than n
#     if len(seq1) < n or len(seq2) < n:
#         return False
    
#     # Get the last n base pairs from each sequence
#     last_n_seq1 = seq1[-n:]
#     last_n_seq2 = seq2[-n:]
    
#     # Compare the last n base pairs
#     return last_n_seq1 == last_n_seq2

# Re-running the code after the execution state reset

# original
# def check_last_n_base_pairs(seq1, seq2, n=5, t=2):
#     """
#     Count the number of mismatches in the last n base pairs of two sequences.

#     Parameters:
#     - seq1: The first sequence (string)
#     - seq2: The second sequence (string)
#     - n: The number of base pairs to compare from the end (integer)

#     Returns:
#     - int: The number of mismatches in the last n base pairs
#     """
    
#     # Initialize mismatch count
#     mismatch_count = 0
    
#     # Check if either sequence is shorter than n
#     if len(seq1) < n or len(seq2) < n:
#         return "One or both sequences are shorter than n."
    
#     # Get the last n base pairs from each sequence
#     last_n_seq1 = seq1[-n:]
#     last_n_seq2 = seq2[-n:]
    
#     # Count mismatches in the last n base pairs
#     for base1, base2 in zip(last_n_seq1, last_n_seq2):
#         if base1 != base2:
#             mismatch_count += 1
            
#     return mismatch_count<2 # returns true if mismatch is less than 2 -> alternative binding - polymerases can go on

# faster
def check_last_n_base_pairs(seq1, seq2, n=5, t=2):
    # Return True if fewer than t mismatches in the last n bases
    mismatches = sum(a != b for a, b in zip(seq1[-n:], seq2[-n:]))
    return mismatches < t
# # Test the function
# seq1 = "ATCGATCG"
# seq2 = "ATCGATCA"
# n = 5

# count_mismatches_in_last_n_base_pairs(seq1, seq2, n)

# Example usage
# result = check_last_n_base_pairs("ATCGGA", "TTTGGA", 2)
# print("Are the last 2 base pairs the same?", result)  # Output should be True

# result = check_last_n_base_pairs("ATCGGA", "TTTGTA", 2)
# print("Are the last 2 base pairs the same?", result)  # Output should be False

#%% origional function
# def has_multiple_binding_sites(sequence, genome, similarity_threshold=90, min_tm=50, max_dg=-10):
#     class TmVar:
#         DivalentSaltConc = 2
#         MonovalentSaltConc = 10
#         dNTPConc = 0.2 
#         OligoConc = 0.5

#     seq_len = len(sequence)
#     genome_len = len(genome)
#     count = 0
#     # flocations = []
#     # rlocations = []
#     n = 5
#     t = 2
#     for i in range(genome_len - seq_len + 1):
#         within_tm_threshold, within_dg_threshold = False, False
#         subseq = genome[i:i + seq_len]
#         reverse_comp_subseq = reverse_complement_sequence(subseq)
#         if calculate_similarity(subseq, sequence)>similarity_threshold:
#             if check_last_n_base_pairs(subseq, sequence, n, t):
#                 # flocations.append(i)
#                 # print('>>>Calculated_similarity')
#                 within_tm_threshold, within_dg_threshold = False, False
#                 (dg, tm) = Thermo.simplified_tm(subseq, complement_sequence(subseq), TmVar)
                
#                 within_tm_threshold = tm >= min_tm
#                 within_dg_threshold = dg['Gibbs'] <= max_dg
#         elif calculate_similarity(reverse_comp_subseq, sequence)>similarity_threshold:
#             if check_last_n_base_pairs(reverse_comp_subseq, sequence, n, t):
#                 # print('====C')
                
#                 # rlocations.append(i)
#                 # print('>>>Calculated_similarity')
#                 within_tm_threshold, within_dg_threshold = False, False
#                 (dg, tm) = Thermo.simplified_tm(reverse_comp_subseq, complement_sequence(reverse_comp_subseq), TmVar)

#                 within_tm_threshold = tm >= min_tm
#                 within_dg_threshold = dg['Gibbs'] <= max_dg
#         if within_tm_threshold & within_dg_threshold: # the alternative binding site is within the threshold
#             # print('====D')
            
#             # print(genome[i:i + seq_len])
#             count += 1
            
#             if count > 1:
#                 # print('>>>within_tm_threshold & within_dg_threshold')
#                 # print("#####")
#                 # print(flocations)
#                 # print(rlocations)
#                 # print('====E')
#                 return True  # Early exit if more than one binding site is found
#     # print(flocations)
#     # print(rlocations)
#     return False  # Return False if only one or no binding sites are found


# # Usage:
# primer = 'CGAACTCGAGGCTGCCTACT'
# # primer = 'GCTCGTCCATGTCCCACCAT'
# genome = primer_selection.extract_sequence_from_fasta(0, genome_size(ref_genome),0)
# result = has_multiple_binding_sites(primer, genome, 81)
# print(result)  # Output: True or False

# new and faster function
def has_multiple_binding_sites(sequence, genome, similarity_threshold=90, min_tm=50, max_dg=-10):
    from functools import lru_cache

    class TmVar:
        DivalentSaltConc = 2
        MonovalentSaltConc = 10
        dNTPConc = 0.2 
        OligoConc = 0.5

    @lru_cache(maxsize=None)
    def rc(seq):
        return reverse_complement_sequence(seq)

    @lru_cache(maxsize=None)
    def comp(seq):
        return complement_sequence(seq)

    seq_len = len(sequence)
    genome_len = len(genome)
    count = 0
    n = 5
    t = 2

    for i in range(genome_len - seq_len + 1):
        subseq = genome[i:i + seq_len]
        reverse_subseq = rc(subseq)

        if subseq == sequence or reverse_subseq == sequence:
            continue
        # Precompute both orientations
        forward_sim = calculate_similarity(subseq, sequence)
        reverse_sim = calculate_similarity(reverse_subseq, sequence)

        if forward_sim > similarity_threshold:
            if not check_last_n_base_pairs(subseq, sequence, n, t):
                continue
            dg, tm = Thermo.simplified_tm(subseq, comp(subseq), TmVar)
        elif reverse_sim > similarity_threshold:
            if not check_last_n_base_pairs(reverse_subseq, sequence, n, t):
                continue
            dg, tm = Thermo.simplified_tm(reverse_subseq, comp(reverse_subseq), TmVar)
        else:
            continue  # Skip this window entirely

        # Thermodynamic thresholds
        if tm >= min_tm and dg['Gibbs'] <= max_dg:
            count += 1
            if count > 1:
                return True  # Found more than one binding site

    return False  # Zero or one valid site



# Check for material binding sites
def check_for_snp(seq, priority, plength, forward, ref_genome, cut_off=50000, spol = False):
    if not spol:
        full_gene = priority.groupby(['genome_pos', 'gene', 'change']).agg({'freq': 'sum'}).reset_index() # this is to group by the same position and sum up the freq
    
    if forward == True:
        ploc = find_sequence_location(seq, ref_genome)[0]
    else:    
        ploc = find_sequence_location(reverse_complement_sequence(seq), ref_genome)[0]
    # print('ploc',ploc)
    snp_freq = full_gene[(full_gene['genome_pos']>= ploc) & (full_gene['genome_pos']<=ploc+plength)]['freq']
    # print('snp_freq',snp_freq)
    # print('full_gene',full_gene[(full_gene['genome_pos']>= ploc) & (full_gene['genome_pos']<=ploc+plength)])
    # print(snp_freq.max())
    if snp_freq.max() > cut_off:
        return True
    else:
        return False
    
# %%
def result_extraction(primer_pool, accepted_primers, sequence, seq, padding, ref_genome, high_b, low_b, read_size, priority, check_snp, global_args, freq_cutoff=50000, _is_recursive='none'):
    # print([len(sequence)-50,len(sequence)+50])
    # print(len(sequence))
    # size_range = f'{int(len(sequence)-padding*1.3)}-{int(len(sequence)-padding*1)}'
    # print('---padding:',padding)
    # print('---sequence:',len(sequence))
    # padding = int(padding/2)
    size_range = f'{len(sequence)-padding*2}-{len(sequence)}'
    with open(global_args, 'r') as file:
        global_args_dict = json.load(file)
        global_args_dict['PRIMER_PRODUCT_SIZE_RANGE'] = size_range
        
    genome = extract_sequence_from_fasta(0, genome_size(ref_genome),padding=0, fasta_file=ref_genome)
    no_primer = []
    ok_region_list = [0, padding,len(sequence)-padding,padding]
    # size_range = f'{len(sequence)-350}-{len(sequence)-250}'
    # print('size_range:',size_range)
    # print('SEQUENCE_INCLUDED_REGION:', [padding-10,len(sequence)-padding+10],)
    try:
        results = bindings.design_primers(
            seq_args={
                'SEQUENCE_ID': 'Amplicon',
                'SEQUENCE_TEMPLATE': sequence,
                # 'SEQUENCE_INCLUDED_REGION': [padding-20,len(sequence)-padding+20],
                # 'SEQUENCE_INCLUDED_REGION': [padding,len(sequence)-(padding*2)],
                # 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': ok_region_list
                'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': f'{ok_region_list[0]},{ok_region_list[1]},{ok_region_list[2]},{ok_region_list[3]}',

                # 'SEQUENCE_INCLUDED_REGION': [(0,len(sequence)),],
                # 'SEQUENCE_INCLUDED_REGION': [(0,padding),(len(sequence)-padding,len(sequence))],
                # 'SEQUENCE_EXCLUDED_REGION':[(padding,len(sequence)-padding)]
                'SEQUENCE_TARGET': [padding,len(sequence)-padding*2]
            },
            global_args=global_args_dict)
    except:
        print('!!!Primer extraction error, trying with increased padding x2')
        padding = int(padding*2)
        size_range = f'{len(sequence)-padding*2}-{len(sequence)}'
        with open(global_args, 'r') as file:
            global_args_dict = json.load(file)
            global_args_dict['PRIMER_PRODUCT_SIZE_RANGE'] = size_range
            
        genome = extract_sequence_from_fasta(0, genome_size(ref_genome),padding=0, fasta_file=ref_genome)
        no_primer = []
        ok_region_list = [0, padding,len(sequence)-padding,padding]
        # size_range = f'{len(sequence)-350}-{len(sequence)-250}'
        # print('size_range:',size_range)
        # print('SEQUENCE_INCLUDED_REGION:', [padding-10,len(sequence)-padding+10],)
        results = bindings.design_primers(
            seq_args={
                'SEQUENCE_ID': 'Amplicon',
                'SEQUENCE_TEMPLATE': sequence,
                # 'SEQUENCE_INCLUDED_REGION': [padding-20,len(sequence)-padding+20],
                # 'SEQUENCE_INCLUDED_REGION': [padding,len(sequence)-(padding*2)],
                # 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': ok_region_list
                'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': f'{ok_region_list[0]},{ok_region_list[1]},{ok_region_list[2]},{ok_region_list[3]}',

                # 'SEQUENCE_INCLUDED_REGION': [(0,len(sequence)),],
                # 'SEQUENCE_INCLUDED_REGION': [(0,padding),(len(sequence)-padding,len(sequence))],
                # 'SEQUENCE_EXCLUDED_REGION':[(padding,len(sequence)-padding)]
                'SEQUENCE_TARGET': [padding,len(sequence)-padding*2]
            },
            global_args=global_args_dict)

    # print(results)
    
    pLeft_ID = []
    pLeft_coord = []
    pLeft_length = []
    pLeft_Tm = []
    pLeft_GC = []
    pLeft_Sequences = []
    pLeft_EndStability = []

    pRight_ID = []
    pRight_coord = []
    pRight_length = []
    pRight_Tm = []
    pRight_GC = []
    pRight_Sequences = []
    pRight_EndStability = []

    Penalty = []
    Product_size = []
    if len(results['PRIMER_PAIR']) == 0:
        print('!!!No primer designed')
        print('!!!Try changing the padding size')
        return 0
    
    for i, primer_num in enumerate(results['PRIMER_PAIR']):
        Product_size.append(primer_num['PRODUCT_SIZE'])
        Penalty.append(np.round(primer_num['PENALTY'],2))
    for i, primer_num in enumerate(results['PRIMER_LEFT']):
        pLeft_ID.append(f'P{seq}-L{i}')
        # pLeft_coord.append(primer_num['COORDS'][0]+low_b)
        # print(primer_num['SEQUENCE'])
        pLeft_coord.append(find_sequence_location(primer_num['SEQUENCE'], ref_genome)[0])
        
        pLeft_length.append(primer_num['COORDS'][1])
        pLeft_Tm.append(np.round(primer_num['TM'],2))
        pLeft_GC.append(primer_num['GC_PERCENT'])
        pLeft_Sequences.append(primer_num['SEQUENCE'])
        pLeft_EndStability.append(primer_num['END_STABILITY'])
        
    for i, primer_num in enumerate(results['PRIMER_RIGHT']):
        pRight_ID.append(f'P{seq}-R{i}')
        # pRight_coord.append(primer_num['COORDS'][0]+low_b)
        pRight_coord.append(find_sequence_location(reverse_complement_sequence(primer_num['SEQUENCE']),ref_genome)[1])
        pRight_length.append(primer_num['COORDS'][1])
        pRight_Tm.append(np.round(primer_num['TM'],2))
        pRight_GC.append(primer_num['GC_PERCENT'])
        # pRight_Sequences.append(reverse_complement_sequence(primer_num['SEQUENCE']))
        pRight_Sequences.append(primer_num['SEQUENCE'])
        pRight_EndStability.append(primer_num['END_STABILITY'])

    df = pd.DataFrame({'pLeft_ID':pLeft_ID, 'pLeft_coord':pLeft_coord, 'pLeft_length':pLeft_length, 'pLeft_Tm':pLeft_Tm, 'pLeft_GC':pLeft_GC, 'pLeft_Sequences':pLeft_Sequences, 'pLeft_EndStability':pLeft_EndStability, 
                    'pRight_ID':pRight_ID, 'pRight_coord':pRight_coord, 'pRight_length':pRight_length, 'pRight_Tm':pRight_Tm, 'pRight_GC':pRight_GC, 'pRight_Sequences':pRight_Sequences, 'pRight_EndStability':pRight_EndStability, 
                    'Penalty':Penalty, 'Product_size':Product_size})
    # print(df)
    # print(df[['pLeft_coord','pRight_coord','Product_size','pLeft_Sequences','pRight_Sequences']])
    # print('original_range:',low_b, high_b)
    # tm_params, dg_params = precompute_dinucleotide_params(genome)
    
    if len(primer_pool) == 0:
        print(f'{df.shape[0]-1} primers designed')
        for i, row in df.iterrows():
            # print(i, df.shape[0]-1)
            left_ok = True
            right_ok = True
            too_far = False
            #runs first to avoid meaningless running of has_multiple_binding_sites which take a long time
            if check_snp:
                if check_for_snp(row['pLeft_Sequences'], priority, row['pLeft_length'], True, ref_genome, cut_off=freq_cutoff) or check_for_snp(row['pRight_Sequences'], priority, row['pRight_length'], False, ref_genome, cut_off=freq_cutoff):
                    print(f'Primer pair #{i+1} binding position has SNP')
                    continue
            # print('low_b',low_b, 'row[pLeft_coord]',row['pLeft_coord'])
            # print('high_b',high_b, 'row[pRight_coord]',row['pRight_coord'])
            if abs(low_b - row['pLeft_coord']) > read_size/2:
            # if abs(low_b - row['pLeft_coord']) > read_size/2 or ~(low_b - row['pLeft_coord'])>0:
                # print('low')
                # print(low_b, row['pLeft_coord'])
                # print(abs(low_b - row['pLeft_coord']))
                left_ok = False
                too_far = True
            if abs(high_b - row['pRight_coord']) > read_size/2:
            # if abs(high_b - row['pRight_coord']) > read_size/2 or ~(high_b - row['pRight_coord'])<0:
                # print('high')
                # print(high_b, row['pRight_coord'])
                # print(abs(high_b - row['pRight_coord']))
                right_ok = False
                too_far = True
            if (not left_ok or not right_ok) and i != df.shape[0]-1:
                # print('==1')
                print(f'Primer pair #{i+1} has alternative binding site')
                continue
            else:
                pass
                # print(f'Primer pair #{i+1} has alternative binding site')
                
            if too_far == False:
                # print(has_multiple_binding_sites(row['pLeft_Sequences'], genome))
                # print(has_multiple_binding_sites(reverse_complement_sequence(row['pRight_Sequences']), genome))
                left_ok = not has_multiple_binding_sites(row['pLeft_Sequences'], genome)
                right_ok = not has_multiple_binding_sites(reverse_complement_sequence(row['pRight_Sequences']), genome)
                if (not left_ok or not right_ok) and i != df.shape[0]-1:
                    # print('==2')

                    print(f'Primer pair #{i+1} has alternative binding site')
                    continue  
                else:
                    pass
                    # print(f'Primer pair #{i+1} has alternative binding site')
            else:
                pass
            
            # print(too_far, left_ok, right_ok)

            #if pass all filtering
            if left_ok == True and right_ok == True:
                primer_pool.append(complement_sequence(row['pLeft_Sequences']))
                primer_pool.append(row['pRight_Sequences'])
                row_df = pd.DataFrame(row).T
                accepted_primers = pd.concat([accepted_primers, row_df],axis=0)
                print(f'-> Primer pair #{i+1} accepted')
                # print(low_b, row['pLeft_coord'])
                # print(abs(low_b - row['pLeft_coord']) > read_size/2)
                no_primer.append('-')
                # print('***')
                break
            else:
                if i == df.shape[0]-1:
                    print(f'!!!No suitable primer found: please manually inspect the sequence')
                    
                    # change = input('Continue with redesign and do not skip?(y/n):')
                    # if change == 'y' or '':
                    #     print('pass')

                    # else:
                    #     print('skip')
                    #     x = [0,row['pLeft_coord'],0,0,0,0,0,0,row['pRight_coord'],0,0,0,0,0,0,0]

                    #     accepted_primers.loc[len(accepted_primers)] = x
                    #     status = 'Skipped'
                    #     no_primer.append(status)
                    #     continue
                    # print('pass2')
                    
                    # skipping the amplicon by default
                    print('Current amplicon skipped')
                    
                    status = 'Redesigned'
                    x = [0,row['pLeft_coord'],0,0,0,0,0,0,row['pRight_coord'],0,0,0,0,0,0,0]
                    accepted_primers.loc[len(accepted_primers)] = x
                    status = 'Skipped'
                    no_primer.append(status)
                    continue
                    status = 'Redesigned'
                    
                    if abs(low_b - row['pLeft_coord']) > read_size/2:
                    # if abs(low_b - row['pLeft_coord']) > read_size/2 or ~abs(low_b - row['pLeft_coord'])>0:
                        print('!Problem with *left(forward)* primer', _is_recursive)
                        if _is_recursive=='right':
                            pass
                        else:
                            left_ok = False
                            print('!Problem with *left(forward)* primer')
                            print('How should I moved the left range? (e.g. -50 = moving start of covered range 50bp upstream)')
                            change = input('Where to move (+/-bps):')
                            if change == 'p' or change == 'pass':
                                x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                                accepted_primers.loc[len(accepted_primers)] = x
                                status = 'Skipped'
                                no_primer.append(status)
                                continue
                                
                            else:
                                seq_template = extract_sequence_from_fasta(low_b, high_b, padding=0, fasta_file=ref_genome)
                                print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} = {low_b-int(change), high_b}')
                                # primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, priority=priority, read_size = read_size)
                                primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome, high_b, low_b, read_size, priority, check_snp, freq_cutoff=50000, _is_recursive='left')

                    if abs(high_b - row['pRight_coord']) > read_size/2:
                    # if abs(high_b - row['pRight_coord']) > read_size/2 or ~(high_b - row['pRight_coord'])<0:
    
                        print('!Problem with *right)* primer', _is_recursive)

                        if _is_recursive=='left':
                            pass
                        else:
                            right_ok = False
                            print('!Problem with *right(backward)* primer')
                            print('How should I moved the right range? (e.g. -50 = moving start of covered range 50bp upstream)')
                            change = input('Where to move (+/-bps):')
                            if change == 'p' or change == 'pass':
                                x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                                accepted_primers.loc[len(accepted_primers)] = x
                                status = 'Skipped'
                                no_primer.append(status)
                                continue
                            else:
                                high_b = high_b+int(change)
                                seq_template = extract_sequence_from_fasta(low_b, high_b+int(change), padding=0, fasta_file=ref_genome)
                                print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} instead of {low_b, high_b-int(change)}')
                                # primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, priority=priority, read_size = read_size)
                                primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome, high_b, low_b, read_size, priority, check_snp, freq_cutoff=50000, _is_recursive='right')
                                continue
                else:
                    print(f'Primer pair #{i+1} has alternative binding site')
                    continue
        # print(primer_pool, accepted_primers)
        # print(accepted_primers)
        # print(df.iloc[0])
        # print(1)
    else:
        #print('Checking for homodimer')
        print(f'{df.shape[0]} primers designed')
        for i, row in df.iterrows():
            # print(row)
            left_ok = True
            right_ok = True
            too_far = False
            hetero = False
            # checking if the primer is too far away from the original range
            if check_snp:
                if check_for_snp(row['pLeft_Sequences'], priority, row['pLeft_length'], True, ref_genome) or check_for_snp(row['pRight_Sequences'], priority, row['pRight_length'], False, ref_genome):
                    print(f'Primer pair #{i+1} binding position has SNP')
                    continue
            if abs(low_b - row['pLeft_coord']) > read_size/2:
            # if abs(low_b - row['pLeft_coord']) > read_size/2 or ~(low_b - row['pLeft_coord'])>0:
                
                # print('low')
                # print(low_b, row['pLeft_coord'])
                # print(abs(low_b - row['pLeft_coord']))
                left_ok = False
                too_far = True
            if abs(high_b - row['pRight_coord']) > read_size/2:
            # if abs(high_b - row['pRight_coord']) > read_size/2 or ~(high_b - row['pRight_coord'])<0:
                # print('high')
                # print(high_b, row['pRight_coord'])
                # print(abs(high_b - row['pRight_coord']))
                right_ok = False
                too_far = True
            if (not left_ok or not right_ok) and i != df.shape[0]-1:
                # print('==4')
                
                print(f'Primer pair #{i+1} has alternative binding site')
                continue
            else:
                pass
                # print(f'Primer pair #{i} has alternative binding site')
#######################
            # if too_far == False:
            #     for x in primer_pool: # heterodimer check
            #         if check_heterodimer(x, complement_sequence(row['pLeft_Sequences'])) == False:
            #             left_ok = False
            #             hetero = True
            #         if check_heterodimer(x, row['pRight_Sequences']) == False:
            #             right_ok = False
            #             hetero = True

            #     if (not left_ok or not right_ok) and i != df.shape[0]-1:
            #         print(f'Primer pair #{i+1} has homodimer')
            #         continue  
            #     else:
            #         pass
            # else:
            #     pass
#######################

            
            for x in primer_pool: # heterodimer check
                if check_heterodimer(x, row['pLeft_Sequences'], global_args_dict) == False:
                    left_ok = False
                    hetero = False
                else:
                    hetero = True
                if check_heterodimer(x, row['pRight_Sequences'], global_args_dict) == False:
                    right_ok = False
                    hetero = False
                else:
                    hetero = True
##############################
            
            if too_far:
                print('Primer binding site too far from target region')
                continue
            elif hetero:
                print(f'Primer pair #{i+1} has heterodimer')
                continue
            else:  
                left_ok = not has_multiple_binding_sites(row['pLeft_Sequences'], genome)
                right_ok = not has_multiple_binding_sites(reverse_complement_sequence(row['pRight_Sequences']), genome)
                if (not left_ok or not right_ok) and i != df.shape[0]-1:
                    # print('==5')
                    print(f'Primer pair #{i+1} has alternative binding site')
                    continue
                else:
                    pass
            # print(too_far, hetero, left_ok, right_ok)
            if left_ok == True and right_ok == True:
                row_df = pd.DataFrame(row).T
                # primer_pool.extend(row[['pLeft_Sequences','pRight_Sequences']].values.tolist())
                primer_pool.append(row['pLeft_Sequences'])
                primer_pool.append(row['pRight_Sequences'])
                # print(row)
                accepted_primers = pd.concat([accepted_primers, row_df],axis=0)
                print(f'-> Primer pair #{i+1} accepted')
                # print(low_b, row['pLeft_coord'])
                # print(abs(low_b - row['pLeft_coord']) > read_size/2)
                no_primer.append('-')
                # print('***')
                break
            else:
                if i == df.shape[0]-1:
                    print(f'!!!No suitable primer found: please manually inspect the sequence')
                    
                    #skipping the amplicon, ask the user
                    # change = input('Continue with redesign and do not skip?(yes / No,skip this amplicon):')
                    # if change == 'y':
                    #     print('pass')
                    # else:
                    #     print('skip')
                    #     x = [0,row['pLeft_coord'],0,0,0,0,0,0,row['pRight_coord'],0,0,0,0,0,0,0]
                    #     accepted_primers.loc[len(accepted_primers)] = x
                        
                    #     status = 'Skipped'
                    #     no_primer.append(status)
                    #     continue
                    # print('pass2')
                    
                    # skipping the amplicon by default
                    print('Current amplicon skipped')
                    
                    status = 'Redesigned'
                    x = [0,row['pLeft_coord'],0,0,0,0,0,0,row['pRight_coord'],0,0,0,0,0,0,0]
                    accepted_primers.loc[len(accepted_primers)] = x
                    status = 'Skipped'
                    no_primer.append(status)
                    continue

                    if abs(low_b - row['pLeft_coord']) > read_size/2:
                    # if abs(low_b - row['pLeft_coord']) > read_size/2 or ~(low_b - row['pLeft_coord'])>=0:
                        
                        left_ok = False
                        print('!Problem with *left* primer', _is_recursive)

                        if _is_recursive=='right':
                            pass
                        else:
                            print('!Problem with *left(forward)* primer')
                            print('How should I moved the left range? (e.g. -50 = moving start of covered range 50bp upstream)')
                            change = input('Where to move (+/-bps):')
                            # if change == 'p' or change == 'pass':
                            #     status = 'Passed'
                            #     pass
                            # else:
                            if change == 'p' or '':
                                x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                                accepted_primers.loc[len(accepted_primers)] = x
                                status = 'Skipped'
                                no_primer.append(status)
                                continue
                            
                            
                            low_b = low_b+int(change)
                            seq_template = extract_sequence_from_fasta(low_b, high_b, padding=0, fasta_file=ref_genome)
                            print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} instead of {low_b-int(change), high_b}')
                            primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome, high_b, low_b, read_size, priority, check_snp, freq_cutoff=50000, _is_recursive='left')
                        
                    if abs(high_b - row['pRight_coord']) > read_size/2:
                    # if abs(high_b - row['pRight_coord']) > read_size/2 or ~(high_b - row['pRight_coord'])<0:
    
                        print('!Problem with *right)* primer', _is_recursive)
                        right_ok = False
                        if _is_recursive=='left':
                            pass
                        else:
                            if change == 'p' or '':
                                x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
                                accepted_primers.loc[len(accepted_primers)] = x
                                status = 'Skipped'
                                no_primer.append(status)
                                continue
                            else:
                                print('!Problem with *Right(backward)* primer')
                                print('How should I moved the right range? (e.g. -50 = moving start of covered range 50bp stream)')
                                change = input('Where to move (+/-bps):')
                                # if change == 'p' or change == 'pass':
                                #     status = 'Passed'
                                #     pass
                                # else:
                                high_b = high_b+int(change)
                                seq_template = extract_sequence_from_fasta(low_b, high_b, padding=0,fasta_file=ref_genome)
                                print(f'Redesigning primers for the new range ({change}bps): {low_b, high_b} instead of {low_b, high_b-int(change)}')
                                # primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome = ref_genome, high_b = high_b, low_b = low_b, priority=priority, read_size = read_size)
                                primer_pool, accepted_primers, no_primer = result_extraction(primer_pool, accepted_primers, seq_template, i+1, padding, ref_genome, high_b, low_b, read_size, priority, check_snp, freq_cutoff=50000, _is_recursive='left')
                                continue
                else:
                    print(f'Primer pair #{i+1} has alternative binding site')
                    continue
            #Alternative binding check
    return primer_pool, accepted_primers, no_primer

# %%
def extract_sequence_from_string(start_pos, end_pos, padding=150, sequence_string='', sequence_id='Chromosome'):
    """
    Extracts a subsequence from a given sequence string based on the start position and end position.
    """
    # Convert the sequence string to a Bio.Seq object
    sequence = Seq(sequence_string)
    
    # Extract the subsequence based on the start and end positions
    subsequence = sequence[start_pos - padding:end_pos + padding]
    
    return str(subsequence)  # Return the subsequence as a string
# %%
