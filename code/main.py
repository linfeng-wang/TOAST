from functools import partial
from random import choices, randint, randrange, random, sample
from typing import List, Optional, Callable, Tuple
import numpy as np
# from geneticalgorithm import geneticalgorithm as ga
from rich_argparse import ArgumentDefaultsRichHelpFormatter
import pandas as pd
from collections import Counter
from tqdm import tqdm
import time
from Bio.SeqUtils import MeltingTemp
from Bio import SeqIO
from plotly import graph_objects as go
import json
from importlib import reload
import primer_selection
reload(primer_selection)
import Amplicon_no
reload(Amplicon_no)
import working_algo_gene_2in1 as w
reload(w)
import argparse
from functools import reduce
import os
import matplotlib.pyplot as plt
import pandas as pd
import plotting1 as p
from icecream import ic
from tabulate import tabulate

def main(args):
    print('>>>Designing Amplicons')
    #test run: python main.py design -s variants.csv -ref MTB-h37rv_asm19595v2-eg18.fa -op ../test -a 400 -p 200 -sn 0 -sg rpoB,katG,embB,pncA,rpsL,rrs,ethA,fabG1,gyrA,gid,inhA,ethR,rpoC,ahpC,gyrB,folC,tlyA,alr,embA,thyA,eis -nn 30 -g -sp -sp_f spacers.bed 
    #test run: python main.py design -op ../test -a 400 -sn 1 -sg rpoB,katG -nn 1
    #test run: python main.py design -op ../test -a 400 -sn  -sg rpoB,katG -nn 35 -sc
    # python main.py design -op ../test -a 400 -sn 1 -sg rpoB,katG -nn 1 -sc
    # python main.py design -op ../test -a 400 -sn 1 -sg rpoB,katG -nn 1 -sc
    # python main.py design -op ../test -a 1000 -sn 2 -sg rpoB,katG -nn 15 -sc
    # python main.py design -op ../test -a 400 -nn 0 -sc
    # Rationale:
    # This program is designed to perform primer design for PCR amplification,
    # taking into consideration various parameters such as amplicon size, SNP priority,
    # reference genome, and more. It also allows users to specify particular amplicons
    # based on gene names and offers graphical output options. Users can also choose to
    # include spoligotyping sequencing information. The program outputs results to a 
    # specified folder.

    # Displaying the User Settings for Verification:
    print("========== Design: User Settings ==========")
    print(f"Amplicon Size: {args.amplicon_size}")
    print(f"SNP Priority: {args.snp_priority}")
    print(f"Reference Genome: {args.reference_genome}")
    if args.padding_size == None:
        print(f"Padding_size: {args.amplicon_size/4}")
    else:
        print(f"Padding_size: {args.padding_size}")
    print(f"Specific Amplicon Number: {args.specific_amplicon_no}")
    if args.specific_amplicon_no > 0:
        print(f"Specific Amplicon Gene: {args.specific_amplicon_gene}")
    print(f"Non-specific Amplicon Number: {args.non_specific_amplicon_no}")
    print(f"Graphic Option: {args.graphic_option}")
    print(f"Spoligo Sequencing: {args.spoligo_coverage}")
    if args.spoligo_coverage:
        print(f"Spoligo Sequencing File: {args.spoligo_sequencing_file}")
    print(f"Output Folder Path: {args.output_folder_path}")
    
    print("===========================================")

    gene_names = [
        "rpoB",
        "katG",
        "embB",
        "pncA",
        "rpsL",
        "rrs",
        "ethA",
        "fabG1",
        "gyrA",
        "gid",
        "inhA",
        "ethR",
        "rpoC",
        "ahpC",
        "gyrB",
        "folC",
        "tlyA",
        "alr",
        "embA",
        "thyA",
        "eis"
    ]

    read_size = args.amplicon_size
    # full data - priority file with weights&frequency for each snp
    full_data= pd.read_csv(args.snp_priority)
    full_data = full_data[~full_data['drugs'].isna()]
    full_gene = full_data[~full_data['type'].isin(['synonymous_variant','non_coding_transcript_exon_variant'])]
    full_data = full_data.sort_values(by=['genome_pos'])
    full_data = full_data.reset_index(drop=True)
    full_data['weight'] = full_data['freq']
    full_data.loc[full_data['gene'].isin(gene_names), 'weight'] += 0.5

    # paddding size
    if args.padding_size == None:    
        padding = int(read_size/4)
    else:
        padding = args.padding_size

    # Reference Genome
    ref_genome = args.reference_genome
    
    #specific_genes
    if args.specific_amplicon_gene:
        specific_gene = args.specific_amplicon_gene.split(',')
        specific_gene = [item.strip() for item in specific_gene]
        for x in specific_gene:
            if x not in full_data['gene'].unique():
                print(f'Error: {x} is not a valid gene name')
                return 0
    else:
        specific_gene = []

    # non_specific_gene = []
    # specific_gene_amplicon
    specific_gene_amplicon = args.specific_amplicon_no
    # non_specific_amplicon
    non_specific_amplicon = args.non_specific_amplicon_no
    # this is way we separate the specific and non-specific amplicons
    # specific_gene_data = full_data[full_data['gene'].isin(specific_gene)]
    # non_specific_gene_data = full_data[~full_data['gene'].isin(specific_gene)]
    # this way we still have non specific amplicons incorporate the specific genes
    specific_gene_data = full_data.copy()
    # Create the condition for rows where 'gene' value is not in specific_gene
    condition = ~specific_gene_data['gene'].isin(specific_gene)
    # Update the 'weight' column in the new DataFrame where the condition is True
    specific_gene_data.loc[condition, 'weight'] = 0
    non_specific_gene_data = full_data.copy()
    ref_size = w.genome_size(ref_genome)
    specific_gene_data_count = 0
    #main output folder path
    output_path = args.output_folder_path

    # Calculating number of amplicon needed
    # target_coverage = 1
    # gene_coverage = Amplicon_no.place_amplicon_search(full_data, target_coverage, read_size, genome_size(ref_genome))

    covered_positions, covered_ranges = [], []
    primer_pool, no_primer_ = [], []
    accepted_primers = pd.DataFrame(columns = ['pLeft_ID', 'pLeft_coord', 'pLeft_length', 'pLeft_Tm', 'pLeft_GC',
        'pLeft_Sequences', 'pLeft_EndStability', 'pRight_ID', 'pRight_coord',
        'pRight_length', 'pRight_Tm', 'pRight_GC', 'pRight_Sequences',
        'pRight_EndStability', 'Penalty', 'Product_size'])

    if len(specific_gene)>0:
        print('=====Specific amplicon=====')
        
        # def place_amplicon(full_data, read_number, read_size, primer_pool, accepted_primers, no_primer_, ref_genome, graphic_output=False, padding=150, output_path = '.'):
    
        covered_positions_sp, covered_ranges_sp, specific_gene_data, primer_pool, accepted_primers, no_primer_ = w.place_amplicon(specific_gene_data, specific_gene_amplicon, read_size, primer_pool, accepted_primers, no_primer_, ref_genome, args.graphic_option, padding=padding, output_path =output_path)
        specific_gene_data_count = accepted_primers.shape[0]                                                  
        print('=====Non-specific amplicon=====')
        for (x, y) in covered_ranges_sp:
            condition = (non_specific_gene_data['genome_pos'] >= x) & (non_specific_gene_data['genome_pos'] <= y)
            non_specific_gene_data.loc[condition, 'weight'] = 0

        # non_specific_gene_data.update(specific_gene_data) # add the specific gene data to the non-specific gene data
        covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = w.place_amplicon(non_specific_gene_data, non_specific_amplicon, read_size, primer_pool, accepted_primers, no_primer_, ref_genome, args.graphic_option, padding=padding, output_path =output_path)
        covered_positions = {**covered_positions_sp, **covered_positions_nosp}
        
        covered_ranges = covered_ranges_sp + covered_ranges_nosp
        non_specific_gene_data_count = accepted_primers.shape[0] - specific_gene_data_count
        
    else:
        if args.non_specific_amplicon_no > 0:
            print('=====Non-specific amplicon=====')
            covered_positions_nosp, covered_ranges_nosp, full_data_cp, primer_pool, accepted_primers, no_primer_ = w.place_amplicon(non_specific_gene_data, non_specific_amplicon, read_size, primer_pool, accepted_primers, no_primer_,ref_genome, args.graphic_option, padding=padding, output_path =output_path)
            covered_positions = covered_positions_nosp
            covered_ranges = covered_ranges_nosp
            non_specific_gene_data_count = accepted_primers.shape[0]
        else:
            pass

    # whether or not you wanna include spoligotyping sequencing
    spoligotype = args.spoligo_coverage
    covered_ranges_spol = []
    if spoligotype:            

        print('=====Spoligotype amplicon=====')

        spacers = pd.read_csv(args.spoligo_sequencing_file, sep='\t', header=None)
        spacers = np.array(spacers)
        spacers = spacers[:, 1:3]
        spacers = spacers.tolist()
        flattened_data = [item for sublist in spacers for item in sublist]
        spacer_max = max(flattened_data)
        spacer_min = min(flattened_data)
        spol_list = np.arange(spacer_min-400,spacer_max+400,1)
        weight = [0.01]*len(spol_list) 
        spol_data = pd.DataFrame({'genome_pos':spol_list,'weight':weight})
        # Create a list of boolean masks, one for each range
        masks = [(spol_data['genome_pos'] >= start) & (spol_data['genome_pos'] <= end) for start, end in spacers]
        # Use reduce and the | operator to combine the masks into a single mask
        combined_mask = reduce(lambda x, y: x | y, masks)
        # Use .loc and the combined mask to update the weight column
        spol_data.loc[combined_mask, 'weight'] = 1
        covered_positions_spol, covered_ranges_spol, full_data_cp, primer_pool, accepted_primers, no_primer_ = Amplicon_no.place_amplicon_spol(spol_data, 1, read_size, ref_genome, primer_pool, accepted_primers, no_primer_, padding=padding, graphic_output=False, check_snp=False)
        covered_ranges = covered_ranges + covered_ranges_spol
        print(covered_ranges)
        print(covered_ranges_spol)
        # covered_ranges_spol = Amplicon_no.place_amplicon_spol(spol_data, 1, read_size, graphic_output=False, ref_size = w.genome_size(ref_genome))
        # covered_ranges.extend(covered_ranges_spol)

    read_number = specific_gene_data_count + non_specific_gene_data_count + len(covered_ranges_spol)
    # output
    primer_label = [f'Gene_specific:{args.specific_amplicon_gene}']*specific_gene_data_count + ['Non_specific']*non_specific_gene_data_count + ['Spoligotype']*len(covered_ranges_spol)
    # print(primer_label)
    # print(len(primer_label))
    # print(accepted_primers)
    # print(accepted_primers.shape)
    # print(covered_ranges)
    # print(no_primer_)
    print(covered_ranges)
    print(accepted_primers)
    print(no_primer_)
    accepted_primers['Amplicon_type'] = primer_label
    # accepted_primers['Redesign'] = no_primer_
    accepted_primers['Designed_ranges'] = covered_ranges
    accepted_primers.reset_index(inplace = True)
    
    #accepted_primers - change primer to iupac
    threshold = 0.001
    all_snps = pd.read_csv(args.all_snps, sep = '\t', header = None)
    all_snps.drop_duplicates(inplace=True)
    for i, row in accepted_primers.iterrows():
        # print(row['pLeft_Sequences'])
        # print(row['pRight_Sequences'])
        #left primer
        primer_seq = ''
        for x,y in zip(range(row['pLeft_coord'], row['pLeft_coord']+len(row['pLeft_Sequences'])), row['pLeft_Sequences']):
            if (all_snps[all_snps[0] == x][[1,2]].shape[0] > 0) and (all_snps[all_snps[0] == x][[3]].values.astype('float').item()>= threshold):
                alleles = ''.join(all_snps[all_snps[0] == x][[1,2]].values[0])
                primer_seq = primer_seq+w.nucleotide_to_iupac(alleles)
            else:
                primer_seq = primer_seq+y
        if row['pLeft_Sequences'] != primer_seq:
            # print('SNP in the Left primer')
            accepted_primers.loc[i, 'pLeft_Sequences'] = primer_seq
            
        #right primer
        primer_seq = ''
        for x,y in zip(range(row['pRight_coord'], row['pRight_coord']+len(row['pRight_Sequences'])), primer_selection.reverse_complement_sequence(row['pRight_Sequences'])):
            if all_snps[all_snps[0] == x][[1,2]].shape[0] > 0 and all_snps[all_snps[0] == x][[3]].values.astype('float').item()>= threshold:
                alleles = ''.join(all_snps[all_snps[0] == x][[1,2]].values[0])
                primer_seq = primer_seq+w.nucleotide_to_iupac(alleles)
            else:
                primer_seq = primer_seq+y
        # print(primer_seq)
        if row['pRight_Sequences'] != primer_seq:
            # print('SNP in the Right primer')
            accepted_primers.loc[i, 'pRight_Sequences'] = primer_seq
    
    # print(accepted_primers[['pLeft_Sequences','pRight_Sequences','pLeft_coord']].head(2))
    if spoligotype:
        sp = '-sp'
    else:
        sp = ''

    # Apply modifications
    accepted_primers['pLeft_ID'] = accepted_primers.apply(lambda x: w.modify_primer_name(x['pLeft_ID'], x['Amplicon_type'], 'L'), axis=1)
    accepted_primers['pRight_ID'] = accepted_primers.apply(lambda x: w.modify_primer_name(x['pRight_ID'], x['Amplicon_type'], 'R'), axis=1)

    op = f'{output_path}/Amplicon_design_output'
    os.makedirs(op, exist_ok=True) #output path
    accepted_primers.to_csv(f'{op}/Primer_design-accepted_primers-{read_number}-{read_size}{sp}.csv',index=False)

    primer_pos = accepted_primers[['pLeft_coord','pRight_coord']].values
    columns = ['pLeft_ID', 'pRight_ID', 'pLeft_coord', 'pRight_coord', 'SNP_inclusion']

    # Create an empty DataFrame with the specified column headings
    primer_inclusion =pd.DataFrame(columns=columns)
    for i, row in accepted_primers.iterrows():
        data = full_data[(full_data['genome_pos']>= row['pLeft_coord']) & (full_data['genome_pos']<= row['pRight_coord'])]    
        info = row[['pLeft_ID', 'pRight_ID', 'pLeft_coord', 'pRight_coord']]
        SNP = data['gene'].str.cat(data['change'], sep='-').unique()
        info['SNP_inclusion'] = ','.join(SNP)
        primer_inclusion.loc[len(primer_inclusion)] = info.tolist()

    primer_inclusion.to_csv(f'{op}/SNP_inclusion-{read_number}-{read_size}{sp}.csv',index=False)


    if specific_gene_amplicon>0 or non_specific_amplicon>0:
        # amp_snp = pd.DataFrame(columns=full_data.columns)
        
        dtypes = {}
    # Iterate over each column in the DataFrame
        for column in full_data.columns:
            # Get the data type of the column
            dtype = str(full_data[column].dtype)
            # Add the column and its data type to the dictionary
            dtypes[column] = dtype
        amp_snp = pd.DataFrame(columns=dtypes.keys()).astype(dtypes)

        _len = accepted_primers[accepted_primers['Amplicon_type']!= 'Spoligotype'].shape[0]
        if _len == 0:
            pass
        else:
            for (x,y) in covered_ranges[:_len]:
                condition = (full_data['genome_pos'] >= x) & (full_data['genome_pos'] <= y)
                filtered_data = full_data[condition]
                # Drop columns that are completely empty or all NA from filtered_data
                filtered_data = filtered_data.dropna(axis=1, how='all')
                amp_snp = pd.concat([amp_snp, filtered_data])
            
            amp_snp = amp_snp.drop_duplicates()
            
            full_data_gene = full_data['gene'].value_counts()
            amp_snp_aligned = amp_snp['gene'].value_counts().reindex(full_data_gene.index).fillna(0)

            # Calculate the ratio and format as string
            result = amp_snp_aligned / full_data_gene
            formatted_result = result.apply(lambda x: f"{x:.0%}")

            # Combine the two series
            combined_series = amp_snp_aligned.astype(str) + "/" + full_data_gene.astype(str) + " (" + formatted_result + ")"

            # print(combined_series)

            # gene_coverage = round(amp_snp['gene'].value_counts()/full_data['gene'].value_counts()*100, 1)

            gene_coverage_df = combined_series.reset_index().fillna(0)
            gene_coverage_df.columns = ['Gene', 'SNP_coverage']
            
            print(tabulate(gene_coverage_df, headers='keys', tablefmt='grid'))

    print('Primer design output files:')
    print(f'{op}/SNP_inclusion-{read_number}-{read_size}.csv')
    print(f'{op}/Primer_design-accepted_primers-{read_number}-{read_size}.csv')
    
    return 0

def main_amplicon_no(args):
    print("========== Amplicon Number: User Setting ==========")
    print(f"Amplicon Size: {args.amplicon_size}")
    print(f"SNP Priority: {args.snp_priority}")
    print(f"Reference Genome: {args.reference_genome}")
    print(f"Target Coverage: {args.target_coverage}")
    print(f"Graphic Option: {args.graphic_option}")
    print(f"Output Folder Path: {args.output_folder_path}")
    print("===================================================")
    # test run: python main.py amplicon_no -s variants.csv -ref MTB-h37rv_asm19595v2-eg18.fa -op ../test -a 400 -c 1 -sp -sp_f spacers.bed
    # test run: python main.py amplicon_no -a 400 -op ../test -g
    # test run: python main.py amplicon_no -a 1000 -op ../test -g
    full_data = pd.read_csv(args.snp_priority)
    full_data = full_data[~full_data['drugs'].isna()]
    full_data = full_data.sort_values(by=['genome_pos'])
    full_data = full_data.reset_index(drop=True)
    full_data['weight'] = full_data['freq']
    ref_genome = args.reference_genome
    target_coverage = args.target_coverage
    read_size = args.amplicon_size
    output_path = args.output_folder_path
    graphics = args.graphic_option
    gene_coverage = Amplicon_no.place_amplicon_search(full_data, target_coverage, read_size, w.genome_size(ref_genome),output_path, graphics)
    return 0
    
def main_plotting(args):
    # test run: python main.py plotting -s variants.csv -gff MTB-h37rv_asm19595v2-eg18.gff -ap ../test/Primer_design-accepted_primers-30-400.csv -rp ../test/Primer_design-accepted_primers-30-400.csv -op ../test
    # test run: python main.py plotting -ap ../test/Amplicon_design_output/Primer_design-accepted_primers-42-400.csv -rp ../db/reference_design.csv -op ../test -r 400
    # test run: python main.py plotting -ap ../test/Amplicon_design_output/Primer_design-accepted_primers-30-1000.csv -rp ../db/reference_design.csv -op ../test -r 1000
    print("========== Plotting: User Settings ==========")
    print(f"SNP Priority: {args.snp_priority}")
    print(f"Accepted Primers: {args.accepted_primers}")
    print(f"GFF Feature file: {args.gff_features}")
    print(f"Reference Design: {args.reference_design}")
    print(f"Read Size: {args.read_size}")
    print(f"Output Folder Path: {args.output_folder_path}")
    print("============================================")
    priority = args.snp_priority
    accepted_primers = args.accepted_primers
    gff = args.gff_features
    reference_design = args.reference_design
    output_dir = args.output_folder_path
    read_size = args.read_size
    p.plotting(priority, read_size, accepted_primers, gff, reference_design, output_dir)
    return 0
    
# %%
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='Amplicon designer for TB', 
                                    description='Amplicon design, list of specific genes that can be priotised: rpoB,katG,embB,pncA,rpsL,rrs,ethA,fabG1,gyrA,gid,inhA,ethR,rpoC,ahpC,gyrB,folC,tlyA,alr,embA,thyA,eis',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(dest="command", help="Task to perform")

    ###### Design pipeline
    parser_sub = subparsers.add_parser('design', help='Run whole design pipeline', formatter_class=ArgumentDefaultsRichHelpFormatter)
    input=parser_sub.add_argument_group("Input options")
    # parser.add_argument('-c','--country_file', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='variants.csv', default=None)
    # in
    input.add_argument('-s','--snp_priority', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='../db/variants.csv')
    input.add_argument('-ref','--reference_genome', type = str, help = 'reference fasta file (default: MTB-h37rv genome)', default='../db/MTB-h37rv_asm19595v2-eg18.fa')
    input.add_argument('-sp_f','--spoligo_sequencing_file', type = str, help = 'Custom spoligotype range files (default: TB spligotype space ranges)', default = '../db/spacers.bed')
    input.add_argument('-as','--all_snps', type = str, help = 'All SNPs in the reference genome', default = '../db/all_snps.csv')
    
    #design
    setting=parser_sub.add_argument_group("Design options")
    setting.add_argument('-a','--amplicon_size', type = int, help = 'Amplicon size', default=400)
    setting.add_argument('-p','--padding_size', type = int, help = 'Size of padding on each side of the target sequence during primer design', default=None)
    setting.add_argument('-sn','--specific_amplicon_no', type = int, help = 'number of amplicon dedicated to amplifying specific genes', default=0 )
    setting.add_argument('-sg','--specific_amplicon_gene', type = str, help = 'give a list of gene names separated by Lineage ', default='')
    setting.add_argument('-nn','--non-specific_amplicon_no', type = int, help = 'number of amplicon dedicated to amplifying all SNPs in all genes according the importants list', default=20)
    setting.add_argument('-g','--graphic_option', action='store_true', help = 'output graphic on amplicon coverage to visualise the running of the algorithm', default = False)
    setting.add_argument('-sc','--spoligo_coverage', action='store_true', help = 'Whether to amplify Spoligotype', default = False)

    # out
    output=parser_sub.add_argument_group("Output options")
    output.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output_folder_path (accepted_primers, SNP_inclusion, gene_covered)', required=True)
    parser_sub.set_defaults(func=main)
    
    ###### Amplicon number estimates
    parser_sub = subparsers.add_parser('amplicon_no', help='Amplicon number estimates', formatter_class=ArgumentDefaultsRichHelpFormatter)
    input=parser_sub.add_argument_group("input options")
    input.add_argument('-s','--snp_priority', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='../db/variants.csv')
    input.add_argument('-sc_f','--spoligo_sequencing_file', type = str, help = 'Custom spoligotype range files (default: TB spligotype space ranges)', default = '../db/spacers.bed')
    input.add_argument('-ref','--reference_genome', type = str, help = 'reference fasta file (default: MTB-h37rv genome)', default='../db/MTB-h37rv_asm19595v2-eg18.fa')
    
    #design
    setting=parser_sub.add_argument_group("Amplicon options")
    setting.add_argument('-a','--amplicon_size', type = int, help = 'Amplicon size', default=400)
    setting.add_argument('-c','--target_coverage', type = int, help = 'target coverage of SNPs default: full coverage(1)', default=1)
    # setting.add_argument('-sc','--spoligotype_sequencing', action='store_true', help = 'Whether to do a separate run on chekcing the number of amplicon needed to cover the spligotypes', default= False)
    # setting.add_argument('-ap','--amplicon_sequencing', action='store_true', help = 'Whether to calculate amplicon number for SNP coverage', default = True)
    setting.add_argument('-g','--graphic_option', action='store_true', help = 'output graphic on amplicon coverage to visualise the running of the algorithm', default = False)
    
    # out
    output=parser_sub.add_argument_group("Output options")
    output.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output_folder_path (Amplicon number estimate graphics)', required=True)
    parser_sub.set_defaults(func=main)
    
    parser_sub.set_defaults(func=main_amplicon_no)

    # output=parser_sub.add_argument_group("Output options")
    # output.add_argument('-op','--output_folder_path', default = './', type = str, help = 'output folder path for covered ranges')
    # output.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output folder path for covered ranges')

    ###### Visualisation of designed amplicons
    parser_sub = subparsers.add_parser('plotting', help='Visualised the designed amplicons', formatter_class=ArgumentDefaultsRichHelpFormatter)
    #input
    input=parser_sub.add_argument_group("input options")
    input.add_argument('-s','--snp_priority', type = str, help = 'SNP priority CSV files (default: collated global 50k clinical TB samples)', default='../db/variants.csv')
    input.add_argument('-gff','--gff_features', type = str, help = 'genomic feature file .gff for the corresponding genome', default='../db/MTB-h37rv_asm19595v2-eg18.gff')
    input.add_argument('-ap','--accepted_primers', type = str, help = 'primer design output file from desgin function', required=True)
    input.add_argument('-rp','--reference_design', type = str, help = '(reference) design that can be plotted against the designed amplicons for comparision', default=None)
    input.add_argument('-r','--read_size', type = str, help = 'size of the designed amplicons', default='unknown-')
    
    #output
    output=parser_sub.add_argument_group("output options")
    output.add_argument('-op','--output_folder_path', default = '', type = str, help = 'output_folder_path (accepted_primers, SNP_inclusion, gene_covered)', required=True)
    parser_sub.set_defaults(func=main_plotting)

    args = parser.parse_args()
    if args.command == 'design':
        main(args)
    if args.command == 'amplicon_no':
        main_amplicon_no(args)
    if args.command == 'plotting':
        main_plotting(args)