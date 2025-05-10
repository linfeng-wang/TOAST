# TOAST (Tuberculosis Optimized Amplicon Sequencing Tool)
<!-- - TB and Other pathogen Amplicon Sequencing Tool
- TB ONT Amplicon Sequencing Tool -->

Here we introduced TOAST software tool aimed at addressing the challenges in amplicon primer design for TB sequencing leveraging the package primer3 for Tm, homopolymers, hairpins, homodimers considerations and in-house pipeline for avoiding heterodimer, alternative binding. This automated tool in takes user defined SNP priority for amplicon coverage and outputs designed amplicon primers with respective Tm and primer coordinates and sequence with capability of focusing on specific genes and taking into account of spoligotypes 

### How to install
Install the required conda environment
```conda env create -n TOAST -f environment.yml```
Clone repository
```git clone <repo html link>```
Install Python package (be in the root directory of the repository)
```pip install .```


It can also be installed through pip at https://pypi.org/project/toast-amplicon/

### Workflow
#### Before running the tool
*Decide on SNP priority by modifying the SNP priority file (variant.csv)
*Decide on amplicon size

#### Installing environment
- Install the required conda environment
    ```conda env create -n TOAST -f environment.yml```
- Install the tool (be in the directory with the .toml file)
    ```pip install .```

1. Estimate amplicon number needed for coverage (*amplicon_no* function)
   - Example: 
   ```
    toast amplicon_no -a 800 -op ./cache/Amplicon_design_output -g
   ```
2. Run amplicon design (*design* function)
    - Example: 
    ```    
    toast design -op ./cache/Amplicon_design_output -a 400 -sn 1 -sg rpoB,katG -nn 40 
    
    toast design -op ./cache/Amplicon_design_output -a 400 -sn 1 -sg rpoB,katG -nn 25

    toast design -op ./cache/output -a 1000 -nn 4 -ud ./cache/test_df.csv

    toast design -op ./cache/output -a 1000 -nn 26

    toast design -op ./cache/Amplicon_design_output -a 400 -sn 1 -sg rpsL -nn 0 -ud ./cache/test_df.csv
    ```
3. Check amplicon design using coverage plot (*plotting* function)
    - Example: 
    ```
    toast plotting -ap ./toast/Amplicon_design_output/Primer_design-accepted_primers-23-400.csv -rp ./toast/db/reference_design.csv -op ./cache/Amplicon_design_output -r 400
    ``` 


## Primer3 Configuration Parameters (default file: db/default_primer_design_setting.txt)

- **PRIMER_NUM_RETURN**: Number of primer pairs to return.
- **PRIMER_PICK_INTERNAL_OLIGO**: Flag to pick internal oligos (0 for no, 1 for yes).
- **PRIMER_INTERNAL_MAX_SELF_END**: Maximum self-complementarity score for internal oligos.
- **PRIMER_MIN_SIZE**: Minimum primer size in bases.
- **PRIMER_MAX_SIZE**: Maximum primer size in bases.
- **PRIMER_MIN_TM**: Minimum melting temperature (Tm) for primers in °C.
- **PRIMER_MAX_TM**: Maximum melting temperature (Tm) for primers in °C.
- **PRIMER_MIN_GC**: Minimum GC content in percent for primers.
- **PRIMER_MAX_GC**: Maximum GC content in percent for primers.
- **PRIMER_MAX_POLY_X**: Maximum length of mononucleotide repeats in primers.
- **PRIMER_INTERNAL_MAX_POLY_X**: Maximum length of mononucleotide repeats in internal oligos.
- **PRIMER_SALT_MONOVALENT**: Concentration of monovalent salts (e.g., Na+, K+) in mM.
- **PRIMER_DNA_CONC**: Concentration of DNA template in nM.
- **PRIMER_MAX_NS_ACCEPTED**: Maximum number of unknown bases (N's) accepted in primers.
- **PRIMER_MAX_SELF_ANY**: Maximum overall self-complementarity score for primers.
- **PRIMER_MAX_SELF_END**: Maximum 3' end self-complementarity score for primers.
- **PRIMER_PAIR_MAX_COMPL_ANY**: Maximum overall complementarity score between primer pairs.
- **PRIMER_PAIR_MAX_COMPL_END**: Maximum 3' end complementarity score between primer pairs.
- **PRIMER_PRODUCT_SIZE_RANGE**: Range of acceptable primer product sizes (e.g., "100-300").

## Example format of the user defined files can be found in ```user_defined_files/``` folder:

  - Configuration Parameters file: ```default_primer_design_setting.txt```
  - User input primer file: ```user_input_primer.csv```



## Ouput file format
```
<filetype>-<number of total amplicon designed>-<minimum amplicon size>-<maximum amplicon size>-<step size>-<number of amplicon for each size>
```

### 5 different files are produced
- Primer_design-accepted_primers: All detailed information about the designed amplicons
- Amplicon_importance: Number of SNP coverd by each amplicon
- Amplicon_mapped: bed file can be used to visualised the amplicon on genome using tools such as igv
- SNP_inclusion: shows SNP covered
- Gene_covereage: show percentage of each gene covered



**Specific mutation (Mutation Priority)file format:**
Essentially all you need would be the genome position (genome_pos). Other columns are needed but you could used imputed values like below if unknown.
The complete example Mutation priority csv can be found in *mutation_priority_example.csv*

| sample_id | **genome_pos** | gene   | change   | freq | type | sublin | drtype | drugs | weight |
|-----------|------------|--------|----------|------|------|--------|--------|-------|--------|
| sample_1  | 321168     | gene_1 | change_1 | 1    | -    | -      | -      | -     | 1      |
| sample_2  | 551767     | gene_2 | change_2 | 1    | -    | -      | -      | -     | 1      |
| sample_3  | 1017188    | gene_3 | change_3 | 1    | -    | -      | -      | -     | 1      |
| sample_4  | 1119158    | gene_4 | change_4 | 1    | -    | -      | -      | -     | 1      |
| sample_5  | 1119347    | gene_5 | change_5 | 1    | -    | -      | -      | -     | 1      |
| sample_6  | 1414872    | gene_6 | change_6 | 1    | -    | -      | -      | -     | 1      |

A script (*mutation_priority_gen.py*) can also be found to generate a file like the above:

example usage: 
```
python mutation_priority_gen.py --positions "322168,553767,1077188" --output <output_path.csv>
```







Test runs:
toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_jp -a 1000 -sn 1 -sg mmpR5 -nn 20
toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_jp -a 1000 -sn 1 -sg mmpR5 -nn 30
toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_jp -a 1000 -sn 1 -sg mmpR5 -nn 40


toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_seg -seg 600,800,50,2 -sg Rv0678
toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_seg -seg 300,800,50,3 -sg Rv0678,katG
toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_seg -seg 400,800,50,1 -sg Rv0678,katG




toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_Dario -a 600 -nn 19

toast design -s /mnt/storage10/lwang/Projects/TOAST_/toast/db/snp_priority_Dario.csv -op /mnt/storage10/lwang/Projects/TOAST_/cache/Amplicon_design_output_Dario -nn 19 -a 300 -p 100

toast design -s /mnt/storage10/lwang/Projects/TOAST_/toast/db/snp_priority_Dario6_1.csv -op /mnt/storage10/lwang/Projects/TOAST_/cache/Amplicon_design_output_Dario6_1 -nn 6 -a 300 -p 100
toast design -s /mnt/storage10/lwang/Projects/TOAST_/toast/db/snp_priority_Dario6_2.csv -op /mnt/storage10/lwang/Projects/TOAST_/cache/Amplicon_design_output_Dario6_2 -nn 6 -a 300 -p 100
toast design -s /mnt/storage10/lwang/Projects/TOAST_/toast/db/snp_priority_Dario6_3.csv -op /mnt/storage10/lwang/Projects/TOAST_/cache/Amplicon_design_output_Dario6_3 -nn 6 -a 300 -p 100


toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output_seg -seg 300,800,50,3 -sg Rv0678,katG

/mnt/storage10/lwang/Projects/TOAST

toast amplicon_no -s /mnt/storage10/lwang/Projects/TOAST/cache/lineage_specific_snp_priority.csv -a 650 -c 1 -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_no_output_lineage


toast design -s /mnt/storage10/lwang/Projects/TOAST/cache/lineage_specific_snp_priority.csv -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_lineage -seg 300,850,50,3 -p 100


toast design -s /mnt/storage10/lwang/Projects/TOAST/cache/lineage_specific_snp_priority_01-29.csv -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_lineage01-29 -seg 200,450,50,2 -p 100


toast design -sg KatG,inhA,FabG1,ahpC,rpoB,embA,embB,ubiA,gyrA,gyrB,pncA,rrl,rplc,thyA,folC,ribD,rpsl,rrs,gid,eis,ethA,Rv0678,pepQ,atpE,mmpl5,ddn,fbiA,fbiB,fbiC,fgd1,Rv2983,alr,Rv0440,Rv3874 -a 280 -p 200 -sn 30 -op /mnt/storage10/lwang/Projects/TOAST/cache/test1