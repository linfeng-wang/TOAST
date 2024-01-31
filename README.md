# TOAST (Tuberculosis Optimized Amplicon Sequencing Tool)
<!-- - TB and Other pathogen Amplicon Sequencing Tool
- TB ONT Amplicon Sequencing Tool -->

Here we introduced TOAST software tool aimed at addressing the challenges in amplicon primer design for TB sequencing leveraging the package primer3 for Tm, homopolymers, hairpins, homodimers considerations and in-house pipeline for avoiding heterodimer, alternative binding. This automated too in takes user defined SNP priority for amplicon coverage and outputs designed amplicon primers with respective Tm and primer coordinates and sequence with capability of focusing on specific genes and taking into account of spoligotypes 

### Workflow
#### Before running the tool
*Decide on SNP priority by modifying the SNP priority file (variant.csv)
*Decide on amplicon size

#### Installing environment
- Install the required conda environment
    ```conda env create -n TOAST -f environment.yml```
    
#### Running the tool
Get to the code directory
```cd code```

1. Estimate amplicon number needed for coverage (*amplicon_no* function)
   - Example: 
   ```
    python main.py amplicon_no -a 400 -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -g

    toast amplicon_no -a 400 -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -g
   ```
2. Run amplicon design (*design* function)
    - Example: 
    ```
    python main.py design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -a 400 -sn 1 -sg rpoB,katG -nn 40 
    
    toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -a 400 -sn 1 -sg rpoB,katG -nn 40 
    
    toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -a 400 -sn 1 -sg rpoB,katG -nn 25

    toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/output -a 400 -nn 40 
    toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/output -a 1000 -nn 26

    toast design -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -a 400 -sn 1 -sg rpsL -nn 0 -ud /mnt/storage10/lwang/Projects/TOAST/cache/test_df.csv
    ```
3. Check amplicon design using coverage plot (*plotting* function)
    - Example: 
    ```
    python main.py plotting -ap /mnt/storage10/lwang/Projects/TOAST/toast/Amplicon_design_output/Primer_design-accepted_primers-23-400.csv -rp /mnt/storage10/lwang/Projects/TOAST/toast/db/reference_design.csv -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -r 400
    
    toast plotting -ap /mnt/storage10/lwang/Projects/TOAST/toast/Amplicon_design_output/Primer_design-accepted_primers-23-400.csv -rp /mnt/storage10/lwang/Projects/TOAST/toast/db/reference_design.csv -op /mnt/storage10/lwang/Projects/TOAST/cache/Amplicon_design_output -r 400
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

## Example format of the user defined files can be found in '''user_defined_files/''' folder:

  - Configuration Parameters file: '''default_primer_design_setting.txt'''
  - User input primer file: '''user_input_primer.csv'''