# TOAST
Name options:
- TB Optimized Amplicon Sequencing Tool
<!-- - TB and Other pathogen Amplicon Sequencing Tool
- TB ONT Amplicon Sequencing Tool -->

Here we introduced Amplicon4TB software tool aimed at addressing the challenges in amplicon primer design for TB sequencing leveraging the package primer3 for Tm, homopolymers, hairpins, homodimers considerations and in-house pipeline for avoiding heterodimer, alternative binding. This automated too in takes user defined SNP priority for amplicon coverage and outputs designed amplicon primers with respective Tm and primer coordinates and sequence with capability of focusing on specific genes and taking into account of spoligotypes 

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
    python main.py amplicon_no -a 400 -op ../test1 -g
   ```
2. Run amplicon design (*design* function)
    - Example: 
    ```
    python main.py design -op ../test1 -a 400 -sn 1 -sg rpoB,katG -nn 40 -sc
    ```
3. Check amplicon design using coverage plot (*plotting* function)
    - Example: 
    ```
    python main.py plotting -ap ../test1/Amplicon_design_output/Primer_design-accepted_primers-42-400.csv -rp ../db/reference_design.csv -op ../test1 -r 400
    ``` 
