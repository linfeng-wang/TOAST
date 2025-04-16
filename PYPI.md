# TOAST (Tuberculosis Optimized Amplicon Sequencing Tool)
<!-- - TB and Other pathogen Amplicon Sequencing Tool
- TB ONT Amplicon Sequencing Tool -->

We present TOAST, a software tool designed to streamline and optimize amplicon primer design for Mycobacterium tuberculosis sequencing. TOAST integrates the robust primer design capabilities of Primer3—accounting for Tm, homopolymers, hairpins, and homodimers—with an in-house pipeline that rigorously filters for heterodimer formation and unintended alternative binding. What sets TOAST apart is its automation and intelligence: it leverages a curated database of over 50 M. tuberculosis genomes to inform amplicon placement, ensuring robust primer performance across strain diversity. Users can prioritize SNPs for coverage, focus on specific resistance genes, and tailor designs to spoligotype backgrounds. The tool outputs primer sequences along with detailed thermodynamic profiles and genomic coordinates, making it an end-to-end solution for targeted TB panel design.

--- 

## Key Functionalities
### 1. **Amplicon Design (`toast design`)**
**Main Inputs:**
- SNP priority lists (`-s`): Default is globally collated clinical TB SNPs
- Reference genome files (`-ref`): Default is MTB-h37rv genome
- Spoligotype sequencing range files (`-sp_f`)
- Optional user-defined primers (`-ud`) and custom genomic features (`-gff`)

**Configurable Settings:**
- Amplicon size (`-a`)
- Padding around target regions (`-p`): Default is amplicon size divided by 6
- Number of specific amplicons (`-sn`) and targeted gene names (`-sg`)
- Number of non-specific amplicons (`-nn`)
- Graphical output option (`-g`) to visualize amplicon coverage
- **All SNPs** in the reference genome (`-all_snp`) uses all frequent SNPs from the default database to design degenerate primers. **For other species**, a custom SNP file in the same format can be provided to override the default.
**Outputs:**  
- Amplicon sequences and primer details organized into a user-specified output directory (`-op`)
- Optional graphical representations of designed amplicons

**Example Usage:**
```bash
toast design -op ./output -a 400 -sn 2 -sg rpoB,katG -nn 20
```
This command designs amplicons of 400 base pairs, including two specifically targeting the `rpoB` and `katG` genes, and 20 additional amplicons for prioritized SNP coverage.

---

### 2. **Amplicon Number Estimation (`toast amplicon_no`)**

**Purpose:**  
Estimate the number of amplicons required to achieve desired SNP coverage in TB genomic studies.

**Main Inputs:**
- SNP priority files (`-s`)
- Reference genome (`-ref`)

**Settings:**
- Desired amplicon size (`-a`)
- Target coverage depth
- Optional graphical output for coverage estimates (`-g`)

**Outputs:**  
Estimates and coverage graphics saved in the specified output directory (`-op`).

---

### 3. **Visualization and Plotting (`toast plotting`)**

**Purpose:**  
Visualize and analyze the coverage and distribution of designed amplicons.

**Main Inputs:**
- SNP priority files (`-s`)
- Genomic feature files (GFF, `-gff`)
- Primer sequences and reference designs

**Settings:**
- Read size specifications

**Outputs:**  
Visualization graphics, including coverage plots, available in the specified output directory (`-op`).

---

## Quick Start

To view available command-line options and their defaults, use:
```bash
toast design -h
toast amplicon_no -h
toast plotting -h
```
---

### Workflow
#### Before running the tool
*Decide on SNP priority by modifying the SNP priority file (mutation_priority_example.csv - can be found in github)
*Decide on amplicon size

#### Dependencies

  - "pandas", 
  - "numpy",
  - "plotly",
  - "rich_argparse",
  - "tabulate",
  - "primer3-py >= 2.0.1"
  - "python=3.11"

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

Segmented amplicon design is used to generate amplicons of varying sizes so they can be easily distinguished on an agarose gel. This provides a quick visual check to confirm successful amplification before sequencing, allowing users to validate that each target produced a distinct band. It helps avoid wasting sequencing resources on failed reactions and serves as a practical sanity check in the experimental workflow.

If issues arise during amplicon design—such as no primers being generated—it is most likely due to insufficient padding. A small padding size can restrict the available sequence context needed for effective primer placement. To resolve this, try increasing the padding value to provide more room for the design algorithm to work with.

---

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

---

**Specific mutation (Mutation Priority)file format:**
Essentially all you need would be the genome position (genome_pos). Other columns are needed but you could used imputed values like below if unknown.
The complete example Mutation priority csv can be found in Github：*mutation_priority_example.csv*

| sample_id | **genome_pos** | gene   | change   | freq | type | sublin | drtype | drugs | weight |
|-----------|------------|--------|----------|------|------|--------|--------|-------|--------|
| sample_1  | 321168     | gene_1 | change_1 | 1    | -    | -      | -      | -     | 1      |
| sample_2  | 551767     | gene_2 | change_2 | 1    | -    | -      | -      | -     | 1      |
| sample_3  | 1017188    | gene_3 | change_3 | 1    | -    | -      | -      | -     | 1      |
| sample_4  | 1119158    | gene_4 | change_4 | 1    | -    | -      | -      | -     | 1      |
| sample_5  | 1119347    | gene_5 | change_5 | 1    | -    | -      | -      | -     | 1      |
| sample_6  | 1414872    | gene_6 | change_6 | 1    | -    | -      | -      | -     | 1      |

You can manually eddit this for though a script (Github：*mutation_priority_gen.py*) can also be found to generate a file like the above:

example usage: 
```
python mutation_priority_gen.py --positions "322168,553767,1077188" --output <output_path.csv>
```

---
**More complete information with example files can be found in the Github repo**
---


REFERENCE:
Wang, L., Naphatcha Thawong, Thorpe, J., Higgins, M., Tan, M., Waritta Sawaengdee, Surakameth Mahasirimongkol, Perdigao, J., Campino, S., Clark, T.G. and Phelan, J.E. (2025). A novel tool for designing targeted gene amplicons and an optimised set of primers for high-throughput sequencing in tuberculosis genomic studies. bioRxiv (Cold Spring Harbor Laboratory). doi:https://doi.org/10.1101/2025.01.13.632698.
