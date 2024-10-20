# OHDLF
a pipeline designed to filter and address orthologous gene heterogeneity, duplication, and loss
## dependencies
- Pyhton 3
- Biopython
- The script execution also requires other software such as Mafft and IQ-TREE. If needed, users can download the OHDLF.yaml file to directly configure the environment. 
## Quick start
```
Concatenation:
python OHDLF.py -l 0.05 -d 6 -s 97 -p 1
Coalescence:
python OHDLF.py -l 0.05 -d 6 -s 97 -p 2
```
## Usage
Use it as a terminal command
```
Program:    OHDLF (orthologous gene filter Workflow)
Version:    1.0

Useage: python OHDLF.py  <command> [options]

Commands:
  -l / --loss Allowable the max missing rate of gene. This option is required.
  -d / --duplication Allowable the max duplication number of gene. This option is required.
  -s / --similarity Allowable the similarity threshold of gene. If you do not set this parameter, the program will use '95' by default
  -p / --process_type process_type: 1 for Concatenation, 2 for Coalescence
```

## Input Setting
**Input** :You need to place this script within the Orthofinder output directory named 'Results_XXX'. The script depends on two directories: 'Orthogroup_Sequences' and 'Orthogroups'.

## Reading the Output
**Output** :If you choose to use the concatenation method, the script will output a result file named 'final_OrthologsAlign_GDL.phy', which will be used for tree construction with RAXML later on. If you opt for the Coalescence method, the script will output a result file named 'all.trees', which will be used for tree construction with ASTRAL subsequently.

Output intermediate files:

- **all_GDL_Orthologue_Sequences**: Generate all orthologous sequences that meet the criteria for gene duplication and loss, with the sequence files in FASTA format.
- **GDL_Orthologue_Sequences**: Generate low-copy orthologous gene sequences with a similarity of over 95% after BLASTP computation, with the sequence files in FASTA format.
- **GDL_Orthologue_Sequences_mafft**: Generate low-copy orthologous gene sequences that meet the similarity criteria and have been aligned using MAFFT, with the sequence files in FASTA format.
- **GDL_Orthologue_Sequences_mafft_fill**: Generate a FASTA file that has been filled and merged after MAFFT alignment, where the filling method is as follows: if multiple sequences have a dash (-) at the same position, that position is filled with a dash; if different amino acids appear at the same position, it is filled with an 'X'.
