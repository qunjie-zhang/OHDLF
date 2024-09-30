#!/bin/bash

# First, use RNA-seq.sh to perform transcript assembly and splicing for the SRR files.
# for example 
bash RNA-seq.sh SRR1170742 

# Then, name the obtained protein files according to the species.
# for example
sed 's/^>/\>bar|/' SRR1170742_trinity_cdhit99.fasta.transdecoder.pep > bar.pep

# After that, place all the named files in the same directory and run OrthoFinder.
# for example
orthofinder -f directory

# finaly, You need to place OHDLF script within the Orthofinder output directory named 'Results_XXX' and run
# for example
python OHDLF.py -l 0.1 -d 10 -s 95 -p 1

# output
# If you choose to use the concatenation method, the script will output a result file named 'final_OrthologsAlign_GDL.phy',which will be used for tree construction with RAXML later on. 
# If you opt for the Coalescence method, the script will output a result file named 'all.trees', which will be used for tree construction with ASTRAL subsequently.