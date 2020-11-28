CS124 Haplotype Phasing Algorithm
Michael Cheng (604928535) and Anweshan Das (940967021)

A haplotype phasing prediction algorithm given a genotype dataset of individuals (missing data or not) and returning a predicted haplotype phase pair for each genotype.

Input: genotype data - Text file with each row a SNP event and each column a different individual. Each column is separated by a ' ' and each row is separated by a '\n' newline.

Output: haplotype phases - Text file in the same format as the genotype data, with twice as many columns because each genotype is composed of a pair of haplotypes.

This script is run with python3

Prerequisites:
Must have fancyimpute, pandas, and numpy installed. Run this in the command line if you do not have the package installed yet:

pip3 install numpy
pip3 install pandas 
pip3 install fancyimpute

Usage/Running the Script:
Go in your command line to the directory where the script is stored. 
Run this line:

python3 haplotype_phasing.py test_data_masked.txt

The output file will be created in the same directory.