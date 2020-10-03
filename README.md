# imputation_sims

Converts fsc26 output to plink files (PED and MAP)

Works for simulated DNA seqs in the form of snps (numeric strings).
Considers individuals to be diploid, so
1_1 and 1_2 individuals from .arp files are considered to belong to the same individual. The plink output will be in the form of:
1_1
1_3 
etc...

Input parameters:
- file.arp (Fsc26 arlequin output file)
- output name
- number of haploid individuals (ie. 1000 haploid individuals = 500 diploid individuals)
