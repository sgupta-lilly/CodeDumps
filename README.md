# CodeDumps
siteseq scripts

** This is a modification of the Camerson 2007 et. al script **
1. Cameron_DumpSites.py

Run command:
 python Cameron_DumpSites.py down_test.bam GRCh38.primary_assembly.genome.fa outfile



Input:
1. aligned bam file
2. reference genome
3. name of output file
Outputs:
1. output file which is in fasta format with all potential offtarget and DBS sites
2. all peak coordinates above 5 reads detected and reported in a json file


** All sites from siteseq **

The purpose of db_secondary.py is to read all the genomics coordinates reported in the experiment by siteseq and macs2. Identify overlaps and create a bed file so we can look and the overall coverage across all samples for all identified sites 
python db_secondary.py path_to_siteseq path_to_macs2_broadpeak

Input:
1. the path to the siteseq directory
2. the path to the macs2 directory

Outputs:
1. record of genomics coordinates from siteseq and macs2. Also, reported if the locations are detected by both methods
2. bed file for all sites. 



The purpose of this file is to read the samplesheet and generate a exhausted list of comparisons between each treatment to all controls 

python macs2_combinations.py samplesheet.csv

Input:
1. samplesheet.csv

Output:
1. macs2combo file which contains trtment_sample,control_sample,trtment_sample_vs_control_sample
