**Z-score Pruning** 

GITHUB Code: https://github.com/neuromol/Z-prune
To account for the SE, summary statistics are converted into Z-scores which is then used in the pruning process (high absolute Z-score in a given correlated pair ) 

This code produces a file like; 
*SNP Pval abs_Z OR SE beta Z
rs2977608 0.2676 1.1086064370584916 0.9859 0.012809189020063547 -0.0142003494011414 -1.1086064370584916
rs112618790 0.8231 0.2235596032380278 1.0045 0.02008370567768239 0.004489905272852001 0.2235596032380278
rs3829740 0.3284 0.9773417305949796 1.0113 0.011497137157079342 0.011236631925987768 0.9773417305949796*

The input file needs to be in the following format; 
*SNP OR P
rs2977608 0.9859 0.2676
rs12562034 1.0028 0.8698
rs112618790 1.0045 0.8231*

Reference files needed; 
chromosomes.ld.h5  This file is just the plink output from doing the LD pruning using an R2 of 0.2 across the 1000G bam files (which can be downloaded from the LDSC ftp site ). h5 format makes it faster to dump the files into RAM per a chromosome. 

python_plink.sif (singularity image with all the python packages etc stored) - will post somewhere

Script INFO: 
> singularity exec ~/python_plink.sif python zprune.py

#########################################################
Pruning based off Z-scores
#########################################################

usage: zprune.py [-h] -i INPUT -d LD_FILE [-o OUTPUT] [-r2 RSQUARED]
[-t THREADS]
zprune.py: error: the following arguments are required: -i/--input, -d/--LD_file

NOTE: lowest R2 value can be 0.2


**TO DO** 
- [x] Upload working code
- [ ] Put reference files up online
- [ ] Redo plink files with lower R2

