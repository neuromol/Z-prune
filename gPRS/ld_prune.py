#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np 
from scipy.stats import norm

def process_chr(df_input, R_squared , LD_file_path, chromosome):
    """ 
    This processes each chromosome to get the LD for each SNP and gets the best SNP
    at a given locus (500KB). set R threahold (original has 0.1) for each chromosome
    and cal the Z-score
    """
    #print ("processing chromosome " + str(chromosome))
    R_squared  = float(R_squared)
    key = "CHR" + str(chromosome)        
    temp = pd.read_hdf("chromsomes.ld.h5", key=key)
    temp.columns = ["CHR_A", "BP_A" , "SNP1", "CHR_B" ,"BP_B", "SNP2" ,"R" ] 
    temp.drop(["CHR_A", "BP_A", "CHR_B" ,"BP_B",], axis=1)
    print ("Getting variants with an R-squared greater than " + str(R_squared) + " for chromosome " + str(chromosome) ) 
    temp = temp[temp["R"] >= R_squared] 
    #### Join in the GWAS results
    merged = df_input.merge(temp , left_index=True , right_on="SNP1")
    merged.rename(columns={'absZ_LC':'Z_SNP1'}, inplace=True)
    merged2 = merged.merge(df_input, left_on="SNP2" , right_index=True)
    merged2.rename(columns={'absZ_LC':'Z_SNP2'}, inplace=True)
    remove_list = merged2[merged2["Z_SNP1"] < merged2["Z_SNP2"]]
    clean = merged2[~merged2.SNP1.isin(remove_list["SNP1"])]
    clean = clean[["SNP1" , "Pval_x" , "Z_SNP1", "OR_x"]].copy()
    clean.columns = ["SNP" , "Pval" , "abs_Z" , "OR"]
    clean = clean.drop_duplicates()
    ### recal the raw Z-scores 
    clean["SE"] = np.absolute(np.log(clean.OR)/norm.ppf(clean.Pval/2))
    clean["beta"] = np.log(clean.OR)
    clean["Z"] = clean.beta / clean.SE
    return clean