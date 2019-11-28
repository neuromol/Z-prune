#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
from datetime import datetime
import numpy as np 
from scipy.stats import norm
from warnings import simplefilter
simplefilter(action='ignore')
from sklearn.utils import shuffle

def open_test(file_input):
    """
    Open the real GWAS results with SNP, OR, Pval,A1 , A2 
    """
    startTime = datetime.now()
    test = pd.read_csv(file_input , sep="\t", header=0 )
    #ref_file = pd.read_csv(ref_file,header=0)
    test.columns=["SNP" , "OR" , "Pval"]
    test = test[test["Pval"] != 1 ] 
    test = test.set_index("SNP")
    test["SE_LC"] = np.absolute(np.log(test.OR)/norm.ppf(test.Pval/2))
    test["beta_LC"] = np.log(test.OR)
    test["Z_LC"] = test.beta_LC / test.SE_LC
    test["absZ_LC"]= np.absolute(test.Z_LC)
    Na_sum = test.absZ_LC.isna().sum() 
    print ("Dropping " + str(Na_sum) + " SNPs due to NaN Z-scores")
    print (" " )
    test = test.drop(["SE_LC", "beta_LC" , "Z_LC"], axis=1)
    test = test.fillna(0)
    #print (test[test.isna().any(axis=1)])
    EndTime = datetime.now()
    Total_time = EndTime - startTime  
    print ("Processed input file in : " + str(Total_time))
    print ("   ")
    return test
    