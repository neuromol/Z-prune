#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z-score pruning 
"""
#########
from gPRS import process_input
from gPRS import ld_prune
import pandas as pd 
import numpy as np 
import os
from scipy.stats import norm
from datetime import datetime
import statsmodels.api as sm
import statsmodels.formula.api as smf
import argparse
import sys
from warnings import simplefilter
simplefilter(action='ignore', category=FutureWarning)
import concurrent.futures as cf
import itertools
#########
sys.path.append(os.getcwd())
startTime = datetime.now()
def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Z-score pruning' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='Input file name ',required='True')
    parser.add_argument('-d', '--LD_file', help='column name ',required='True')
    parser.add_argument('-o', '--output', help='output directory, default current',default='"z-score.output"')
    parser.add_argument('-r2', '--rsquared', help='r2 cutoff default= 0.3',default=0.3)
    parser.add_argument('-t', '--threads', help='number of threads to use',default=1)
    results = parser.parse_args(args)
    return (results.input , results.LD_file, results.output, results.rsquared, results.threads )

def main(file_input , LD_directory, output, R_squared, threads):
    threads = int(threads)
    real_data =(process_input.open_test(file_input) ) 
    SNP_list = pd.DataFrame()
    chromosomes = [1, 2 , 3 , 4, 5 ,6 ,7 , 8 ,9 , 10, 11, 12, 13, 14 ,15 ,16 ,17 ,18, 19, 20, 21 , 22]
    #chromosomes = [15]
    with cf.ProcessPoolExecutor(max_workers=int(threads)) as executor:
        for dataframe in executor.map(ld_prune.process_chr, itertools.repeat(real_data), itertools.repeat(R_squared) , itertools.repeat(LD_directory) , chromosomes):
            ### append to dataframe from concurrent futures
            SNP_list= SNP_list.append(dataframe)
    print (" " )
    print ("There are " + str(len(SNP_list)) + " rows in the prune SNP list")
    print (" " )
    SNP_list.to_csv(output, sep="\t" , index=None )
    EndTime = datetime.now()
    Total_time = EndTime - startTime 
    print ("Total analysis time = " + str(Total_time))
    
if __name__ == '__main__':
    print( """
    #########################################################
    Pruning based off Z-scores 
    #########################################################
    """)
    file_input , LD_file, output, R_squared, threads = check_arg(sys.argv[1:])
    main (file_input , LD_file, output, R_squared, threads)
        
        
        
        
        
        
        
        
        
