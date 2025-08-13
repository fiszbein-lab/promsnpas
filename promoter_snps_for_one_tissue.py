import argparse
import gzip
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import pickle
import time
import os 
import subprocess

# from tqdm import tqdm
from matplotlib.lines import Line2D
from scipy.stats import ks_2samp

import sys
sys.path.insert(-1, '/projectnb/encore/sofyaga/promoters_position_dependency/' )
import gtf_parser as gtfp
import rmats_parser as rmatsp
from vcf_parser import collect_snps_from_vcf, snp_lists_to_df



def one_tissue_iteration(tissue, genes_with_exons, samples_to_tissues, 
                        path_to_01_file, path_to_12_file, columns_01, columns_12):
    
    # load all the rmats samples into gene entries
    print(f'loading all rmats data for {tissue} into gene entries')
    tt = time.time()
    for sample_id in os.listdir(f'/projectnb/encore/gtex/data/{tissue}/'):
        if '.' in sample_id:
            continue
        if "gen3-client" in sample_id:
            continue
            
        path_to_sample_dir = f'/projectnb/encore/gtex/data/{tissue}/{sample_id}/rmats/'
        rmatsp.parse_rmats_sample(genes_with_exons, path_to_sample_dir, sample_id)
    print(time.time() - tt)
    tt = time.time()
    
    # calculate statistics with snp data
    print(f'calculating correlation stats of snps with rmats')
    for gene in genes_with_exons:
        
        if len(genes_with_exons[gene].snp_data) == 0:
            continue
            
        genes_with_exons[gene].dump_rmats_data_buffer(samples_to_tissues)    
        genes_with_exons[gene].merge_snp_w_rmats()
        genes_with_exons[gene].get_snp_correlation(tissue= tissue)
    print(time.time() - tt)
        
    open_01_file = open(path_to_01_file, 'a')
    open_12_file = open(path_to_12_file, 'a')
        
        
    # save only the very significant into a file 
    print('saving significant results')
    tt = time.time()
    for gene in genes_with_exons:
        
        if len(genes_with_exons[gene].snp_data) == 0:
            continue
            
        temp = genes_with_exons[gene].snp_correlation[tissue]
        temp = temp[np.abs(temp['diff_0_1']) > DIFFTHRES]
        temp = temp[temp['pval_0_1'] < PVALTHRES]
        temp = temp[temp['len_1'] > SAMPLETHRESH]
        temp = temp[temp['len_0'] > SAMPLETHRESH]
        temp['gene'] = gene
        temp['tissue'] = tissue
        temp[columns_01].to_csv(path_or_buf=open_01_file, sep='\t', columns=None, header=False, index=False) 
        
        temp = genes_with_exons[gene].snp_correlation[tissue]
        temp = temp[np.abs(temp['diff_1_2']) > DIFFTHRES]
        temp = temp[temp['pval_1_2'] < PVALTHRES]
        temp = temp[temp['len_1'] > SAMPLETHRESH]
        temp = temp[temp['len_2'] > SAMPLETHRESH]
        temp['gene'] = gene
        temp['tissue'] = tissue
        temp[columns_12].to_csv(path_or_buf=open_12_file, sep='\t', columns=None, header=False, index=False)
     
    print(time.time() - tt)
    
    open_01_file.close()
    open_12_file.close()
    
    
    
#     clean the rmats data from gene entries
    print('cleaning gene rmats data')
    tt = time.time()
    for gene in genes_with_exons:
        
        genes_with_exons[gene].clean_rmats_data()
        genes_with_exons[gene].clean_snp_correlation(tissue)
    
    print(time.time() - tt)
        
#     return genes_with_exons
        
        


# Define the parser
parser = argparse.ArgumentParser(description='Short sample app')
parser.add_argument('--tissue', action="store", dest='TISSUE', default='NONE')
args = parser.parse_args()
TISSUE = args.TISSUE
print('will be processing: ', TISSUE)

# get the tissue list
tissues = os.listdir('/projectnb/encore/gtex/data/')

if TISSUE not in tissues:
    print('tissue is wrong')
    exit()

lengths = list()
for tis in tissues:
    lengths.append(len(os.listdir(f'/projectnb/encore/gtex/data/{tis}')))
    
    
# And the list of samples analyzed    
tissues_to_samples = dict()
samples_to_tissues = dict()

for tissue in tissues:
    tissues_to_samples[tissue] = list()
    runs = [run for run in os.listdir(f'/projectnb/encore/gtex/data/{tissue}/') if ('.json' not in run)]
    for run_id in runs:
        tissues_to_samples[tissue].append(run_id)
        assert (run_id in samples_to_tissues) == False
        samples_to_tissues[run_id] = tissue
        
        
        
        
        
# # promoter_db = pd.read_csv('promoters_from_hit.csv', sep = '\t')
# promoter_db = pd.read_csv('promoters_from_hit.csv', sep = '\t')
# promoter_db['promoter_id'] = promoter_db.apply(lambda row: 
#                              row['chromosome']+':'+str(row['genome_start'])+'-'+str(row['genome_end']), 
#                              axis = 1)
        
        
        
        
# load the genes with exons data 
tt = time.time()
# with open('genes_with_promoter_snps.pickle', 'rb') as handle:
with open('genes_with_encodepromoter_snps.pickle', 'rb') as handle:
    genes_with_exons = pickle.load(handle)
print('Loaded the genes with snp data, ', time.time() - tt)

      
    
DIFFTHRES = 0.1
PVALTHRES = 0.01
SAMPLETHRESH = 5



for gene in genes_with_exons:
    sample_columns = genes_with_exons[gene].sample_columns
    num_snps = (genes_with_exons[gene].snp_data[sample_columns] > 0).sum(axis = 1)
    enough_snps = num_snps > SAMPLETHRESH*2
    genes_with_exons[gene].snp_data = genes_with_exons[gene].snp_data.loc[enough_snps].reset_index()
    if len(genes_with_exons[gene].snp_data) > 0:
        genes_with_exons[gene].snp_ids = genes_with_exons[gene].snp_data['ID'].values



columns_01  = ['exon_start', 'exon_end', 'upes', 'upee', 'does', 'doee', 'snpID', 
              'diff_0_1', 'pval_0_1', 'len_0', 'len_1', 'inc_0', 'inc_1', 'inc_2', 'gene', 'tissue']
path_to_01_file = f'encodePLS_SNP_data_01_{TISSUE}.tsv'

open_file = open(path_to_01_file, 'w')
open_file.write('\t'.join(columns_01))
open_file.write('\n')
open_file.close()


path_to_12_file = f'encodePLS_SNP_data_12_{TISSUE}.tsv'
columns_12 = ['exon_start', 'exon_end', 'upes', 'upee', 'does', 'doee', 'snpID', 
              'diff_1_2', 'pval_1_2', 'len_1', 'len_2', 'inc_0', 'inc_1', 'inc_2', 'gene', 'tissue']

open_file = open(path_to_12_file, 'w')
open_file.write('\t'.join(columns_12))
open_file.write('\n')
open_file.close()


one_tissue_iteration(TISSUE, genes_with_exons, samples_to_tissues, 
                        path_to_01_file, path_to_12_file, columns_01, columns_12)

# one_tissue_iteration(TISSUE, genes_with_exons, samples_to_tissues, 
#                         path_to_0_12_file, path_to_01_2_file, columns_0, columns_01)

        
        
        
        
        
        