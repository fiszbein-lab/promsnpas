import os 
import numpy as np 
import pandas as pd
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import scipy.stats as sps

import promoter_df_collapse as prom_collapse


def intersection(ls, le, rs, re):
    assert ls <= le, f'error in intersection function: {ls} > {le}'
    assert rs <= re, f'error in intersection function: {rs} > {re}'
    inter =  min(le, re) - max(ls, rs)
    if inter < 0:
        return False
    return True

def covers(s, e, to_cover):
    assert s <= e
    
    if (s <= to_cover) and (to_cover <= e):
        return True
    
    assert (to_cover < s) or (e < to_cover)
    
    return False

# this is a class to organize together the annotation and the hit index output data 

# Each gene contains exon entries, that have unique start and end coordinates. 
# the exon entries are super small and meaningless, now they are a class just in case, 
# if we will need to speed this up, we might switch to arrays or smth. 
# 
# Now we have two things to combine: 
#    1) the genome annotation, that has exon numbers, and out of which we actually need all the first exons
#    2) the hit index output, from which we get the metaexon coordinates, if it is first of hybrid or internal, 
#       the number of upstream and downstream reads, hit index and things like that. 
# For the 2 we will create a MetaExon category. 
# When processing annotation we will add to the Exon entries if it is first in gencode annotation.
# Importantly, each gene will have an array of exons, and this array will be sorted by the exon starts. 


class Exon:
    def __init__(self, start, end, ori):
        self.s = int(start)
        self.e = int(end)
        self.ori = ori 
        self.positions = list()
    
    def is_first(self):
        if 1 in self.positions:
            return True
        return False
    
    def start(self):
        assert self.s <= self.e
        if self.ori == '+':
            return self.s
        else:
            return self.e
            
    
class MetaExon(Exon):
    def __init__(self, start, end, ori):
        super().__init__(start, end, ori)
        self.sample_data = list()
        
        
    def add_sample_data(self, sample, up_reads=0, down_reads=0,  
                 HITindex=None, pval_internal=None, ID_position=None, reads_in_sample=None):
        
        self.sample_data.append([sample, up_reads, down_reads, 
                                 HITindex, pval_internal, ID_position, reads_in_sample])

    def create_dataframe(self):
        self.sample_df = pd.DataFrame(self.sample_data, 
                                     columns = ['sample', 'U', 'D', 
                                                'HIT', 'p_internal', 'ID_pos', 'reads_in_sample'])
        
        

class GTFGene: 
    def __init__(self, gene_id, chromosome, start, end, ori): 
        self.gene_id = gene_id
        self.chromosome = chromosome
        self.gene_start = int(start)
        self.gene_end = int(end)
        self.ori = ori
        self.exons = dict()
        
    def add_exon_with_number(self, start, end, ori, num):
        start, end = int(start), int(end)
        if (start in self.exons) == False:
            self.exons[start] = dict()
            
        if (end in self.exons[start]) == False:
            self.exons[start][end] = Exon(start, end, ori)
                
        assert type(num) == int, 'Please provide exon number as an integer'
        assert num > 0, 'Please provide exon number as > 1'
        self.exons[start][end].positions.append(num)
    
    
class HITGene(GTFGene):
    def __init__(self, gene_id, chromosome, start, end, ori): 
        super().__init__(gene_id, chromosome, start, end, ori)
        
        self.metaexons = dict()
        
    
    def add_metaexon_data(self, start, end, ori):
        
        if (start in self.metaexons) == False:
            self.metaexons[start] = dict()
            
        if (end in self.metaexons[start]) == False:
            self.metaexons[start][end] = MetaExon(start, end, ori)
            
    def get_promoters(self):
        promoters = list()
        
        for s in self.exons:
            for e in self.exons[s]:
                exon = self.exons[s][e]
                if exon.is_first():
                    promoters.append(exon.start())
        
        promoters = sorted(promoters)
        if self.ori == '-':
            promoters = promoters[::-1]
        self.promoters = promoters
        
        return promoters

    def promoter_number(self):
        return len(self.promoters)
        
    
    def get_promoter_data(self, samples_in_tissue = 3, be_strict = True):
        promoters = self.promoters
        promoter_data = list()
        
        for start in self.metaexons:
            for end in self.metaexons[start]:
                
                
                promoters_in_exon = list()
                for ip, p in enumerate(promoters):
                    if covers(start, end, to_cover=p):
                        promoters_in_exon.append(p)
                
                self.metaexons[start][end].create_dataframe()
                df = self.metaexons[start][end].sample_df
                df['tissue'] = df.apply(lambda row: row['sample'].split('GTEX')[0], axis = 1)

                if be_strict:
                    df['isf'] = df['ID_pos'] == 'first'
                    if max(df.groupby(['tissue'])['isf'].sum()) < samples_in_tissue:
                        continue
#                     if ('first' in df['ID_pos'].values) != True:
#                         continue
                    if sum(df['D'] - df['U'] > 0) <= samples_in_tissue: # this would also mean the hit index to be more than 0 
                        continue
                      

                df['exon_starts'] = ','.join(map(str, promoters_in_exon))
                df['exon_position_mean'] = np.mean(promoters_in_exon) if len(promoters_in_exon) else np.nan
                
                df['all_promoters'] = ','.join(map(str, promoters_in_exon))
                
                if self.ori == '-':
                    df['promoter'] = end
                else:
                    df['promoter'] = start
                df['meta_exon_s'] = start
                df['meta_exon_e'] = end
                promoter_data.append(df)
        
        if len(promoter_data) == 0:
            self.promoter_data_df = pd.DataFrame([], 
                                                 columns = ['sample', 'U', 'D', 'HIT', 'p_internal', 
                                                            'ID_pos', 'reads_in_sample', 'tissue',
                                                            'exon_starts',
                                                            'exon_position_mean', 'promoter',
                                                            'all_promoters', 'meta_exon_s', 
                                                            'meta_exon_e'])
        else:
            self.promoter_data_df =  pd.concat(promoter_data)  

        self.promoter_data_df = prom_collapse.collapse_promoters(self.promoter_data_df, self.ori, distance = 500)
        
        return self.promoter_data_df
    
    
    def get_by_tissue_promoter_data(self, tissues_to_samples, fill_negatives = False, read_normalization_recovery_coef = 17733575):
        
        promoter_by_tissue = list()
        if len(self.promoter_data_df) == 0:
            self.promoter_by_tissue = pd.DataFrame([], columns = ['tissue', 'promoter', 'mean_reads', 'std_reads', 'num_samples_expressed', 'num_samples_first'])
            return self.promoter_by_tissue
        
        assert 'tissue' in self.promoter_data_df.columns
        
        difference = self.promoter_data_df['D'] - self.promoter_data_df['U']
        if fill_negatives:
            difference[difference < 0] = 0
            
#         self.promoter_data_df['D-U'] = difference
        self.promoter_data_df['D-U'] = difference/(self.promoter_data_df['reads_in_sample']/read_normalization_recovery_coef)

        
        for tissue, tissue_group in self.promoter_data_df.groupby('tissue'):
            for promoter, tissue_promoter_group in tissue_group.groupby('promoter'):
                
                n_expressed = len(tissue_promoter_group)
#                 n_samples_first =  sum(tissue_promoter_group['ID_pos'].isin(['first', 'firstInternal', 'FirstInternal']))
                first_filter = tissue_promoter_group['ID_pos'].isin(['first', 'firstInternal', 'FirstInternal']).values
                print('here is the first filter: ', first_filter)
                n_samples_first = sum(first_filter)
                freads = tissue_promoter_group['D-U'].values[first_filter]
                print('here is the reads_first:\n', freads)
                print( tissue_promoter_group['D'].values[first_filter])
                print( tissue_promoter_group['U'].values[first_filter])
                print()
                
                n = len(tissues_to_samples[tissue]) - len(tissue_promoter_group)
                reads = list(tissue_promoter_group['D-U'].values) + [0] * n
                promoter_by_tissue.append([tissue, promoter, np.mean(reads), np.var(reads),
                                           np.mean(freads), np.var(freads),
                                           n_expressed, n_samples_first])
        
        self.promoter_by_tissue = pd.DataFrame(promoter_by_tissue, 
                                               columns = ['tissue', 'promoter', 
                                                          'mean_reads', 'std_reads', 
                                                          'mean_freads', 'std_freads',
                                                         'num_samples_expressed', 'num_samples_first'])
               
        return self.promoter_by_tissue        
   
        
def process_gtf_for_first_exons(gtf_file_path, gene_types = ['protein_coding'], num_genes = 3000, 
                               geneclass = HITGene):
    genes_with_exons = dict()
    curr_gene = None

    with open(gtf_file_path, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            splitted = line.split('\t')
            
            if splitted[2] == 'gene':
                if len(genes_with_exons) == num_genes:
                    return genes_with_exons
                
                
                
                chromosome, start, end, ori = splitted[0], splitted[3], splitted[4], splitted[6]
                annot = splitted[8].split(';')
                for elem in annot:
                    if 'gene_type' in elem:
                        gene_type =  elem.split('"')[1]
                    if 'gene_id' in elem:
                        curr_gene = elem.split('"')[1]
                if gene_type not in gene_types:
                    continue 
                    
                genes_with_exons[curr_gene] = geneclass(curr_gene, chromosome, start, end, ori)
                        
            if gene_type not in gene_types:
                continue
                        
            if splitted[2] == 'exon':
                num = -1
                chromosome, start, end, ori = splitted[0], splitted[3], splitted[4], splitted[6]
                annot = splitted[8].split(';')

                for elem in annot:
                    if 'gene_id' in elem:
                        assert curr_gene == elem.split('"')[1]
                    if 'exon_number' in elem:
                        num = elem.split(' ')[2]
                        if num.startswith('"'):
                            num = num[1:-1]
                        num = int(num)
                genes_with_exons[curr_gene].add_exon_with_number(start, end, ori, num)


    return genes_with_exons
