import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
import time

import gzip

import sys
sys.path.insert(-1, '/projectnb/encore/sofyaga/promoters_position_dependency/' )
import gtf_parser as gtfp


# return difference, pval, lengths
def get_influence_stats(data1, data2):
    l1, l2 = len(data1), len(data2)
    if l1*l2 == 0:
        return [np.nan, np.nan, l1, l2]
    difference = np.mean(data1) - np.mean(data2)
    return [difference, ks_2samp(data1, data2).pvalue, l1, l2]


    

class RmatsGene(gtfp.GTFGene):
    def __init__(self, gene_id, chromosome, start, end, ori): 
        super().__init__(gene_id, chromosome, start, end, ori)
        
        self.rmats_data = pd.DataFrame(columns = ['sample_id',
                                                'exonStart_0base', 'exonEnd', 
               'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
               'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IncLevel1', 'tissue'])

        self.rmats_row_buffer = list()
    
    def add_rmats_row_to_buffer(self, rmats_row):
        self.rmats_row_buffer.append(rmats_row)
        
    def dump_rmats_data_buffer(self, samples_to_tissues):
        
        new_rmats_data = pd.DataFrame(self.rmats_row_buffer, 
                                     columns = ['sample_id', 'exonStart_0base', 'exonEnd', 
                   'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
                   'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 
                   'IncLevel1'])
        
        if len(new_rmats_data) > 0:
            new_rmats_data['tissue'] = new_rmats_data.apply(lambda row: samples_to_tissues[row['sample_id']], 
                                                            axis = 1)
        
        self.rmats_data = pd.concat([self.rmats_data, new_rmats_data], ignore_index = True)
        self.rmats_row_buffer = list()
        
        integer_cols = ['exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 
                        'downstreamES', 'downstreamEE', 'IJC_SAMPLE_1', 'SJC_SAMPLE_1']
        float_cols = ['IncLevel1']
        str_cols = ['sample_id']
        self.rmats_data[integer_cols] = self.rmats_data[integer_cols].astype(int)
        self.rmats_data[float_cols] = self.rmats_data[float_cols].astype(float)
        self.rmats_data[str_cols] = self.rmats_data[str_cols].astype(str)

        
    def clean_rmats_data(self):
        del self.rmats_data
        del self.rmats_row_buffer
        self.rmats_data = pd.DataFrame(columns = ['sample_id',
                                                'exonStart_0base', 'exonEnd', 
               'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
               'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IncLevel1', 'tissue'])

        self.rmats_row_buffer = list()

        
class SNPRmatsGene(RmatsGene):
    def __init__(self, gene_id, chromosome, start, end, ori): 
        super().__init__(gene_id, chromosome, start, end, ori)
        self.snp_data = pd.DataFrame()
        self.snp_correlation = dict()
        self.sample_columns = list()
        
    def init_snp_data(self, snp_data):

        self.sample_columns = list()
        
        for col in snp_data.columns:
            if col.startswith('GTEX'):
                self.sample_columns.append(col)

        self.snp_data = snp_data.reset_index(drop = True)
    

    def merge_snp_w_rmats(self):
        if len(self.rmats_data) == 0:
            return 

        sample_columns = self.sample_columns
        values_transposed = self.snp_data[sample_columns].values.T
        self.snp_ids = self.snp_data['ID'].values
        
        snp_df_transposed = pd.DataFrame(values_transposed, 
                                           columns = self.snp_ids)
        snp_df_transposed['short_sample_name'] = self.sample_columns
        
        self.rmats_data['short_sample_name'] = self.rmats_data.apply(lambda row:
                                               '-'.join(row['sample_id'].split('-')[:2]), 
                                                        axis = 1)

        self.rmats_data = self.rmats_data.merge(snp_df_transposed, how = 'left', 
                           on = 'short_sample_name')


    def get_snp_correlation(self, tissue, specific_exons = None):

        snp_change = list()
        
        tissue_rmats = self.rmats_data.loc[self.rmats_data['tissue'] == tissue]

        for (exon_start, exon_end, upes, upee, does, doee), group in tissue_rmats.groupby(['exonStart_0base', 'exonEnd', 'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE']):
            if specific_exons != None:
                if (exon_start, exon_end, upes, upee, does, doee) not in specific_exons:
                    continue
            y = group['IncLevel1'].values

            for snp in self.snp_data['ID'].values:

                x = group[snp].values
                if len(x.shape) == 2:
                    x = x.T[0]

                yna = y[~np.isnan(x)]
                xna = x[~np.isnan(x)]


                if sum(xna) == 0:
                    continue

                y0snp = yna[xna == 0]
                y1snp = yna[xna == 1]
                y2snp = yna[xna == 2]


                result = [exon_start, exon_end, upes, upee, does, doee, snp]
                result.extend( get_influence_stats(y0snp, y1snp) )
                result.extend( get_influence_stats(y1snp, y2snp) )
                
                
                result.extend([';'.join(["{:.3f}".format(elem) for elem in y0snp]),
                               ';'.join(["{:.3f}".format(elem) for elem in y1snp]),
                               ';'.join(["{:.3f}".format(elem) for elem in y2snp])])

                snp_change.append(result)


        self.snp_correlation[tissue] = pd.DataFrame(snp_change, 
                                                    columns = ['exon_start', 'exon_end', 
                                                               'upes', 'upee', 'does', 'doee',
                                                    'snpID', 
                                                    'diff_0_1', 'pval_0_1', 'len_0', 'len_1_01',
                                                    'diff_1_2', 'pval_1_2', 'len_1', 'len_2', 
                                                    'inc_0', 'inc_1', 'inc_2'])

    def clean_snp_correlation(self, tissue):
        if tissue in self.snp_correlation:
            del self.snp_correlation[tissue]
        
def get_sorted_gene_coordinates(genes_with_exons, promoter_buffer = 1000): 
    
    gene_coordinates = list()
    for gene in genes_with_exons:
        ent = genes_with_exons[gene]
        chromosome, s, e, ori, gene_id = ent.chromosome, ent.gene_start, ent.gene_end, ent.ori, ent.gene_id
        assert s < e
        if ori == '+':
            s -= promoter_buffer
        else:
            e += promoter_buffer
            
        gene_coordinates.append([chromosome, s, e, gene_id])
    gene_coordinates = sorted(gene_coordinates, key = lambda elem: (elem[0], elem[1]))
    
    return gene_coordinates
    
        
def parse_rmats_sample(genes_with_exons, path_to_sample_dir, sample_id, onlychr1 = False, one_gene = False, 
                      gene_subset: list = None):
    
    temp = pd.read_csv(f'{path_to_sample_dir}SE.MATS.JCEC.txt', sep = '\t')
    temp['sample_id'] = sample_id
    
    if onlychr1:
        temp = temp.loc[temp['chr'] == 'chr1']
        
    if one_gene:
        temp = temp.loc[temp['GeneID'] == one_gene]
    
    if gene_subset:
        temp = temp.loc[temp['GeneID'].isin(gene_subset)]
        
        
    
    for gene, group in temp.groupby('GeneID'):
        if gene in genes_with_exons:
            
            for rmats_row in group[['sample_id', 'exonStart_0base', 'exonEnd', 
                   'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE',
                   'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 
                   'IncLevel1']].values:
                genes_with_exons[gene].add_rmats_row_to_buffer(rmats_row)
                

                
                
                