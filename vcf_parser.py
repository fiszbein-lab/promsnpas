import gzip
import subprocess

import numpy as np
import pandas as pd



def get_column_names(full_path_to_vcf):

    vcf_columns = list()
    with gzip.open(full_path_to_vcf,'r') as fin:        
        for line in fin:  
            decoded = line.decode()
            if decoded.startswith('#CHROM'):
                vcf_columns = decoded[1:].strip().split()
                break
                
    return vcf_columns



def collect_snps_from_vcf(full_path_to_vcf, chromosome = 'chr5', region_start = 1295038, region_end = 1295289):
    
    vcf_columns = get_column_names(full_path_to_vcf)
    
    cmd = ['tabix', 
           full_path_to_vcf, 
           f'{chromosome}:{region_start}-{region_end}']
    
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    snp_lists = list()
    
    for elem in stdout.decode().split('\n'):
        one_snp = elem.split('\t')
        if len(one_snp) == 1:
            continue
            
        snp_lists.append(one_snp)
        
    return snp_lists, vcf_columns


string_to_copies_dict = {'0|0': 0, '0|1': 1, '1|0': 1, '1|1': 2, '.|.': 0}

def snp_lists_to_df(snp_lists, vcf_columns):
    snp_data = pd.DataFrame(snp_lists, columns = vcf_columns)
    snp_data['POS'] = snp_data['POS'].astype(int)
    
    sample_columns = list()

    for col in vcf_columns:
        if col.startswith('GTEX'):
            sample_columns.append(col)

    if len(snp_data) == 0:
        return snp_data

    snp_data[sample_columns] = np.array(
        [[string_to_copies_dict[v] for v in varr] 
         for varr in snp_data[sample_columns].values])
    
    return snp_data


