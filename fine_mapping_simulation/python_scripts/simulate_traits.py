### simulate phenotypes
import os
from pandas_plink import read_plink
import numpy as np
import sys, argparse


### parsing arguements ###
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-s", "--start", help = "start index")
parser.add_argument("-f", "--finish", help = "end index")
 
# Read arguments from command line
args = parser.parse_args()
start, finish = args.start, args.finish


snp_info_output_path = '/home/patrickgibbs/dsmwpred/patrick/Patrick_Aarhus/Aarhus-Visit/fine_mapping_simulation/one_marker_phenotypes/casual_snp_info'
phenotypes_output_path = '/home/patrickgibbs/dsmwpred/patrick/Patrick_Aarhus/Aarhus-Visit/fine_mapping_simulation/one_marker_phenotypes/phenotypes'
var_g = 0.005

for sample in range(200):

    # read in genomic regions
    ukbb_path = f'/home/patrickgibbs/dsmwpred/patrick/Patrick_Aarhus/Aarhus-Visit/fine_mapping_simulation/genome_regions/genome_sample_{sample}'
    bim_df = pd.read_csv(f'{ukbb_path}.bim', delimiter='\t', header = None)
    bim_df.columns = ['Chromosome code', 'Variant identifier', 'Position in morgans or centimorgans', 'Base-pair coordinate ', 'Allele 1', 'Allele 2']
    fam_df = pd.read_csv(f'{ukbb_path}.fam', delimiter=' ', header = None)

    # read plink into python to simulate phenotype
    (bim, fam, bed) = read_plink(ukbb_path, verbose=False)
    X = bed.compute().T
    num_markers = X.shape[1]

    # such that the casual snp will never have a missing values across all individuals
    non_nan_columns = ~np.any(np.isnan(X), axis=0)
    def get_casual(snp_range):
        casual = randint(0,snp_range)
        if not non_nan_columns[casual]:
            return get_casual(snp_range)
        return casual
    casual = get_casual(num_markers)

    # one snp with effect size one so phenotypes start as being the same as allele state for 1 locus
    y_g = 1 * X[:,casual] 

    # add noise and rescale to predermined heritiability
    var_e = 1 - var_g
    y = (y_g) * (0.05 / np.std(y_g)**2) + np.random.normal(0,var_e, size= (len(y_g)))

    # save the results
    casual_info = bim_df.iloc[casual].copy()
    casual_info['heritability'] = var_g
    casual_info['snp_index_in_chunk'] = casual
    casual_info = pd.DataFrame(casual_info)
    casual_info.to_csv(f'{snp_info_output_path}/snp_info_sample_{sample}', sep='\t', header = False)

    phenotype_df = pd.DataFrame({'acc1': fam_df[0], 'acc2': fam_df[1], 'pheno': y})
    phenotype_df.to_csv(f'{phenotypes_output_path}/pheno_{sample}', sep=' ', header=None, index=None)
