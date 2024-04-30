### simulate phenotypes
import os
import numpy as np
import pandas as pd
import sys, argparse, random
from scipy.stats import mode
from bed_reader import open_bed, sample_file
from math import sqrt

### parsing arguements ###
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-n", "--n-casual-snps", help = "The number of casual snps that the trait has")
parser.add_argument("-he", "--trait-heritability", help = "The heritability of the trait i.e. v/v_g")
parser.add_argument("-b", "--bfile", help = "path to bineries")
parser.add_argument("-f", "--output-path", help = "path to output the phenotypes")
parser.add_argument("-i", "--info", help = "path to output infomation e.g. casual snps")
 
# Read arguments from command line
args = parser.parse_args()
n_snps, heritability = int(args.n_casual_snps), float(args.trait_heritability)
bfile, output_path, info_path = args.bfile, args.output_path, args.info
bim = pd.read_csv(f'{bfile}.bim', header=None, sep='\t')
bim.columns = ['Chromosome code', 'Variant identifier', 
               'Position in morgans or centimorgans', 'Base-pair coordinate ', 'Allele 1', 'Allele 2']

fam = pd.read_csv(f'{bfile}.fam', header=None, sep='\t')


casual_snp_indexes = np.random.randint(0, len(bim)-1, size=n_snps)


bed = open_bed(f'{bfile}.bed')
X = bed.read(index=np.s_[:,list(casual_snp_indexes)])

# impute NaN to the mode
modes = mode(X, nan_policy='omit').mode
for i in range(X.shape[1]):
    X[np.isnan(X[:,i])] = modes[i]


beta = np.random.normal(0, 1, n_snps)
y_g = X @ beta
v_g = np.var(y_g)


var_e = 1-heritability
y = (y_g) * (heritability / np.std(y_g)**2) + np.random.normal(0,sqrt(var_e), size= (len(y_g)))

casual_info = bim.iloc[casual_snp_indexes].copy()
casual_info['snp_effect_size'] = beta
casual_info['snp_index_in_chunk'] = casual_snp_indexes

casual_info.to_csv(f'{info_path}', sep='\t', header = True)
print(y)
phenotype_df = pd.DataFrame({'acc1': fam[0], 'acc2': fam[1], 'pheno': y})
phenotype_df.to_csv(f'{output_path}', sep=' ', header=None, index=None)
