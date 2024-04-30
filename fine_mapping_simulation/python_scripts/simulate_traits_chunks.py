### simulate phenotypes
import os
from bed_reader import open_bed, sample_file
import numpy as np
import sys, argparse
import pandas as pd
import statistics
from scipy.stats import mode

### parsing arguements ###
parser = argparse.ArgumentParser()
 

parser.add_argument("-sp", "--snp-info-out", help = "where infomation about each SNP is to be saved")
parser.add_argument("-b", "--bfile", help = "binary of genomic region")
parser.add_argument("-pp", "--phenotype-path", help = "where the phenotypes are saved")
parser.add_argument("-i", "--simulation-number", help = "used for naming the file, e.g. 10th simulation")
parser.add_argument("-he", "--heritability", help = "SNP Heritability")
parser.add_argument("-nc", "--num-casual", help="the number of SNPs that are casual in the simulation")


# Read arguments from command line
args = parser.parse_args()

bfile = args.bfile 
phenotypes_output_path = args.phenotype_path
snps_info = args.snp_info_out
sample = int(args.simulation_number)
heritability = float(args.heritability)
var_g = heritability # 0.005
num_casual = int(args.num_casual)

# read in genomic regions
ukbb_path = f'{bfile}_{sample}'
bim_df = pd.read_csv(f'{ukbb_path}.bim', delimiter='\t', header = None)
bim_df.columns = ['Chromosome code', 'Variant identifier', 'Position in morgans or centimorgans', 'Base-pair coordinate ', 'Allele 1', 'Allele 2']
fam_df = pd.read_csv(f'{ukbb_path}.fam', delimiter=' ', header = None)

# read plink into python to simulate phenotype

bed = open_bed(f'{ukbb_path}.bed')
X = bed.read()
num_markers = X.shape[1]

# impute NaN to the mode
modes = mode(X, nan_policy='omit').mode
for i in range(X.shape[1]):
    X[np.isnan(X[:,i])] = modes[i]


if num_casual == 1:
    # one snp with effect size one so phenotypes start as being the same as allele state for 1 locus
    casual = np.array([num_markers // 2])
    beta = np.array([1])
    y_g =  X[:,casual] @ beta
else:
    casual = np.random.randint(num_markers//10,num_markers-num_markers//10, size=num_casual)
    beta = np.random.normal(0,2,size=num_casual)
    y_g = X[:,casual] @ beta
    
# add noise and rescale to predermined heritiability

y = (y_g) + np.random.normal(0,statistics.variance(y_g)/heritability - statistics.variance(y_g), size= (len(y_g)))

# save the results
casual_info = bim_df.iloc[casual].copy()
if len(casual) == 1:
    casual_info = pd.DataFrame(casual_info)

casual_info.to_csv(f'{snps_info}/snp_info_sample_{sample}', sep='\t', header = True)
phenotype_df = pd.DataFrame({'acc1': fam_df[0], 'acc2': fam_df[1], 'pheno': y})
phenotype_df.to_csv(f'{phenotypes_output_path}/pheno_{sample}', sep=' ', header=None, index=None)
