import pandas as pd
from random import randint
from tqdm import tqdm
import os, sys, getopt, argparse
CHROMOSOME = 0
START = 1
FINISH = 2
NUM_HUMAN_AUTOSOMES = 22

### parsing arguements ###
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-s", "--sample-size", help = "the number of base pairs to sample")
parser.add_argument("-b", "--bfile", help = "path to genome binary")
parser.add_argument("-o", "--output-path", help = "locatation to write the output")
 
# Read arguments from command line
args = parser.parse_args()
sample_size, genome_path, output_path = args.sample_size, args.bfile, args.output_path 

### sample genome ###
ukbb_bim = pd.read_csv("/faststorage/project/dsmwpred/data/ukbb/geno.bim", sep='\t', header=None)

# making another column containing the position of each SNP on  a chromosome as an integer
ukbb_bim['chrom_pos'] = [int(item[item.index(':') + 1:item.index('_')]) for item in ukbb_bim[1]]
ukbb_bim.rename(columns = {0:'chrom_number'}, inplace = True)

approx_chromosome_size = [max(ukbb_bim[ukbb_bim['chrom_number'] == chrom]['chrom_pos']) for chrom in range(1,NUM_HUMAN_AUTOSOMES + 1)]

def is_overlap(p_bounds, c_bound):
    for p_bound in p_bounds:
        if (c_bound[CHROMOSOME] == p_bound[CHROMOSOME]
            and pd.Interval(p_bound[START], p_bound[FINISH]).overlaps(pd.Interval(c_bound[START], c_bound[FINISH]))):
            return True
    return False 

def sample_snp(previous_bounds, size = 10**6):
    """
    Selects a random genome location of `size`. Also ensures that selected region has not before been selected
    and is thus not in previous_bounds
    """

    # selected a random snp from the bim
    snp = ukbb_bim.iloc[randint(0,len(ukbb_bim))]
    
    # extract snp infomation
    snp_chromosome = snp['chrom_number']
    snp_pos = snp['chrom_pos']

    # extract += 1mb from the randomly seleted SNP
    start = int(snp_pos - size/2)
    end = int(snp_pos + size/2 - 1)


    # if the the bound is outside the chromosome try again
    if (end < 0 or start > approx_chromosome_size[snp_chromosome - 1] 
             or is_overlap(previous_bounds, (snp_chromosome, start, end))):
        return sample_snp(bounds)
    

    return snp_chromosome, start, end
    

# select 200 1mb regions
bounds = []
for i in tqdm(range(200)):
    bounds.append(sample_snp(bounds))

    with open(f'genome_regions_map/genome_sample_{i}', 'w') as f:
        f.write(f'{bounds[i][CHROMOSOME]} {bounds[i][START]} {bounds[i][FINISH]} sample_{i}')
   