
from utils import *
import argparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""takes info file and makes a figure of each casual variant."""

### parsing arguements ###
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-p", "--posterior", help = "The heritability of the trait i.e. v/v_g")
parser.add_argument("-g", "--p-value", help = "path to bineries")
parser.add_argument("-f", "--output-path", help = "path to output the phenotypes")
parser.add_argument("-i", "--info", help = "path to output infomation e.g. casual snps")


# Read arguments from command line
args = parser.parse_args()

gwas_p_path, probs_path  = args.p_value, args.posterior
output_path, info_path = args.output_path, args.info

info = pd.read_csv(info_path, sep='\t')
markers_examined = list(info['snp_index_in_chunk'])

for marker in sorted(markers_examined):
    plot_region(marker, gwas_p_path, probs_path, f'{output_path}_{marker}', markers_examined = 1000, point_size = 3.5)


