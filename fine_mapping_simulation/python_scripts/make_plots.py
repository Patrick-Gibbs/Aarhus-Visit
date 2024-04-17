import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


### parsing arguements ###
parser = argparse.ArgumentParser()
 
# Adding optional argument
parser.add_argument("-s", "--sample", help = "simulation number for plots")
 
# Read arguments from command line
args = parser.parse_args()
sample = args.sample 
real_snp_location_df = pd.read_csv(f'one_marker_phenotypes/casual_snp_info/snp_info_sample_{sample}', sep='\t',header=None)
casual_i = int(list(real_snp_location_df.set_index(0).loc['snp_index_in_chunk'])[0])

gwas_path = f'fine_map_gwas/gwas_sample_{sample}.pvalues'
bayesr_path = f'fine_map_bayesr/bayesr_sample_{sample}.probs'
output_path = f'figures/sample_{sample}'

df = pd.read_csv(bayesr_path, sep=' ')
probs = np.array(df['Probability'])

marker_chrom = [int(loci[:loci.index(':')]) for loci in df['Predictor']]
colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

#### make a figure with GWAS and Post.Prob ####

# Post. Prob Plot 
fig, axs = plt.subplots(2, 1)
axs[0].scatter(np.arange(len(probs))[casual_i], probs[casual_i], s=casual_size, marker = '*', c='red')
axs[0].scatter(np.delete(np.arange(len(probs)), casual_i), np.delete(probs, casual_i), s=size)

axs[0].set_ylabel("Post. Prob")


# Basic GWAS
p_critical = 5*10**(-8)
df = pd.read_csv(gwas_path, sep=' ')
p_vals = -np.log(np.array(df['P']))

inf_positions = np.isinf(p_vals)

# Step 2: Find the maximum of non-infinite elements
max_value = np.max(p_vals[~inf_positions])

# Step 3: Replace inf with max_value + 100
p_vals[inf_positions] = max_value + 100

axs[1].scatter([casual_i],p_vals[casual_i], s=casual_size, marker = '*', c='red')
axs[1].scatter(np.delete(np.arange(len(probs)), casual_i),np.delete(p_vals, casual_i), s=size)
axs[1].axhline(y=-np.log(p_critical), color='black', linestyle='-')

axs[1].set_xlabel("Marker")
axs[1].set_ylabel("-Log p wald")
plt.savefig(f'{output_path}.png')