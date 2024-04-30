mkdir fine_mapping
cd fine_mapping

wget http://dougspeed.com/wp-content/uploads/highld.txt
ldak --calc-cors cors --bfile /faststorage/project/dsmwpred/data/ukbb/geno --keep ../LDAK_practice/gwas_kinship/rand.5000
ldak --mega-prs prs_bayesr_height --model bayesr --summary ../LDAK_practice/basic_gwas/basic_gwas_result.txt.summaries --cors cors --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4


### plotting in python ####
python3 << EOF
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv('prs_bayesr_height.probs', sep=' ')
probs = np.array(df['Probability'])

marker_chrom = [int(loci[:loci.index(':')]) for loci in df['Predictor']]
colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

#### make a figure with GWAS and Post.Prob ####

# Post. Prob Plot 
fig, axs = plt.subplots(2, 1)
axs[0].scatter(np.arange(len(probs)), probs, c=colors, s=0.5)
axs[0].set_ylabel("Post. Prob")


# Basic GWAS
p_critical = 5*10**(-8)
df = pd.read_csv('../LDAK_practice/basic_gwas/basic_gwas_result.txt.pvalues', sep=' ')
p_vals = -np.log(np.array(df['P']))
marker_chrom = [int(loci[:loci.index(':')]) for loci in  df['Predictor']]
colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

axs[1].scatter(np.arange(len(p_vals)),p_vals, c=colors, s=0.5)
axs[1].axhline(y=-np.log(p_critical), color='black', linestyle='-')
axs[1].set_xlabel("Marker")
axs[1].set_ylabel("-Log p wald")
plt.savefig('BayesR_Post.png')
EOF