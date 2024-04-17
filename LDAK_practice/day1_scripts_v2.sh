
# basic gwas
# this is not corrected for population structure meaning that we might
# find inflated  peaks of assiocation which are correlated with eachother

mkdir basic_gwas
cd basic_gwas

# the gwas
ldak --linear ./basic_gwas_result.txt --bfile /faststorage/project/dsmwpred/data/ukbb/geno --pheno /faststorage/project/dsmwpred/data/ukbb/height.train --max-threads 8

# taking a look at the manhattan plot using python
python3 << EOF

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

p_critical = 5*10**(-8)

df = pd.read_csv('basic_gwas_result.txt.pvalues', sep=' ')
p_vals = -np.log(np.array(df['P']))

marker_chrom = [int(loci[:loci.index(':')]) for loci in  df['Predictor']]
colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

plt.scatter(np.arange(len(p_vals)),p_vals, c=colors, s=0.5)

# corrected p val
plt.axhline(y=-np.log(p_critical), color='black', linestyle='-')

plt.savefig('gwas_height_basic.png')

EOF
# inflated assioactions it seems from figure


# correcting for population structure/kinship
# to fix this problem
cd ..
ldak --linear ./corrected_gwas_result.txt --bfile /faststorage/project/dsmwpred/data/ukbb/geno --pheno /faststorage/project/dsmwpred/data/ukbb/height.train --covar /home/patrickgibbs/snpher/faststorage/biobank/genotypes/cauc.sex.town.age.pcs.covar --max-threads 4

# reuse python3 script above for plot


# now generating PRS using GWAS summary statistics for ridge and bolt.
# these represent two linear models assuming different prior distributions
# of effect sizes

# I have not used summary statistics for a model before. I would normally
# simply fit to markers directly, I assume using summary statistics
# is computationally faster. Nevertheless, the rationale is not yet completely
# clear to me

# finding correlated predictors, I assume we are doing this due to helps with
# shrinkage / is used in optimization, this is also a little unfamiliar as
# I assumed penalized models automatically apply shrinkage between correlated
# predictors implicitly
ldak --calc-cors cors --bfile /faststorage/project/dsmwpred/data/ukbb/geno --keep ../gwas_kinship/rand.5000

# the software recommended to include regions in high LD, I assume these are
# treated differently by the model, the exact details of which I am so far
# unsure.
wget http://dougspeed.com/wp-content/uploads/highld.txt
ldak --mega-prs prs_ridge_height --model ridge --summary ../basic_gwas/basic_gwas_result.txt.summaries --cors cors --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 8
ldak --mega-prs prs_bolt_height --model bolt --summary ../basic_gwas/basic_gwas_result.txt.summaries --cors cors --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 8


# use test set to get held out r2 values
ldak --calc-scores scores_ridge_height --scorefile prs_ridge_height.effects --bfile /faststorage/project/dsmwpred/data/ukbb/geno --power 0
ldak --calc-scores scores_bolt_height --scorefile prs_bolt_height.effects --bfile /faststorage/project/dsmwpred/data/ukbb/geno --power 0

# find r2 score with python
# (there might be a bash way to do this which could be more elegant)
python3 << EOF

from sklearn.metrics import r2_score
import pandas as pd

Y_PHENO_COL = 2
Y_HAT_PHENO_COL = 'Profile_1'
Y_HAT_ID =  'ID1'
Y_ID = 0

def get_score(pheno_y_hat):
    pheno_true_path = '/faststorage/project/dsmwpred/data/ukbb/height.train'

    y = pd.read_csv(pheno_true_path, header=None, sep = ' ')
    y_hat = pd.read_csv(pheno_y_hat, sep = '\t')
    y_hat = y_hat[y_hat[Y_HAT_ID].isin(set(y[Y_ID]))]

    return (r2_score(y[Y_PHENO_COL], y_hat[Y_HAT_PHENO_COL]))

print('ridge r2')
print(get_score('scores_ridge_height.profile'))
print('bolt r2')
print(get_score('scores_bolt_height.profile'))
EOF
