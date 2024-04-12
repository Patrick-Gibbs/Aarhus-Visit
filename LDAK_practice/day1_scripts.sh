
# basic gwas
# this is not corrected for population structure meaning that we might
# find inflated  peaks of assiocation which are correlated with eachother

mkdir basic_gwas
cd basic_gwas

# the gwas
ldak --linear ./basic_gwas_result.txt --bfile /faststorage/project/dsmwpred/data/ukbb/geno --pheno /faststorage/project/dsmwpred/data/ukbb/height.train --max-threads 8

# taking a look at the manhattan plot using python
python3
```
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np

    df = pd.read_csv('basic_gwas_result.txt.pvalues', sep=' ')
    p_vals = -np.log(np.array(df['P']))

    marker_chrom = [int(loci[:loci.index(':')]) for loci in  df['Predictor']]
    colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

    plt.scatter(np.arange(len(p_vals)),p_vals, c=colors)

    # corrected p val
    plt.axhline(y=-np.log(0.05/len(p_vals)), color='black', linestyle='-')

    plt.savefig('gwas_height_basic.png')
```
# inflated assioactions it seems from figure


# correcting for population structure/kinship
# to fix this problem
cd ..
mkdir kinship_gwas
cd kinship_gwas

# making kinship matrix
ldak --thin le --bfile /faststorage/project/dsmwpred/data/ukbb/geno --window-prune .05 --window-cm 1 --keep rand.5000

# having trouble here with memory, might be because square kinship matrix is large
# choosing to move on for now
ldak --calc-kins-direct le --bfile /faststorage/project/dsmwpred/data/ukbb/geno --power -1 --extract le.in
ldak --linear single4 --bfile human --pheno quant.pheno --grm le

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

# with validation dataset
ldak --validate prs_ridge_height.val --bfile /faststorage/project/dsmwpred/data/ukbb/geno --pheno /faststorage/project/dsmwpred/data/ukbb/height.test --scorefile prs_ridge_height.effects --max-threads 8
ldak --validate prs_bolt_height.val --bfile /faststorage/project/dsmwpred/data/ukbb/geno --pheno /faststorage/project/dsmwpred/data/ukbb/height.test --scorefile prs_bolt_height.effects --max-threads 8
