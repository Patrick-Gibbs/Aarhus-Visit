#!/bin/bash
#SBATCH -c 4
#SBATCH --mem 32g
#SBATCH --account dsmwpred
alias ldak='~/snpher/faststorage/ldak5.2.linux'
ldak --linear ./basic_gwas_result.txt --bfile /faststorage/project/dsmwpred/data/ukbb/geno --pheno /faststorage/project/dsmwpred/data/ukbb/height.train --covar /home/patrickgibbs/snpher/faststorage/biobank/genotypes/cauc.sex.town.age.pcs.covar --max-threads 8

