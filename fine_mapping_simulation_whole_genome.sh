#ldak --calc-cors whole_genome/cor --bfile /faststorage/project/dsmwpred/data/ukbb/geno --keep ../LDAK_practice/gwas_kinship/rand.5000
#wget http://dougspeed.com/wp-content/uploads/highld.txt

# simulate trait with 50 markers corosponding to the whole genome # n is number of snps and he is trait heritability
for i in {1..2} do;
python3 python_scripts/simulate_genome_wide_trait.py  -n 50 -he 0.15 --bfile /faststorage/project/dsmwpred/data/ukbb/geno --output-path whole_genome/phenotypes/phenotypes/pheno_${i} -i whole_genome/phenotypes/infomation/info_${i}
ldak --linear whole_genome/gwas/gwas_${i} --bfile /faststorage/project/dsmwpred/data/ukbb/geno --pheno whole_genome/phenotypes/phenotypes/pheno_$i --max-threads 4
ldak --mega-prs whole_genome/bayesr/bayesr_${i} --model bayesr --summary whole_genome/gwas/gwas_${i}.summaries --cors whole_genome/cors --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4
mkdir whole_genome/figures/simulation_${i}
python3 python_scripts/view_regions.py --p-value whole_genome/gwas/gwas_${i}.pvalues --posterior whole_genome/bayesr/bayesr_${i}.probs --output-path whole_genome/figures/simulation_${i}/whole_genome --info whole_genome/phenotypes/infomation/info_${i}
done