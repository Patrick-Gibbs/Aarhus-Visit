
#### Constants #####
num_simulations=1
random_seed=42
number_of_indiduals=50000 # max is 392214
sample_size=2000000
output_path=small_chunks
###################

for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do

experiment_name=${output_path}/chucks_1_marker_her_${he}_ind${number_of_indiduals}_casual_${casual_markers}
genome_sampled_marks_dir=$experiment_name/genome_regions_map
genome_sample_binary_dir=$experiment_name/genome_sample
phenotype_path=$experiment_name/phenotypes_casual
gwas_sum_stats_path=$experiment_name/fine_map_gwas
casual_snp_info=$experiment_name/casual_snp_info

genome_sampled_markers_path=$genome_sampled_marks_dir/genome_sample
genome_sample_binary_path=$genome_sample_binary_dir/genome_sample

python -c "import python_scripts.utils as utils; utils.random_indivdual_subset(50000,42,'/faststorage/project/dsmwpred/data/ukbb/geno.fam', 'rand50000.fam')"
ldak='~/snpher/faststorage/ldak5.2.linux'
plink='/home/patrickgibbs/dsmwpred/patrick/plink/plink'

echo """#"'!'"/bin/bash
#SBATCH --mem 32G
#SBATCH -c 4
#SBATCH -t 00:12:00

#### Setting up files ####
mkdir $experiment_name
mkdir genome_regions_map
for path in $genome_sampled_marks_dir $genome_sample_binary_dir $phenotype_path $gwas_sum_stats_path $casual_snp_info; do
    mkdir \$path
done


#### SIMULATE TRAIT ####
ubb_binary_path=/faststorage/project/dsmwpred/data/ukbb/geno
python3 -c \"import python_scripts.utils as utils; utils.random_indivdual_subset($number_of_indiduals,42,'/faststorage/project/dsmwpred/data/ukbb/geno.fam', 'rand${number_of_indiduals}.fam')\"

# sample size is the size of each genomic region
python3 -c \"import python_scripts.utils as utils; utils.sample_genome_domains(sample_size=${sample_size}, sample_number=${num_simulations}, genome_path=${ubb_binary_path}, output_path=${genome_sampled_markers_path})\"

for ((i=0; i<$num_simulations; i++)); do
# make the bed using the choosen genome sample regions
$plink --bfile $ubb_binary_path --extract ${genome_sampled_markers_path}_\${i} --range --make-bed --out ${genome_sample_binary_path}_\${i} --keep rand${number_of_indiduals}.fam 

# simulate traits
#python3 python_scripts/simulate_traits_chunks.py --b $genome_sample_binary_path -i \$i -he $he --phenotype-path $phenotype_path --snp-info-out $casual_snp_info --num-casual $casual_markers
python3 -c \"import python_scripts.utils as utils; utils.simulate_traits_for_genome_chunks(bfile=$genome_sample_binary_path, phenotypes_output_path=$phenotype_path, snps_info=$casual_snp_info, heritability=$he, num_casual=$casual_markers, sample=\$i)\"


#### LDAK GWAS####
$ldak --linear $gwas_sum_stats_path/gwas_sample_\${i} --bfile ${genome_sample_binary_path}_\${i} --pheno $phenotype_path/pheno_\${i} --max-threads 4

done""" > scripts/job.mk_genome.${he}.${casual_markers}

done; done


for i in ./scripts/job*; do sbatch -A dsmwpred $i --job-name job.mk_genome.${he}.${casual_markers} --#SBATCH --output scripts/std-out; done
mv ./scripts/job* mv ./scripts/past-jobs 



