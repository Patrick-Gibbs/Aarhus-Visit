####     Example Simulation Variables     ####
#### which are usually set using python constructor ####
num_simulations=1
he=0.1
random_seed=42
casual_markers=3
number_of_indiduals=50000 # max is 392214
sample_size=2000000
model=bolt
default_model=1
model_label=bolt_cheat
params_path=params
############################


################################
# Code to Run Simulation Below #
################################

# alias
alias ldak='~/snpher/faststorage/ldak5.2.linux'
alias plink='/home/patrickgibbs/dsmwpred/patrick/plink/plink'

#### Setting up files ####
experiment_name=small_chunks/chucks_her_${he}_ind${number_of_indiduals}_casual_${casual_markers}
genome_sampled_markers_path=$experiment_name/genome_regions_map/genome_sample
genome_sample_binary_path=$experiment_name/genome_regions/genome_sample
phenotype_path=$experiment_name/phenotypes_casual
gwas_sum_stats_path=$experiment_name/fine_map_gwas
bayesr_output_path=$experiment_name/$model_label/model
figure_output_class=$experiment_name/$model_label/figures
casual_snp_info=$experiment_name/casual_snp_info
cors=$experiment_name/cors
mkdir $experiment_name
mkdir $experiment_name/$model_label
for path in $genome_sampled_markers_path $genome_sample_binary_path $phenotype_path $gwas_sum_stats_path $bayesr_output_path $figure_output_class $casual_snp_info $cors; do
    mkdir $path
done


###########################
#### SIMULATE TRAIT ####
ubb_binary_path=/faststorage/project/dsmwpred/data/ukbb/geno
python -c "import python_scripts.utils as utils; utils.random_indivdual_subset($number_of_indiduals,42,'/faststorage/project/dsmwpred/data/ukbb/geno.fam', 'rand${number_of_indiduals}.fam')"

# sample size is the size of each genomic region
python3 python_scripts/sample_genome_domains.py --sample-size $sample_size --bfile $ubb_binary_path --output-path ${genome_sampled_markers_path} --sample-number $num_simulations 

for ((i=0; i<$num_simulations; i++)); do
# make the bed using the choosen genome sample regions
plink --bfile $ubb_binary_path --extract ${genome_sampled_markers_path}_${i} --range --make-bed --out ${genome_sample_binary_path}_${i} --keep rand${number_of_indiduals}.fam 

# simulate traits
python3 python_scripts/simulate_traits_chunks.py --b $genome_sample_binary_path -i $i -he $he --phenotype-path $phenotype_path --snp-info-out $casual_snp_info --num-casual $casual_markers
#########################


#### LDAK GWAS####
ldak --linear $gwas_sum_stats_path/gwas_sample_${i} --bfile ${genome_sample_binary_path}_${i} --pheno $phenotype_path/pheno_${i} --max-threads 4
##################

















#### LDAK FINE MAP ####
# wget http://dougspeed.com/wp-content/uploads/highld.txt
# fit bayesR mega PRS to each
ldak --calc-cors $cors/cor_sample_${i} --bfile ${genome_sample_binary_path}_${i} --keep rand.5000

if [[ $default_model -eq 0 ]]; then
    params=()
else
    params=("--parameters" "${params_path}")
fi

ldak --mega-prs $bayesr_output_path/${model_label}_sample_${i} --model $model --summary $gwas_sum_stats_path/gwas_sample_${i}.summaries --cors $cors/cor_sample_${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 "${params[@]}"
#######################


#### Figures for each simulation ####
python -c "import python_scripts.utils as utils; utils.make_plots_for_genome_region('$gwas_sum_stats_path/gwas_sample_${i}.pvalues', '$bayesr_output_path/${model_label}_sample_${i}.probs', '$figure_output_class/sample_${i}', '$casual_snp_info/snp_info_sample_${i}')"
done
#####################################



#### Generating Statisitcs ####
python -c "import python_scripts.utils as utils; utils.get_stats_chunk_1_snp('$experiment_name', None,range($num_simulations))"
python -c "import python_scripts.utils as utils; utils.get_stats_chunk_1_snp('$experiment_name', None,range($num_simulations), ignore_bad=True)"
################################

