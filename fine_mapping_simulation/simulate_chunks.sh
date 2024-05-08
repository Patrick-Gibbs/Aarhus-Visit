######################
######################
# Simulating Genomes #
######################
######################

source simulate_chunks_const.sh

#### SIMULATE TRAIT ####
# sample size is the size of each genomic region
python3 -c "import python_scripts.utils as utils; utils.sample_genome_domains(sample_size=${sample_size}, sample_number=${num_simulations}, genome_path='${ubb_binary_path}', output_path='${genome_sampled_markers_path}/sample')"
python3 -c "import python_scripts.utils as utils; utils.random_indivdual_subset(50000,'/faststorage/project/dsmwpred/data/ukbb/geno.fam', '$subset_individuals', seed=42)"
for ((i=0; i<$num_simulations; i++)); do
plink --bfile $ubb_binary_path --extract ${genome_sampled_markers_path}/sample_${i} --range --make-bed --out ${genome_sample_binary_path}/sample_${i} --keep $subset_individuals 
# make the bed using the choosen genome sample regions
done

for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do

source simulate_chunks_loop_const.sh

ldak='~/snpher/faststorage/ldak5.2.linux'
plink='/home/patrickgibbs/dsmwpred/patrick/plink/plink'

echo """$(cat base_slurm.sh)

#### Setting up files ####
mkdir $experiment_name
mkdir genome_regions_map
for path in $phenotype_path $gwas_sum_stats_path $casual_snp_info; do
    mkdir \$path
done

for ((i=0; i<$num_simulations; i++)); do
python3 -c \"import python_scripts.utils as utils; utils.simulate_traits_for_genome_chunks(bfile='$genome_sample_binary_path/sample', phenotypes_output_path='$phenotype_path', snps_info='$casual_snp_info', heritability=$he, num_casual=$casual_markers, sample=\$i)\"

#### LDAK GWAS####
$ldak --linear $gwas_sum_stats_path/sample_\${i} --bfile ${genome_sample_binary_path}/sample_\${i} --pheno $phenotype_path/pheno_\${i} --max-threads 4

done""" > scripts/job.mk_genome.${he}.${casual_markers}

done; done

mkdir $snps_cors_path

for i in ./scripts/job*; do sbatch -A dsmwpred $i --job-name job.mk_genome.${he}.${casual_markers} --output scripts/std-out; done
mv ./scripts/job* ./scripts/past-jobs
# submit genome sample jobs
python3 -c "import python_scripts.utils as utils; utils.make_ld_matrix(bfile='$genome_sample_binary_path', sample_number=$num_simulations, output_path='$snps_cors_path')"

mkdir $cors_path
# make cors
for ((i=0; i<$num_simulations; i++)); do
    ldak --calc-cors $cors_path/sample_${i} --bfile ${genome_sample_binary_path}/sample_${i} --keep $output_path/rand.5000
done



#######################
#######################
#### FINE  MAPPING ####
#######################
#######################
wget http://dougspeed.com/wp-content/uploads/highld.txt
shuf -n 5000 $subset_individuals > $output_path/rand.5000

for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
mkdir $models_path
echo """$(cat base_slurm.sh)

# fit default LDAK
mkdir $models_path/LDAK_bayesr_default
for ((i=0; i<$num_simulations; i++)); do
    $ldak --mega-prs $models_path/LDAK_bayesr_default/sample_\${i} --model bayesr --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4
done
"""
done; done


for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh

echo """$(cat base_slurm.sh)

# bolt cheating
mkdir $models_path/bolt_cheat
for ((i=0; i<$num_simulations; i++)); do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; $casual_markers/\$num_lines\" | bc)
    echo \"$he 0\$casual_prop 0\" > $models_path/bolt_cheat/params_\$i
    $ldak --mega-prs $models_path/bolt_cheat/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_cheat/params_\$i
done
"""
done;done


# bolt cheat he
mkdir $models_path/bolt_cheat_he
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)

for ((i=0; i<$num_simulations; i++)); do
    for ((as_c=1; as_c<=10; as_c++ )); do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; \$as_c/\$num_lines\" | bc)
    echo \"$he 0\$casual_prop 0\" >> $models_path/bolt_cheat_he/params_\$i
    done
    $ldak --mega-prs $models_path/bolt_cheat_he/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_cheat_he/params_\$i
done
""" > scripts/job.test_models.bolt_cheat_he.${he}.${casual_markers}
done; done

# bolt cheat casual
mkdir $models_path/bolt_cheat_casual
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh

echo """$(cat base_slurm.sh)

for ((i=0; i<$num_simulations; i++)); do
    for he_t in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; $casual_markers/\$num_lines\" | bc)
    echo \"\$he_t 0\$casual_prop 0\" >> $models_path/bolt_cheat_casual/params_\$i
    done
    $ldak --mega-prs $models_path/bolt_cheat_casual/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_cheat_casual/params_\$i
done
""" > scripts/job.test_models.bolt_cheat_casual.${he}.${casual_markers}
done; done

# bolt no cheat
mkdir $models_path/bolt_no_cheat
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)

for ((i=0; i<$num_simulations; i++)); do
    for he_t in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
    for ((as_c=1; as_c<=10; as_c++ )); do
        # make the params file for bolt
        num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
        casual_prop=\$(echo \"scale=7; \$as_c/\$num_lines\" | bc)
        echo \"\$he_t 0\$casual_prop 0\" >> $models_path/bolt_no_cheat/params_\$i
    done; done
    $ldak --mega-prs $models_path/bolt_no_cheat/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_no_cheat/params_\$i   
done
""" > scripts/job.test_models.bolt_no_cheat.${he}.${casual_markers}
done; done


# susieR default
mkdir $models_path/susieR
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)

for ((i=0; i<$num_simulations; i++)); do
    Rscript python_scripts/run_sussie.R $gwas_sum_stats_path/sample_\${i}.summaries $snps_cors_path/sample_\${i} $models_path/susieR/sample_\${i} $number_of_indiduals
done

""" > scripts/job.test_models.susieR.${he}.${casual_markers}
done;done

# susieR default paper settings

for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
mkdir $models_path/susieR_paper_setting
echo """$(cat base_slurm.sh)

for ((i=0; i<$num_simulations; i++)); do
    Rscript python_scripts/run_sussie.R $gwas_sum_stats_path/sample_\${i}.summaries $snps_cors_path/sample_\${i} $models_path/susieR_paper_setting/sample_\${i} $number_of_indiduals 10 1000 TRUE
done

""" > scripts/job.test_models.susieR_paper_setting.${he}.${casual_markers}
done; done



# submit genome sample jobs
for i in ./scripts/job*; do sbatch -A dsmwpred $i --job-name job.mk_genome.${he}.${casual_markers} --output scripts/std-out; done
mv ./scripts/job* ./scripts/past-jobs


#######################
#######################
####   Analyisis   ####
#######################
#######################

for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh

echo """$(cat base_slurm.sh)

for model in \$(find $models_path -mindepth 1 -maxdepth 1 -type d); do
    python3 -c \"import python_scripts.utils as utils; utils.get_stats_chunk_1_snp(pvals_file='$gwas_sum_stats_path', probs_file='\$model', casual_file='$casual_snp_info', model='\$model',sim_numbers=range($num_simulations))\"
done

mkdir $experiment_name/stacked_figures
for ((i=0; i<$num_simulations; i++)); do
python3 -c \"import python_scripts.utils as utils; utils.make_stacked_plots('$gwas_sum_stats_path', '$models_path', '$experiment_name/stacked_figures', '$casual_snp_info', \$i)\"
done

python3 -c \"import python_scripts.utils as utils; utils.get_json_summary('$gwas_sum_stats_path', '$models_path', '$experiment_name/summary.json', '$casual_snp_info', sim_numbers=range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.get_json_summary('$gwas_sum_stats_path', '$models_path', '$experiment_name/summary_ignore_bad.json', '$casual_snp_info', sim_numbers=range($num_simulations), ignore_bad=True)\"
mkdir $experiment_name/summary_plots
python3 -c \"import python_scripts.utils as utils; utils.model_compare_figures('$experiment_name/summary.json', '$experiment_name/summary_plots/summary_plots', statistics=None)\"
python3 -c \"import python_scripts.utils as utils; utils.model_compare_figures('$experiment_name/summary_ignore_bad.json', '$experiment_name/summary_plots/summary_plots_summary_ignore_bad', statistics=None)\"
""" > scripts/job.sum_stats.${he}.${casual_markers}
done; done

# submit genome sample jobs
for i in ./scripts/job*; do sbatch -A dsmwpred $i --job-name job.mk_genome.${he}.${casual_markers} --output scripts/std-out; done
mv ./scripts/job* ./scripts/past-jobs




# makes figures 

    #mkdir \$model/figures
    #for ((i=0; i<$num_simulations; i++)); do
    #python3 -c \"import python_scripts.utils as utils; utils.make_plots_for_genome_region('$gwas_sum_stats_path/sample_\${i}.pvalues', '\$model/sample_\${i}.probs', '\$model/figures/sample_\${i}', '$casual_snp_info/snp_info_sample_\${i}')\"
    #python3 -c \"import python_scripts.utils as utils; utils.make_stacked_plots('$gwas_sum_stats_path', '$models_path', '$experiment_name/stacked_figures', '$casual_snp_info', \$i)\"
    #done; done



