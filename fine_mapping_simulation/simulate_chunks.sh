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
python3 -c "import python_scripts.utils as utils; utils.make_ld_matrix(bfile='$genome_sample_binary_path', sample_number=$num_simulations, output_path='$snps_cors_path')"

# make cors
mkdir $cors_path
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
mkdir $models_path/bayesr_default
echo """$(cat base_slurm.sh)
# fit default LDAK
mkdir $models_path/bayesr_default
for ((i=0; i<$num_simulations; i++)); do
    $ldak --mega-prs $models_path/bayesr_default/sample_\${i} --model bayesr --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4
done
""" > scripts/job.test_models.bayesr_default.${he}.${casual_markers}
done; done

# bolt anti cheat
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/bolt_anti_cheat
for ((i=0; i<$num_simulations; i++)); do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; 10/\$num_lines\" | bc)
    echo \"0.5 0\$casual_prop 0\" > $models_path/bolt_anti_cheat/params_\$i
    $ldak --mega-prs $models_path/bolt_anti_cheat/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_anti_cheat/params_\$i
done
""" > scripts/job.test_models.bolt_anti_cheat.${he}.${casual_markers}
done;done


# bolt cheating
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/bolt_cheat
for ((i=0; i<$num_simulations; i++)); do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; $casual_markers/\$num_lines\" | bc)
    echo \"$he 0\$casual_prop 0\" > $models_path/bolt_cheat/params_\$i
    $ldak --mega-prs $models_path/bolt_cheat/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_cheat/params_\$i
done
""" > scripts/job.test_models.bolt_cheat.${he}.${casual_markers}
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
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
mkdir $models_path/bolt_no_cheat
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

# bolt no cheat
# same as above but power is 0
model_name=bolt_no_cheat_uniform
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
mkdir $models_path/$model_name
echo """$(cat base_slurm.sh)
for ((i=0; i<$num_simulations; i++)); do
    for he_t in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
    for ((as_c=1; as_c<=10; as_c++ )); do
        # make the params file for bolt
        num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
        casual_prop=\$(echo \"scale=7; \$as_c/\$num_lines\" | bc)
        echo \"\$he_t 0\$casual_prop 0\" >> $models_path/$model_name/params_\$i
    done; done
    $ldak --mega-prs $models_path/$model_name/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/$model_name/params_\$i   
done
""" > scripts/job.test_models.$model_name.${he}.${casual_markers}
done; done

# susieR default
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
mkdir $models_path/susieR
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


# susieR cheat casusal
# susieR default paper settings
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
model_name=susie_cheat_casusal
mkdir $models_path/$model_name
echo """$(cat base_slurm.sh)
for ((i=0; i<$num_simulations; i++)); do
    Rscript python_scripts/run_sussie.R $gwas_sum_stats_path/sample_\${i}.summaries $snps_cors_path/sample_\${i} $models_path/$model_name/sample_\${i} $number_of_indiduals $casual_markers 1000 TRUE
done
""" > scripts/job.test_models.$model_name.${he}.${casual_markers}
done; done


# fine map default
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
model_name=finemap
file_name=finemap_file
echo """$(cat base_slurm.sh)
mkdir $experiment_name/finemap_z_files
mkdir $models_path/$model_name
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_z_files(dir='$experiment_name/finemap_z_files', gwas_base='$gwas_sum_stats_path', binary_base='$genome_sample_binary_path', samples=range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
$finemap --sss --in-files $experiment_name/$file_name --n-threads 4
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
""" > scripts/job.test_models.$model_name.${he}.${casual_markers}
done; done

# fine map step wise conditioning
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
model_name=finemap_step_wise
file_name=finemap_step_wise_file
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
#mkdir $experiment_name/finemap_z_files
mkdir $models_path/$model_name
#python3 -c \"import python_scripts.utils as utils; utils.make_finemap_z_files(dir='$experiment_name/finemap_z_files', gwas_base='$gwas_sum_stats_path', binary_base='$genome_sample_binary_path', samples=range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
$finemap --cond --in-files $experiment_name/$file_name
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
""" > scripts/job.test_models.$model_name.${he}.${casual_markers}
done; done


# fine map searching effect sizes
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
model_name=finemap_search_effect_size
file_name=finemap_search_effect_size_file
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
#mkdir $experiment_name/finemap_z_files
mkdir $models_path/$model_name
#python3 -c \"import python_scripts.utils as utils; utils.make_finemap_z_files(dir='$experiment_name/finemap_z_files', gwas_base='$gwas_sum_stats_path', binary_base='$genome_sample_binary_path', samples=range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
$finemap --sss --in-files $experiment_name/$file_name --n-threads 4 --prior-std 0.0001,0.00025,0.0005,0.001,0.0025,0.005,0.01,0.05,0.1,0.5,1
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
""" > scripts/job.test_models.$model_name.${he}.${casual_markers}
done; done


# fine map cheat casual
for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
model_name=finemap_cheat_casual
file_name=finemap_cheat_casual_file
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
#mkdir $experiment_name/finemap_z_files
mkdir $models_path/$model_name
#python3 -c \"import python_scripts.utils as utils; utils.make_finemap_z_files(dir='$experiment_name/finemap_z_files', gwas_base='$gwas_sum_stats_path', binary_base='$genome_sample_binary_path', samples=range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
$finemap --sss --in-files $experiment_name/$file_name --n-threads 4 --n-causal-snps $casual_markers
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
""" > scripts/job.test_models.$model_name.${he}.${casual_markers}
done; done



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
python3 -c 

""" > scripts/job.sum_stats.${he}.${casual_markers}
done; done

for he in 0.02 0.05 0.1; do
for casual_markers in 1 3; do
source simulate_chunks_loop_const.sh
python3 -c "import python_scripts.utils as utils; utils.plot_power_fdr('$models_path', '$casual_snp_info', '$output_path/power_fdr_he_${he}_casual_${casual_markers}.png', range(150), title='he_${he}_casual_${casual_markers}')"
done;done



# makes figures 

    #mkdir \$model/figures
    #for ((i=0; i<$num_simulations; i++)); do
    #python3 -c \"import python_scripts.utils as utils; utils.make_plots_for_genome_region('$gwas_sum_stats_path/sample_\${i}.pvalues', '\$model/sample_\${i}.probs', '\$model/figures/sample_\${i}', '$casual_snp_info/snp_info_sample_\${i}')\"
    #python3 -c \"import python_scripts.utils as utils; utils.make_stacked_plots('$gwas_sum_stats_path', '$models_path', '$experiment_name/stacked_figures', '$casual_snp_info', \$i)\"
    #done; done



squeue -u $USER -h | awk '{print $1}' | xargs scancel