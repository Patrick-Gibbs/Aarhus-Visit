######################
######################
# Simulating Genomes #
######################
######################

source simulate_chunks_const.sh
echo """$(cat base_slurm.sh)
source simulate_chunks_const.sh
mkdir \$output_path

mkdir \$genome_sample_binary_path
mkdir \$genome_sampled_markers_path
mkdir \$chunked_simulations
mkdir \$combined_simulations
mkdir \$chunked_compheno_simulations
mkdir \$combined_samples
mkdir \$snps_cors_path
mkdir \$cors_path
for he in \$HE_VALUES; do
for casual_markers in \$CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
mkdir \$experiment_name
for path in \$phenotype_path \$gwas_sum_stats_path \$casual_snp_info \$gwas_sum_stats_path_combined_pheno \$experiment_name_chunk_compheno \$gwas_chunked_combine_pheno \$experiment_name_combined \$gwas_sum_stats_combined_path; do
    mkdir \$path
done
done; done
 
#### SIMULATE TRAIT ####
# sample size is the size of each genomic region
python3 -c \"import python_scripts.utils as utils; utils.sample_genome_domains(sample_size=\${sample_size}, sample_number=\${num_simulations}, genome_path='\${ukbb_binary_path}', output_path='\${genome_sampled_markers_path}/sample')\"
python3 -c \"import python_scripts.utils as utils; utils.random_indivdual_subset(50000,'/faststorage/project/dsmwpred/data/ukbb/geno.fam', '\$subset_individuals', seed=42)\"
for ((i=0; i<\$num_simulations; i++)); do
$plink --bfile \$ukbb_binary_path --extract \${genome_sampled_markers_path}/sample_\${i} --range --make-bed --out \${genome_sample_binary_path}/sample_\${i} --keep \$subset_individuals 
# make the bed using the choosen genome sample regions
done

# make genome combine
source simulate_chunks_const.sh
rm \$combined_samples_regions
touch \$combined_samples_regions
for ((i=0; i<\$num_simulations; i++)); do
    cat \$genome_sampled_markers_path/sample_\${i} >> \$combined_samples_regions
    echo '' >> \$combined_samples_regions
done
$plink --bfile \$ukbb_binary_path --extract \$combined_samples_regions --range --make-bed --out \$combined_samples/combined_samples --keep \$subset_individuals 

for he in \$HE_VALUES; do
for casual_markers in \$CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
# simulate phenotypes
for ((i=0; i<\$num_simulations; i++)); do
python3 -c \"import python_scripts.utils as utils; utils.simulate_traits_for_genome_chunks(bfile='\$genome_sample_binary_path/sample', phenotypes_output_path='\$phenotype_path', snps_info='\$casual_snp_info', heritability=\$he, num_casual=\$casual_markers, sample=\$i)\"
done; done; done

# make average phenotypes
for he in \$HE_VALUES; do
for casual_markers in \$CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
python3 -c \"import python_scripts.utils as utils; utils.combine_phenotypes('\$phenotype_path', '\$combined_pheno', range(\$num_simulations))\"
done; done

shuf -n 5000 $subset_individuals > $output_path/rand.5000
""" > scripts/job.0.simulate_chunks


#### LDAK GWAS####
source simulate_chunks_const.sh
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh

echo """$(cat base_slurm.sh)
for ((i=0; i<$num_simulations; i++)); do
$ldak --linear $gwas_sum_stats_path/sample_\${i} --bfile ${genome_sample_binary_path}/sample_\${i} --pheno $phenotype_path/sample_\${i} --max-threads 4
$ldak --linear $gwas_chunked_combine_pheno/sample_\${i} --bfile ${genome_sample_binary_path}/sample_\${i} --pheno $combined_pheno --max-threads 4
done

$ldak --linear $gwas_sum_stats_combined_path/gwas --bfile $combined_samples/combined_samples --pheno $combined_pheno --max-threads 4
""" > scripts/job.1.gwas_combined.${he}_${casual_markers}
done;done

# make snp ld matrix for sussie and finemap
echo """$(cat base_slurm_small.sh)
python3 -c \"import python_scripts.utils as utils; utils.make_ld_matrix(bfile='$genome_sample_binary_path', sample_number=$num_simulations, output_path='$snps_cors_path')\"
""" > scripts/job.1.make_ld_matrix.${he}_${casual_markers}



# make cors for ldak
echo """$(cat base_slurm.sh)
for ((i=0; i<$num_simulations; i++)); do
    $ldak --calc-cors $cors_path/sample_\${i} --bfile ${genome_sample_binary_path}/sample_\${i} --keep $output_path/rand.5000
done
""" > scripts/job.1.make_cors_ldak


# make cors for ldak combined
echo """$(cat base_slurm.sh)
$ldak --calc-cors $combined_samples/cors --bfile $combined_samples/combined_samples --pheno $combined_pheno --max-threads 4 --keep $output_path/rand.5000

""" > scripts/job.1.make_cors_ldak_combined

echo """$(cat base_slurm_small.sh)
source simulate_chunks_const.sh
for he in \$HE_VALUES; do
for casual_markers in \$CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
mkdir \$models_path
mkdir \$combine_models
mkdir \$chunked_compheno_models
mkdir \$experiment_name/finemap_z_files
mkdir \$experiment_name_chunk_compheno/finemap_z_files
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_z_files(dir='\$experiment_name/finemap_z_files', gwas_base='\$gwas_sum_stats_path', binary_base='\$genome_sample_binary_path', samples=range(\$num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_z_files(dir='\$experiment_name_chunk_compheno/finemap_z_files', gwas_base='\$gwas_chunked_combine_pheno', binary_base='\$genome_sample_binary_path', samples=range(\$num_simulations))\"
done; done 

""" > scripts/job.2.sim_dirs

#######################
#######################
#### FINE  MAPPING ####
#######################
#######################

wget http://dougspeed.com/wp-content/uploads/highld.txt
source simulate_chunks_const.sh


# fit default LDAK
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/bayesr_default
for ((i=0; i<$num_simulations; i++)); do
    $ldak --mega-prs $models_path/bayesr_default/sample_\${i} --model bayesr --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4
done
""" > scripts/job.3.test_models.bayesr_default.${he}.${casual_markers}
done; done

# bolt anti cheat
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
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
""" > scripts/job.3.test_models.bolt_anti_cheat.${he}.${casual_markers}
done;done

# bolt cheating
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/bolt_cheat
mkdir $chunked_compheno_models/bolt_cheat
for ((i=0; i<$num_simulations; i++)); do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; $casual_markers/\$num_lines\" | bc)
    echo \"$he 0\$casual_prop 0\" > $models_path/bolt_cheat/params_\$i
    $ldak --mega-prs $models_path/bolt_cheat/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_cheat/params_\$i
    $ldak --mega-prs $chunked_compheno_models/bolt_cheat/sample_\${i} --model bolt --summary $gwas_chunked_combine_pheno/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_cheat/params_\$i
done
""" > scripts/job.3.test_models.bolt_cheat.${he}.${casual_markers}
done;done

# bolt cheat he
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
model_name=bolt_cheat_he
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/$model_name
rm $models_path/$model_name/params_\$i
touch $models_path/$model_name/params_\$i
for ((i=0; i<$num_simulations; i++)); do
    for ((as_c=1; as_c<=10; as_c++ )); do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; \$as_c/\$num_lines\" | bc)
    echo \"$he 0\$casual_prop 0\" >> $models_path/$model_name/params_\$i
    done
    $ldak --mega-prs $models_path/$model_name/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/$model_name/params_\$i
done
""" > scripts/job.3.test_models.$model_name.${he}.${casual_markers}
done; done

# bolt cheat casual
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/bolt_cheat_casual
for ((i=0; i<$num_simulations; i++)); do
    for he_t in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
    # make the params file for bolt
    num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
    casual_prop=\$(echo \"scale=7; $casual_markers/\$num_lines\" | bc)
    echo \"\$he_t 0\$casual_prop 0\" >> $models_path/bolt_cheat_casual/params_\$i
    done
    $ldak --mega-prs $models_path/bolt_cheat_casual/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_cheat_casual/params_\$i
done
""" > scripts/job.3.test_models.bolt_cheat_casual.${he}.${casual_markers}
done; done

# bolt no cheat
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/bolt_no_cheat
mkdir $chunked_compheno_models/bolt_no_cheat
for ((i=0; i<$num_simulations; i++)); do
    for he_t in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
    for ((as_c=1; as_c<=10; as_c++ )); do
        # make the params file for bolt
        num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
        casual_prop=\$(echo \"scale=7; \$as_c/\$num_lines\" | bc)
        echo \"\$he_t 0\$casual_prop 0\" >> $models_path/bolt_no_cheat/params_\$i
    done; done
    $ldak --mega-prs $models_path/bolt_no_cheat/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_no_cheat/params_\$i   
    $ldak --mega-prs $chunked_compheno_models/bolt_no_cheat/sample_\${i} --model bolt --summary $gwas_chunked_combine_pheno/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/bolt_no_cheat/params_\$i   
done
""" > scripts/job.3.test_models.bolt_no_cheat.${he}.${casual_markers}
done; done


# bolt non-zero f1
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
model_name=bolt_no_cheat_non_zero_f1
echo """$(cat base_slurm.sh)
mkdir $models_path/$model_name
mkdir $chunked_compheno_models/$model_name
for ((i=0; i<$num_simulations; i++)); do
    for he_t in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
    for ((as_c=1; as_c<=10; as_c++ )); do
        # make the params file for bolt
        num_lines=\$(wc -l < $gwas_sum_stats_path/sample_\${i}.summaries)
        casual_prop=\$(echo \"scale=7; \$as_c/\$num_lines\" | bc)
        echo \"\$he_t 0\$casual_prop 0.0001\" >> $models_path/$model_name/params_\$i
    done; done
    $ldak --mega-prs $models_path/$model_name/sample_\${i} --model bolt --summary $gwas_sum_stats_path/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/$model_name/params_\$i   
    $ldak --mega-prs $chunked_compheno_models/$model_name/sample_\${i} --model bolt --summary $gwas_chunked_combine_pheno/sample_\${i}.summaries --cors $cors_path/sample_\${i} --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters $models_path/$model_name/params_\$i   
done
""" > scripts/job.3.test_models.$model_name.${he}.${casual_markers}
done; done

for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
rm $models_path/sample_*
rm $models_path/params_*
done; done

# bolt no cheat
# same as above but power is 0
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
model_name=bolt_no_cheat_uniform
echo """$(cat base_slurm.sh)
mkdir $models_path/$model_name
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
""" > scripts/job.3.test_models.$model_name.${he}.${casual_markers}
done; done

# susieR default
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/susieR
mkdir $chunked_compheno_models/susieR

for ((i=0; i<$num_simulations; i++)); do
    Rscript python_scripts/run_sussie.R $gwas_sum_stats_path/sample_\${i}.summaries $snps_cors_path/sample_\${i}.ld $models_path/susieR/sample_\${i} $number_of_indiduals
    Rscript python_scripts/run_sussie.R $gwas_chunked_combine_pheno/sample_\${i}.summaries $snps_cors_path/sample_\${i}.ld $chunked_compheno_models/susieR/sample_\${i} $number_of_indiduals
done

""" > scripts/job.3.test_models.susieR.${he}.${casual_markers}
done;done

## susieR default paper settings
#for he in $HE_VALUES; do
#for casual_markers in $CASUAL_VALUES; do
#source simulate_chunks_loop_const.sh
#model_name=susieR_paper_setting
#
#echo """$(cat base_slurm.sh)
#mkdir $models_path/susieR_paper_setting
#mkdir $chunked_compheno_models/susieR_paper_setting
#
#for ((i=0; i<$num_simulations; i++)); do
#    Rscript python_scripts/run_sussie.R $gwas_sum_stats_path/sample_\${i}.summaries $snps_cors_path/sample_\${i}.ld $models_path/$model_name/sample_\${i} $number_of_indiduals 10 1000 TRUE
#    Rscript python_scripts/run_sussie.R $gwas_chunked_combine_pheno/sample_\${i}.summaries $snps_cors_path/sample_\${i}.ld $chunked_compheno_models/susieR_paper_setting/sample_\${i} $number_of_indiduals 10 1000 TRUE
#done
#
#""" > scripts/job.3.test_models.${model_name}.${he}.${casual_markers}
#done; done

# susieR cheat casusal
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
model_name=susie_cheat_casusal
echo """$(cat base_slurm.sh)
mkdir $models_path/$model_name
mkdir $chunked_compheno_models/$model_name
for ((i=0; i<$num_simulations; i++)); do
    Rscript python_scripts/run_sussie.R $gwas_sum_stats_path/sample_\${i}.summaries $snps_cors_path/sample_\${i}.ld $models_path/$model_name/sample_\${i} $number_of_indiduals $casual_markers 1000 TRUE
    Rscript python_scripts/run_sussie.R $gwas_chunked_combine_pheno/sample_\${i}.summaries $snps_cors_path/sample_\${i}.ld $chunked_compheno_models/$model_name/sample_\${i} $number_of_indiduals $casual_markers 1000 TRUE
done
""" > scripts/job.3.test_models.$model_name.${he}.${casual_markers}
done; done

# fine map default
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
model_name=finemap
file_name=finemap_file
echo """$(cat base_slurm.sh)
mkdir $models_path/$model_name
mkdir $chunked_compheno_models/$model_name
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name_chunk_compheno/$file_name', z_file_base='$experiment_name_chunk_compheno/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$chunked_compheno_models/$model_name')\"
$finemap --sss --in-files $experiment_name_chunk_compheno/$file_name --n-threads 4
$finemap --sss --in-files $experiment_name/$file_name --n-threads 4
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$chunked_compheno_models/$model_name', range($num_simulations))\"
""" > scripts/job.3.test_models.$model_name.${he}.${casual_markers}
done; done


# fine map cheat casual
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
model_name=finemap_cheat_casual
file_name=finemap_cheat_casual_file
echo """$(cat base_slurm.sh)
mkdir $models_path/$model_name
mkdir $chunked_compheno_models/$model_name
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name_chunk_compheno/$file_name', z_file_base='$experiment_name_chunk_compheno/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$chunked_compheno_models/$model_name')\"
$finemap --sss --in-files $experiment_name_chunk_compheno/$file_name --n-threads 4 --n-causal-snps $casual_markers
$finemap --sss --in-files $experiment_name/$file_name --n-threads 4 --n-causal-snps $casual_markers
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$chunked_compheno_models/$model_name', range($num_simulations))\"
""" > scripts/job.3.test_models.$model_name.${he}.${casual_markers}
done; done


#### genome wide models. #####
# fit bolt no cheat.
# same as above but power is 0
model_name=bolt_no_cheat
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
### fill in below
echo """$(cat base_slurm.sh)
mkdir $combine_models
mkdir $combine_models/$model_name
rm $combine_models/${model_name}/params
touch $combine_models/${model_name}/params
for he_t in 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
for casual_prop in 0.0001 0.0005 0.001 0.002 0.004 0.008 0.016 0.032 0.064 0.128 0.256; do
    # make the params file for bolt
    echo \$he_t \$casual_prop 0 >> $combine_models/${model_name}/params
done; done
$ldak --mega-prs $combine_models/$model_name/sample_0 --model bolt --summary $gwas_sum_stats_combined_path/gwas.summaries --cors $combined_samples/cors --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4 --parameters  $combine_models/${model_name}/params
""" > scripts/job.3.test_models.${model_name}_combine.${he}.${casual_markers}
done; done

model_name=bayesr_default
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $combine_models/$model_name
$ldak --mega-prs $combine_models/$model_name/sample_0 --model bayesr --summary $gwas_sum_stats_combined_path/gwas.summaries --cors $combined_samples/cors --power 0.25 --high-LD highld.txt --allow-ambiguous YES --max-threads 4
""" > scripts/job.3.test_models.${model_name}_combine.${he}.${casual_markers}
done; done



#######################
#######################
####   Analyisis   ####
#######################
#######################

for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)

for model in \$(find $models_path -mindepth 1 -maxdepth 1 -type d); do
    python3 -c \"import python_scripts.utils as utils; utils.get_stats_chunk_1_snp(pvals_file='$gwas_sum_stats_path', probs_file='\$model', casual_file='$casual_snp_info', model='\$model',sim_numbers=range($num_simulations))\"
done

for model in \$(find $chunked_compheno_models -mindepth 1 -maxdepth 1 -type d); do
    python3 -c \"import python_scripts.utils as utils; utils.get_stats_chunk_1_snp(pvals_file='$gwas_chunked_combine_pheno', probs_file='\$model', casual_file='$casual_snp_info', model='\$model',sim_numbers=range($num_simulations))\"
done

mkdir $experiment_name/stacked_figures
for ((i=0; i<$num_simulations; i++)); do
python3 -c \"import python_scripts.utils as utils; utils.make_stacked_plots('$gwas_sum_stats_path', '$models_path', '$experiment_name/stacked_figures', '$casual_snp_info', \$i)\"
done

python3 -c \"import python_scripts.utils as utils; utils.get_json_summary('$gwas_sum_stats_path', '$models_path', '$experiment_name/summary.json', '$casual_snp_info', sim_numbers=range($num_simulations))\"
python3 -c \"import python_scripts.utils as utils; utils.get_json_summary('$gwas_sum_stats_path', '$models_path', '$experiment_name/summary_ignore_bad.json', '$casual_snp_info', sim_numbers=range($num_simulations), ignore_bad=True)\"

python3 -c \"import python_scripts.utils as utils; utils.get_json_summary('$gwas_chunked_combine_pheno', '$chunked_compheno_models', '$experiment_name_chunk_compheno/summary.json', '$casual_snp_info', sim_numbers=range($num_simulations))\"
#python3 -c \"import python_scripts.utils as utils; utils.get_json_summary('$gwas_sum_stats_combined_path', '$combine_models', '$experiment_name_combined/summary.json', '$casual_snp_info', sim_numbers=range($num_simulations))\"

mkdir $experiment_name/summary_plots
mkdir $experiment_name_chunk_compheno/summary_plots
python3 -c \"import python_scripts.utils as utils; utils.model_compare_figures('$experiment_name/summary.json', '$experiment_name/summary_plots/summary_plots', statistics=None)\"
python3 -c \"import python_scripts.utils as utils; utils.model_compare_figures('$experiment_name/summary_ignore_bad.json', '$experiment_name/summary_plots/summary_plots_summary_ignore_bad', statistics=None)\"

python3 -c \"import python_scripts.utils as utils; utils.model_compare_figures('$experiment_name_chunk_compheno/summary.json', '$experiment_name_chunk_compheno/summary_plots/summary_plots', statistics=None)\"
""" > scripts/job.4.sum_stats.${he}.${casual_markers}
done; done

source simulate_chunks_const.sh

# need to make this a job and to add genome-wide lines
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
echo """$(cat base_slurm.sh)
source simulate_chunks_loop_const.sh
python3 -c \"import python_scripts.utils as utils; utils.plot_power_fdr('$models_path', '$casual_snp_info', '$chunked_simulations/power_fdr_he_${he}_casual_${casual_markers}.png', range(${num_simulations}), title='he_${he}_casual_${casual_markers}')\"
python3 -c \"import python_scripts.utils as utils; utils.plot_power_fdr('$chunked_compheno_models', '$casual_snp_info', '$chunked_compheno_simulations/power_fdr_he_${he}_casual_${casual_markers}.png', range(${num_simulations}), title='he_${he}_casual_${casual_markers}', genome_wide_models='$combine_models')\"
""" > scripts/job.4.fdr.${he}.${casual_markers}
done;done







squeue -u $USER -h | awk '{print $1}' | xargs scancel