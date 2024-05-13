source simulate_chunks_const.sh
experiment_name=${chunked_simulations}/chucks_1_marker_her_${he}_ind${number_of_indiduals}_casual_${casual_markers}
genome_sampled_marks_dir=$experiment_name/genome_regions_map
genome_sample_binary_dir=$experiment_name/genome_sample
phenotype_path=$experiment_name/phenotypes
gwas_sum_stats_path=$experiment_name/fine_map_gwas
casual_snp_info=$experiment_name/casual_snp_info
models_path=$experiment_name/models



### Combine Geno and Pheno ###
experiment_name_combined=$combined_simulations/combined_her_${he}_ind${number_of_indiduals}_casual_${casual_markers}
combined_pheno=$combined_samples/combined_pheno_her_${he}_ind${number_of_indiduals}_casual_${casual_markers}
gwas_sum_stats_combined_path=$experiment_name_combined/fine_map_gwas
combine_models=$experiment_name_combined/models


## combine pheno chunked simulations ##
experiment_name_chunk_compheno=$chunked_compheno_simulations/combined_her_${he}_ind${number_of_indiduals}_casual_${casual_markers}
gwas_chunked_combine_pheno=$experiment_name_chunk_compheno/fine_map_gwas
chunked_compheno_models=$experiment_name_chunk_compheno/models
