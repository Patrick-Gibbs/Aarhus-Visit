#### Constants #####
HE_VALUES="0.02"
CASUAL_VALUES="1 3"

num_simulations=10
random_seed=42
number_of_indiduals=50000 # max is 392214
sample_size=2000000
output_path=small_chunks_current
ukbb_binary_path=/faststorage/project/dsmwpred/data/ukbb/geno
ukbb_cors=$output_path/ukbb_cors

subset_individuals=${output_path}/rand${number_of_indiduals}_included.fam
genome_sampled_markers_path=$output_path/genome_sample
genome_sample_binary_path=$output_path/genome_binary
snps_cors_path=$output_path/snps_cors
cors_path=$output_path/cors
chunked_simulations=${output_path}/chunked_sims

finemap='/home/patrickgibbs/dsmwpred/patrick/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64'
plink='/home/patrickgibbs/dsmwpred/patrick/plink/plink'
ldak='~/dsmwpred/patrick/ldak'

### Combine ###
combined_samples=$output_path/genome_sample_combine
combined_samples_regions=$combined_samples/_combined_samples
combined_simulations=${output_path}/combined_sims


### chunked combine pheno
chunked_compheno_simulations=$output_path/chunk_compheno_sims





### Slurm ###
mem=32G
cores=4
time=01:00:00
