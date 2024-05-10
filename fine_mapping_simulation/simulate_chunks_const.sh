#### Constants #####
num_simulations=200
random_seed=42
number_of_indiduals=50000 # max is 392214
sample_size=2000000
output_path=small_chunks_current
ubb_binary_path=/faststorage/project/dsmwpred/data/ukbb/geno
subset_individuals=${output_path}/rand${number_of_indiduals}_included.fam
genome_sampled_markers_path=$output_path/genome_sample
genome_sample_binary_path=$output_path/genome_binary
snps_cors_path=$output_path/snps_cors
cors_path=$output_path/cors

finemap='/home/patrickgibbs/dsmwpred/patrick/finemap_v1.4.2_x86_64/finemap_v1.4.2_x86_64'
plink='/home/patrickgibbs/dsmwpred/patrick/plink/plink'
ldak='~/dsmwpred/patrick/ldak'

mem=32G
cores=4
time=01:00:00
