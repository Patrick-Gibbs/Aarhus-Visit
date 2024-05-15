#### Constants #####

HE_COMBINE="0.1 0.2 0.5 0.75"
HE_CHUNKED="0.02 0.05 0.1"
CASUAL_VALUES="1 3"

HE_VALUES=""

for value in $HE_COMBINE $HE_CHUNKED; do
    if [[ ! $value =~ (^|[[:space:]])$HE_VALUES($|[[:space:]]) ]]; then
        HE_VALUES="$HE_VALUES $value"
    fi
done

HE_VALUES=${HE_VALUES:1}  # Remove the leading space

num_simulations=20
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
