
# fine map searching effect sizes
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
model_name=finemap_search_effect_size
file_name=finemap_search_effect_size_file
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
#mkdir $experiment_name/finemap_z_files
mkdir $models_path/$model_name
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
$finemap --sss --in-files $experiment_name/$file_name --n-threads 4 --prior-std 0.0001,0.00025,0.0005,0.001,0.0025,0.005,0.01,0.05,0.1,0.5,1
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
""" > scripts/job.test_models.$model_name.${he}.${casual_markers}
done; done



# fine map step wise conditioning
for he in $HE_VALUES; do
for casual_markers in $CASUAL_VALUES; do
model_name=finemap_step_wise
file_name=finemap_step_wise_file
source simulate_chunks_loop_const.sh
echo """$(cat base_slurm.sh)
mkdir $models_path/$model_name
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name/$file_name', z_file_base='$experiment_name/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$models_path/$model_name')\"
python3 -c \"import python_scripts.utils as utils; utils.make_finemap_file(f_name='$experiment_name_chunk_compheno/$file_name', z_file_base='$experiment_name_chunk_compheno/finemap_z_files', cors_base='$snps_cors_path', num_samples=$number_of_indiduals,  samples=range($num_simulations), model_path='$chunked_compheno_models/$model_name')\"
$finemap --cond --in-files $experiment_name/$file_name
python3 -c \"import python_scripts.utils as utils; utils.make_fine_map_probs('$models_path/$model_name', range($num_simulations))\"
""" > scripts/job.3.test_models.$model_name.${he}.${casual_markers}
done; done