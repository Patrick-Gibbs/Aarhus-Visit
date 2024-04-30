"""Constructor program to execude `fine_mapping_chunks_base.sh with different parameters"""

import pandas as pd

cores = 4
mem = 4*cores

SLURM=f"""#!/bin/bash
#SBATCH -c {cores}
#SBATCH --mem {mem}g\n\n
"""

def generate_simuluation(
    script_name,
    base_path,
    num_simulations=1,
    he=0.1,
    random_seed=42,
    casual_markers=3,
    number_of_indiduals=50000, # max is 392214
    sample_size=2000000,
    model='bolt',
    default_model=True,
    params_path='params'
    ):

    f = open(f"scripts/{script_name}", 'w')

    f.write(SLURM)

    experiment_name=f'small_chunks/chucks_1_marker_her_{he}_ind_{number_of_indiduals}_casual_{casual_markers}'
    f.write(f"""            
#### Simulation Variables ####
num_simulations={num_simulations}
he={he}
random_seed={random_seed}
casual_markers={casual_markers}
number_of_indiduals={number_of_indiduals} # max is 392214
sample_size={sample_size}
model={model}
default_model={int(bool(default_model))}
params_path={params_path}
experiment_name={experiment_name}
params_path={params_path}
############################\n\n\n""")

  

    base_code_import = open(base_path, 'r').read()
    start_base_code = '# Code to Run Simulation Below #'
    base_code_clean = base_code_import[base_code_import.index(start_base_code):]
    f.write(base_code_clean)


# test
generate_simuluation('build.test', 'fine_mapping_chunks_base.sh', params_path='params')