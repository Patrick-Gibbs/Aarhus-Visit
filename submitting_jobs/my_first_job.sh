#!/bin/bash
#SBATCH --account my_project
#SBATCH -c 4
#SBATCH --mem 32g


echo hello world
echo some data output > an_output_file.txt

