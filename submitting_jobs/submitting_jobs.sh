#!/bin/bash
#SBATCH --account my_project
#SBATCH -c 1
#SBATCH --mem 1g


echo hello world
echo some data output > an_output_file.txt

