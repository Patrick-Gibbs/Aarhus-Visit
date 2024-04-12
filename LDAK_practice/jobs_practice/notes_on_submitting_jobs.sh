# to submit a job
sbatch my_first_job.sh

# to check on my jobs
squeue -u patrickgibbs

# also can use to list all jobs in the que
squeue

# jobs can be cancelled using the scancel command:
scancel 17129500

# 
srun --account dsmwpred -c 4 --mem 32g --time 01:00:00 --pty bash
