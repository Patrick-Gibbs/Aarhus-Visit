declare -a jid
declare -a jid_prev
n=5
jid=()
jid_prev=()
for i in $(seq 0 $n); do
    echo $i
    for job in $(ls scripts/job.$i.*); do
        if [ $i -eq 0 ] || [ ${#jid[@]} -eq 0 ]; then
            jid_prev+=("$(sbatch -A dsmwpred $job --job-name=$job --output=scripts/std-out | awk '{print $4}')")
        else
            jid_prev+=("$(sbatch -A dsmwpred --dependency=afterok:$(IFS=","; echo "${jid[*]}"; unset IFS) $job --job-name=$job --output=scripts/std-out | awk '{print $4}')")
        fi
        mv $job scripts/past-jobs/
    done
    jid=("${jid_prev[@]}")
done

#start any remaining jobs
for job in $(ls scripts/job.*); do
   sbatch -A dsmwpred $job --job-name=$job --output=scripts/std-out
   mv $job scripts/past-jobs/
done


# test case
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.0.test
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.1.test
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.2.test
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.2.test
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.2.test
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.3.test
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.4.test
echo "#"'!'"/bin/bash
#SBATCH --mem 1G
#SBATCH -c 1
#SBATCH -t 00:05:00
sleep 10s
echo hi" > scripts/job.5.test