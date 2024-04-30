# continueing on from day1 scripts
# here I am using jackknife to compare ridge and bolt prs for height

cd prs_practice
# we only want see how good the model is at the held out data
awk '{print $1 ','  $2}' /faststorage/project/dsmwpred/data/ukbb/height.test > height_testing_accessions

# held out predictions
ldak --calc-scores scores_ridge_height --scorefile prs_ridge_height.effects --bfile /faststorage/project/dsmwpred/data/ukbb/geno --power 0 --keep height_testing_accessions
ldak --calc-scores scores_bolt_height --scorefile prs_bolt_height.effects --bfile /faststorage/project/dsmwpred/data/ukbb/geno --power 0 --keep height_testing_accessions

# exctracting predictions
cut -f5 scores_ridge_height.profile | tail -n +2 > ridge_pred_height_test
cut -f5 scores_bolt_height.profile | tail -n +2 > bolt_pred_height_test
cut -f3 /faststorage/project/dsmwpred/data/ukbb/height.test -d ' ' > obs_height_test



# csv with predictions and observations
paste ridge_pred_height_test obs_height_test > ridge_height_pred_obs
paste bolt_pred_height_test obs_height_test > bolt_height_pred_obs

# measures how accurate the model was 
ldak --jackknife jackknife_height_ridge --data-pairs ridge_height_pred_obs --num-blocks 200
ldak --jackknife jackknife_height_bolt --data-pairs bolt_height_pred_obs --num-blocks 200

## Ridge
## Measure Estimate SD
## Correlation 0.530686 0.004642
## Squared_correlation 0.281627 0.004927
## Mean_squared_error 0.727765 0.006006
## Mean_absolute_error 0.677080 0.003111
## 
## Bolt
## Measure Estimate SD
## Correlation 0.573852 0.004680
## Squared_correlation 0.329306 0.005372
## Mean_squared_error 0.673201 0.006051
## Mean_absolute_error 0.648416 0.003020