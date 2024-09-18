PLINKPATH=~/lstm_data/cease/variants_bycohort/combined_cohorts/plink

for S in 1 2 3 4 5
do
        SEED=$RANDOM
        mkdir ${SEED}_run
        cd ${SEED}_run
        for K in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
                do
                        ~/software/admixture_linux-1.3.0/admixture --cv --seed $SEED -j10 $PLINKPATH/CM023248.thin0.1.allinds_50.10.0.1.pruned.bed $K | tee log${K}.${SEED}.${S}.out
                done
        cd ..
done
