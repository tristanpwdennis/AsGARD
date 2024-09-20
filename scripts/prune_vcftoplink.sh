VCFPATH=~/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/
PLINKPATH=~/lstm_data/cease/variants_bycohort/combined_cohorts/plink


~/software/plink --vcf $VCFPATH/CM023248.noindels.filtered.ann.vcf.gz --allow-extra-chr --chr CM023248 --indep-pairwise 50 10 0.1 --thin 0.1 --out $PLINKPATH/CM023248.thin0.1.allinds_50.10.0.1 --threads 20 --make-bed --double-id --set-missing-var-ids @:#

~/software/plink --vcf $VCFPATH/CM023248.noindels.filtered.ann.vcf.gz --allow-extra-chr --chr CM02324 --set-missing-var-ids @:# --extract $PLINKPATH/CM023248.thin0.1.allinds_50.10.0.1.prune.in --make-bed --threads 20 --double-id --out $PLINKPATH/CM023248.thin0.1.allinds_50.10.0.1.pruned

sed -i s/CM023248/2/g $PLINKPATH/$PLINKPATH/CM023248.thin0.1.allinds_50.10.0.1.pruned.bim
