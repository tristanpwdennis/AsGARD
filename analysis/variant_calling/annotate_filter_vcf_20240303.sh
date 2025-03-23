#annotate and filter vcf
SNPEFF=~/lstm_data/cease/variants_bycohort/annotation_db/snpEff/snpEff.jar
VCFPATH=~/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/

#java -jar $SNPEFF anstep_uci_2018 $VCFPATH/CM023248.vcf.gz  | bcftools view --threads 10 -Oz -S ^samples_ex.txt --max-alleles 2 --exclude-types indels  -i 'F_MISSING < 0.1 && QD > 5 && MQ > 40 && FS < 20 && SOR < 3 && INFO/DP > 6000 && INFO/DP < 13000' > $VCFPATH/CM023248.noindels.filtered.ann.vcf.gz && tabix CM023248.noindels.filtered.ann.vcf.gz &

#java -jar $SNPEFF anstep_uci_2018 $VCFPATH/CM023249.vcf.gz  |  bcftools view  --threads 10 -Oz -S ^samples_ex.txt --max-alleles 2 --exclude-types indels  -i 'F_MISSING < 0.1 && QD > 5 && MQ > 40 && FS < 20 && SOR < 3 && INFO/DP > 6000 && INFO/DP < 13000' > $VCFPATH/CM023249.noindels.filtered.ann.vcf.gz && tabix CM023249.noindels.filtered.ann.vcf.gz &

#java -jar $SNPEFF anstep_uci_2018 $VCFPATH/CM023250.vcf.gz  | bcftools view --threads 10 -Oz -S ^samples_ex.txt --max-alleles 2 --exclude-types indels  -i 'F_MISSING < 0.1 && QD > 5 && MQ > 40 && FS < 20 && SOR < 3 && INFO/DP > 6000 && INFO/DP < 13000' > $VCFPATH/CM023250.noindels.filtered.ann.vcf.gz &&  tabix CM023250.noindels.filtered.ann.vcf.gz &

bcftools view -Oz --threads 5 -R CM023248.accessiblepositions.txt $VCFPATH/CM023248.noindels.filtered.ann.vcf.gz > $VCFPATH/CM023248.noindels.filtered.ann.acc.vcf.gz && tabix $VCFPATH/CM023248.noindels.filtered.ann.acc.vcf.gz &
bcftools view -Oz --threads 5 -R CM023249.accessiblepositions.txt $VCFPATH/CM023249.noindels.filtered.ann.vcf.gz > $VCFPATH/CM023249.noindels.filtered.ann.acc.vcf.gz && tabix $VCFPATH/CM023249.noindels.filtered.ann.acc.vcf.gz &
bcftools view -Oz --threads 5 -R CM023250.accessiblepositions.txt $VCFPATH/CM023250.noindels.filtered.ann.vcf.gz > $VCFPATH/CM023250.noindels.filtered.ann.acc.vcf.gz && tabix $VCFPATH/CM023250.noindels.filtered.ann.acc.vcf.gz &

wait

python3 sgkit_bychr_vcftozarr.py

python3 allel_ann_vcf_to_zarr.20240108.py
