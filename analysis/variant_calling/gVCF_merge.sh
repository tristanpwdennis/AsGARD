COMBINED_VCF_OUTDIR=/my_data/Clara_Parabricks/gVCF_merge/variant_calling_workspace/cease/variants_bycohort/combined_cohorts/
REF=/my_data/cease/genomes/anstep/VectorBase-61_AstephensiUCISS2018_Genome.fasta

INTERVAL=$1

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx10g -Xms10g' GenomicsDBImport --genomicsdb-workspace-path /my_data/AnstepGenomicsDB_Dec24/$INTERVAL.db --sample-name-map map_20241228.txt --reader-threads 5 --genomicsdb-shared-posixfs-optimizations true --intervals $INTERVAL --consolidate TRUE --batch-size 50

gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx10g -Xms10g' GenotypeGVCFs -R $REF -V gendb://AnstepGenomicsDB_Dec24/$INTERVAL.db -O file:////my_data/cease_variants/combined_cohorts/vcf_intervals/$INTERVAL.vcf.gz
