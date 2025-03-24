#!/bin/bash

OUTDIR=/lstm_data/cease/analysis/phasing_20240702
IBDDIR=/lstm_data/cease/analysis/ibd_20240702
VCFPATH=/lstm_data/cease/variants_bycohort/combined_cohorts/vcf/
BEAGLEBIN=/lstm_data/cease/analysis/scripts/beagle.22Jul22.46e.jar
IBDJAR=/home/software/hap-ibd/hap-ibd.jar

for CHROM in CM023248 CM023249 CM023250
do
	java -Xms500g -Xmx500g -jar $BEAGLEBIN nthreads=20 gt=$VCFPATH/$CHROM.noindels.filtered.ann.vcf.gz out=$OUTDIR/$CHROM.phased.vcf.gz window=10 overlap=1
done
