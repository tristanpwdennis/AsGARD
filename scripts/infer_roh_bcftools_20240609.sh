VCFPATH=/mnt/temp_scratch/vcf/

for chrom in CM023250 CM023248 CM023249
	do
		bcftools roh $VCFPATH/$chrom'.noindels.filtered.ann.vcf.gz' --threads 20 -G30 -M 1.36e-9 --output-type r -o $chrom'.bcftools.gt.roh.txt'
		bcftools roh $VCFPATH/$chrom'.noindels.filtered.ann.vcf.gz' --threads 20 -M 1.36e-9 --output-type r -o $chrom'.bcftools.pl.roh.txt'
	done
