#!/bin/sh
#gatk pipeline
####This WORKS for loopimng over files - it was temperamental to setup so duplicate and do not modify
INPUT_DIR=/mnt/hp_scratch/Users/tristand/Clara_Parabricks/Germline/data
OUTPUT_DIR=/mnt/hp_scratch/Users/tristand/Clara_Parabricks/Germline/ancestral_rec_21.6.24
REFERENCE_FILE=VectorBase-61_AstephensiUCISS2018_Genome.fasta

for f in `ls $INPUT_DIR/*.fq.gz | sed 's/_\([12]\)\.fq\.gz$//' | uniq | xargs -n1 basename
` #finish this expression - finished by SW
	do
 	docker run --rm --gpus all --volume $INPUT_DIR:/workdir --volume $OUTPUT_DIR:/outputdir \
    -w /workdir \
    nvcr.io/nvidia/clara/clara-parabricks:4.3.0-1 \
    pbrun germline \
    --ref /workdir/$REFERENCE_FILE \
    --in-fq /workdir/"$f"_1.fq.gz /workdir/"$f"_2.fq.gz \
    --out-bam /outputdir/$f.bam \
    --out-variants /outputdir/$f.g.vcf.gz \
    --gvcf \
    --logfile /outputdir/$f.log
    mv  $INPUT_DIR/"$f"_1.fq.gz $INPUT_DIR/done/"$f"_1.fq.gz
    mv  $INPUT_DIR/"$f"_2.fq.gz $INPUT_DIR/done/"$f"_2.fq.gz
    #--in-fq /workdir/$f_1.fq.gz /workdir/$f_2.fq.gz \
    #--out-variants /outputdir/$f.gvcf.gz
	done
