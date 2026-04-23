#!/bin/bash

if [ $# -ne 1 ];then
        echo "Usage: sh $0 <capture_name>"
        exit 1
fi

capture=$1

cellrangerDIR="/group/canc2/anson/working/cellranger/${capture}"
REGIONVCF="/group/canc2/anson/reference/region/cellsnplite/genome1K.phase3.SNP_AF5e4.chr1toX.hg38.vcf.gz"
outDIR=`pwd`/data/cellsnp-lite

mkdir -p ${outDIR}

echo "Start calling SNPs from ${capture} ..."

/group/canc2/anson/miniconda3/bin/cellsnp-lite \
--samFile ${cellrangerDIR}/outs/per_sample_outs/${capture}/count/sample_alignments.bam \
--outDir ${outDIR}/${capture} \
--regionsVCF ${REGIONVCF} \
--barcodeFile `pwd`/data/cellbender/${capture}/${capture}.cellbender_cell_barcodes.csv \
--nproc 16 \
--minMAF 0.1 \
--minCOUNT 20 \
--gzip

echo "SNP calling from ${capture} done!"
