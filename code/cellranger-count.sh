#!/bin/bash

if [ $# -ne 2 ];then
        echo "Usage: sh $0 <capture> <workDIR>"
        exit 1
fi

capture=$1
workDIR=$2

module load cellranger/7.1.0
cd $workDIR
cellranger count --sample=${capture}-GEX --id=${capture} --fastqs /ngs/canc2-ngs/Anson_Wong_scRNA/240423_A00692_0408_AH3YGKDSXC_DI --transcriptome /group/canc2/anson/reference/genome/cellranger/refdata-gex-GRCh38-2020-A --localcores=12 --localmem=48
