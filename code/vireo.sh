#!/bin/bash

if [ $# -ne 2 ];then
        echo "Usage: sh $0 <capture_name> <n_donor>"
        exit 1
fi

capture=$1
n_donor=$2

cellSNPDIR=`pwd`/data/cellsnp-lite/${capture}
outDIR=`pwd`/data/vireo/${capture}

mkdir -p ${outDIR}

echo "Start assigning barcode to donors for ${capture} ..."

source /group/canc2/anson/miniconda3/etc/profile.d/conda.sh
conda activate cellbender

vireo \
--cellData=${cellSNPDIR} \
--outDir=${outDIR} \
--nDonor=${n_donor} \
--nInit=200 \
--randSeed=1990 \
--nproc=12
