#!/bin/bash

if [ $# -ne 1 ];then
        echo "Usage: sh $0 <capture_name>"
        exit 1
fi

capture=$1

input="/group/canc2/anson/working/cellranger/${capture}/outs/raw_feature_bc_matrix.h5"
outDIR="/group/canc2/anson/working/cf-pbmc-bal/data/cellbender/${capture}"

mkdir -p ${outDIR} && cd ${outDIR}

echo "Start cellbender for ${capture} ..."

source /group/canc2/anson/miniconda3/etc/profile.d/conda.sh
conda activate cellbender
python --version
cellbender remove-background --cuda --input $input --output ${outDIR}/${capture}.cellbender.h5

