#!/bin/bash

if [ $# -ne 0 ];then
        echo "Usage: sh $0"
        exit 1
fi

echo [MSG] Start submitting cellsnp-lite jobs...
while read line
do
	qsub cellsnp-lite.pbs -v capture=$line,workDIR=`pwd`
done < "capture.list" #"capture.list"

echo [MSG] All jobs submitted! 
