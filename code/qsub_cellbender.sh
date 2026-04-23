#!/bin/bash

if [ $# -ne 0 ];then
        echo "Usage: sh $0"
        exit 1
fi

echo [MSG] Start submitting cellbender jobs...
while read line
do
	qsub cellbender.pbs -v capture=`echo $line | cut -f1 -d','`,workDIR=`pwd`
done < "capture3.list" #"capture.list"

echo [MSG] All jobs submitted! 
