#!/usr/bin/env bash

if [ -z "$1" ] || [ -z "$2" ] || [ -z "$3" ]; then
    echo "Usage: $0 <input_directory> <output_directory> <number_of_cores>"
    exit 1
fi

RUN_PATH=$1
OUT_PATH=$2
N_CORES=$3

mkdir -p $OUT_PATH
cd $RUN_PATH
for file in $(ls $RUN_PATH)
do
    SAMPLE=`basename $file`
    fastqc -q -t $N_CORES -o ${OUT_PATH} ${SAMPLE}
    echo "Processed $SAMPLE"
done

