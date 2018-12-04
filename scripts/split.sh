#/usr/bin/env bash

OUT=split
mkdir -p keys ${OUT}
awk -vOFS="\t" '{gsub(/\-I./,"",$2);gsub(/Kontrola/,"Kontrola_"$1"_"$2,$6);print > "keys/"$1"_"$2}' meta/reformatted.tsv

time for run in $(seq 1 6); do
  for index in 1 2 3; do
    R1=data/Hordeum_run${run}-*/FASTQ_Generation_*/${index}_*/*R1_001.fastq.gz;
    R2=${R1/R1/R2};
    # OUT=out_split_run${run}_index${index}
    KEY=keys/key_${run}_${index}
    cat keys/Run_index keys/${run}_${index} | cut -f3- > ${KEY}
    echo -e "Processing run ${run}, index ${index}"
    ls ${R1} ${R2}
    paste <(pigz -dcp2 ${R1} | paste - - - - ) <(pigz -dcp2 ${R2} | paste - - - - ) \
      | java -jar scripts/merutensils.jar split \
          --key-file ${KEY} \
          --min-length-single-read 50 \
          --min-length-pair-each 50 \
          --min-length-pair-sum 100 \
          --print-user-settings \
          --out-dir ${OUT} \
          --blank-samples-name Kontrola \
          2>&1 | tee ${OUT}/summary_${run}_${index}.txt
  done
done