#!/usr/bin/env bash

biokanga align \
  --sfx ${idxmeta.target}.sfx \
  --mode 0 \
  --format 5 \
  --pemode 2 \
  --in ${r1} \
  --pair ${r2}  \
  --out ${sample}.bam \
  --snpreadsmin 6 \
  --snpfile ${sample}.vcf \
  --threads ${task.cpus} \
  --substitutions 3 \
  --pairminlen 50
  #--minchimeric 50
  # --snpfile ${sample}
  # --snpcentroid ${sample}