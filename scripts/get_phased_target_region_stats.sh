#!/usr/bin/bash

## Get bedcov stats for hs37d5 HG002 for phase eval targeted regions ###########
BAMDIR=/Volumes/giab/data/alignment/AshkenazimTrio/Ultralong_OxfordNanopore/ULbams_guppy3.2.4
TARGETSBED=resources/phase_eval_targets_subset_GRCh37.bed.gz

## Subset bam to target regions
samtools --version

samtools view \
  -L ${TARGETSBED}  \
  ${BAMDIR}/phased_hs37d5.bam -O bam \
  > tmp/phase_targets_hs37d5.bam

samtools index tmp/phase_targets_hs37d5.bam

## Split bam by haplotype using `bedtools filter`
bamtools filter -in tmp/phase_targets_hs37d5.bam \
    -out tmp/phase_targets_subset_GRCh37_HP1.bam \
    -tag HP:1
samtools index tmp/phase_targets_subset_GRCh37_HP1.bam

bamtools filter -in tmp/phase_targets_hs37d5.bam \
    -out tmp/phase_targets_subset_GRCh37_HP2.bam \
    -tag HP:2
samtools tmp/phase_targets_subset_GRCh37_HP2.bam

## Compute number of reads mapped to each region total (HP1, HP2, and unphased), HP1, and HP2. 
samtools bedcov \
    ${TARGETSBED}  \
    tmp/phase_targets_hs37d5.bam \
    > data/guppy_V3.2.4/phase_eval/phase_targets_subset_GRCh37_bedcov.tsv
    
samtools bedcov \
    ${TARGETSBED} \
    tmp/phase_targets_subset_GRCh37_HP1.bam \
    > data/guppy_V3.2.4/phase_eval/phase_targets_subset_GRCh37_HP1_bedcov.tsv
    
samtools bedcov \
    ${TARGETSBED} \
    tmp/phase_targets_subset_GRCh37_HP2.bam \
    > data/guppy_V3.2.4/phase_eval/phase_targets_subset_GRCh37_HP2_bedcov.tsv
    
    
## GRCh38 - code moved from Rmarkdown - can convert to script if deemed necessary
# __Calculating Read Counts by haplotype and region__ 
# Subsetting bam to target regions
# ```
# samtools view \
#   -L ../resources/phase_eval_targets_HG002_GRCh38.bed  \
#   /Volumes/giab/data/alignment/AshkenazimTrio/Ultralong_OxfordNanopore/ULbams_guppy3.2.4/phased_GRCh38.bam -O bam \
#   > phase_targets_GRCh38.bam
# ```
# 
# Using `bedtools filter` to seperate haplotypes
# `bamtools filter -in phase_subset_target_regions.bam -out phase_subset_target_regions_HP1.bam -tag HP:1`  
# `bamtools filter -in phase_subset_target_regions.bam -out phase_subset_target_regions_HP2.bam -tag HP:2`  
# `samtools index phase_subset_target_regions_HP1.bam` 
# `samtools index phase_subset_target_regions_HP2.bam`
# 
# Compute number of reads mapped to each region total (HP1, HP2, and unphased), HP1, and HP2. 
# `samtools bedcov ../resources/phase_eval_targets_subset_GRCh38.bed.gz phase_subset_target_regions.bam > phase_subset_target_regions_bedcov.tsv`
# `samtools bedcov ../resources/phase_eval_targets_subset_GRCh38.bed.gz phase_subset_target_regions_HP1.bam > phase_subset_target_regions_HP1_bedcov.tsv`
# `samtools bedcov ../resources/phase_eval_targets_subset_GRCh38.bed.gz phase_subset_target_regions_HP2.bam > phase_subset_target_regions_HP2_bedcov.tsv`