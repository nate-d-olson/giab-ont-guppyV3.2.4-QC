#!/usr/bin/bash
## Generate targeted regions bed for evaluating phased bams
# Analysis performed on 11/18/2019 by Nate Olson based on conversation with Justin Zook.


## Input files
### hs37d5
hs37d5_cds_bed=/Volumes/giab/analysis/benchmarking-tools/resources/stratification-bed-files/FunctionalRegions/refseq_union_cds.sort.bed.gz
hs37d5_v3_bed=/Volumes/giab/analysis/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed 

### GRCh38
GRCh38_cds_bed=/Volumes/giab/analysis/GRCh38_stratification_bed_files/FunctionalRegions/GRCh38_refseq_cds_merged.bed.gz
GRCh38_v3_bed=/Volumes/giab/analysis/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh38/HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_noinconsistent.bed 

## Getting bedtools version info

bedtools --version

## Generating Interset
gunzip -ck ${hs37d5_cds_bed} \
    | bedtools intersect -a stdin -b ${hs37d5_v3_bed}\
    | bedtools sort \
    | bedtools merge -d 50 \
    > resources/phase_eval_targets_HG002_GRCh37.bed 

gunzip -ck ${GRCh38_cds_bed} \
    | bedtools intersect -a stdin -b ${GRCh38_v3_bed}\
    | bedtools sort \
    | bedtools merge -d 50 \
    > resources/phase_eval_targets_HG002_GRCh38.bed 