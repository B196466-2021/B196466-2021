#!/bin/bash

source ~/.bashrc
# get SRR ID from command line argument
SRR=$1

cd /home/spdx/scRNA-seq/scRNA/Cellranger_result
# step 1: prefetch
/home/spdx/scRNA-seq/sratoolkit.2.10.7-ubuntu64/bin/prefetch.2.10.7 $SRR --max-size 20G

cd $SRR
# step 2: fastq-dump
/home/spdx/scRNA-seq/sratoolkit.2.10.7-ubuntu64/bin/fastq-dump.2.10.7 --gzip --split-files ${SRR}.sra

# step 3: rename files
mv ${SRR}_1.fastq.gz ${SRR}_S1_L001_R1_001.fastq.gz
mv ${SRR}_2.fastq.gz ${SRR}_S1_L001_R2_001.fastq.gz

# step 4: run cellranger count
/home/spdx/scRNA-seq/cellranger-7.1.0/bin/cellranger count --id=${SRR}_out \
                 --transcriptome=/home/spdx/scRNA-seq/refdata-gex-GRCh38-2020-A \
                 --fastqs=/home/spdx/scRNA-seq/scRNA/Cellranger_result/$SRR \
                 --sample=$SRR