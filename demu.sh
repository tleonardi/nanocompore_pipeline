#!/bin/bash

DATA_DIR="/nfs/leia/research/enright/nanopore/datasets/7SK_gurdon_runs/20190220_DKC1_MOLM_CTRL-sh94_7SK/"
OUTDIR="/hps/nobackup/enright/tom/new_nanocompore/7SK_DKC1/demultiplexed"
mkdir -p $OUTDIR/logs/

bgadd -L 1 /demuxing

for fastfive in $DATA_DIR/*fast5; do
	samp_name=$(basename $fastfive .fast5)
	bsub -g /demuxing -q standard -eo $OUTDIR/logs/${samp_name}.err -oo $OUTDIR/logs/${samp_name}.out -M 8000 -n 2 -R 'rusage[mem=8000]' ./poreplex_demux.sh $fastfive $OUTDIR
done



