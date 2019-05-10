#!/bin/bash
set -e -o pipefail
# Runs poreplex to demultiplex multi-fast5 files.
# Usage: poreplex_demux.sh /path/to/file.fast5 /out/path 
# The same out path can be used for multiple input files

THREADS=2
INFILE=$1 
DEST_DIR=$2 
INFILE_BASENAME=$(basename $INFILE .fast5)

OUT=$(mktemp -d -p /dev/shm)
mkdir "${OUT}/split"
mkdir "${OUT}/demux"
function cleanup {
	rm -rf "$OUT"
}
trap cleanup EXIT

multi_to_single_fast5 -t $THREADS -i "$INFILE" -s "${OUT}/split"
poreplex -i "${OUT}/split" -o "${OUT}/demux" --trim-adapter --barcoding --fast5 --basecall --parallel $THREADS
   
for BC in BC1 BC2 BC3 BC4 undetermined; do
	mkdir -p ${DEST_DIR}/${BC}/${INFILE_BASENAME}
	[[ -d "${OUT}/demux/fast5/pass/$BC" ]] && find  "${OUT}/demux/fast5/pass/$BC" -name "*.fast5" | xargs -I{} cp --no-clobber {} "${DEST_DIR}/${BC}/${INFILE_BASENAME}"
done
echo "$INFILE_BASENAME" >> ${DEST_DIR}/demuxed.txt
 
