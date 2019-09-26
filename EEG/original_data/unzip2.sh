#!/bin/bash
cd /eeg_full/temp
for D in `ls`
do
    cd /eeg_full/temp/${D}
    for f in *.gz
    do
	gunzip -v "$f"
    done
done