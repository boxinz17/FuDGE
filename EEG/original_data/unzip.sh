#!/bin/bash

for f in *.tar.gz; do tar xf "$f"; done

for f in *.tar.gz
do
    tar -xvzf "$f" -C /eeg_full/temp
done