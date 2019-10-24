#!/bin/bash

samtools view -o to_count.txt input..bam
wc -l to_count.txt > the_count.txt
echo "SCRIPT COMPLETED"
