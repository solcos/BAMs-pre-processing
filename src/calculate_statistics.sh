#!/bin/bash

# Calculate final statistics for pre-process BAMs

samtools flagstat $1.bam > statistics/$1.flagstat.statistics
samtools stats $1.bam | grep ^SN | cut -f 2- > statistics/$1.stats

