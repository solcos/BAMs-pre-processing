#!/bin/bash

# Calculate depth of coverage

samtools depth  $1.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > statistics/$1.depth_of_coverage
