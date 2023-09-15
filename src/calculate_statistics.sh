#!/bin/bash

# Calculate depth of coverage

samtools depth  $1.bam  |  awk '{sum+=$3} END { print "Average = ",sum/NR}' > statistics/$1.depth_of_coverage

slc00023@sllogin2:/slgpfs/projects/slc00/sv_data> cat src/calculate_statistics.sh 
#!/bin/bash

# Calculate final statistics for pre-process BAMs

samtools flagstat $1.bam > statistics/$1.flagstat.statistics
samtools stats $1.bam | grep ^SN | cut -f 2- > statistics/$1.stats
