#!/bin/bash


# Author : cbortoluzzi@ethz.ch


#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=40000
#SBATCH --job-name=bamstats
#SBATCH --output=output_%J
#SBATCH --error=error_%J


if [ $# -ne 1 ]
then
	echo -e "\nusage: `basename $0` <bam>\n"
	echo -e "DESCRIPTION:Use Calculate general statistics of aligned sequences in BAM format using bamtools and qualimap\n\n"
	echo -e "INPUT:		<bam>	Aligned sequences in BAM format\n"

	echo -e "REQUIRES: The script requires samtools (v1.2)"
	exit
fi


module load samtools/1.2

bam=$1


samtools depth $bam | awk '{sum+=$3} END { print "Average = ",sum/NR}' > $bam.cov



