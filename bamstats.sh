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

	echo -e "OUTPUT:	<bamtools>	A bamtools output"
	echo -e "		<qualimap>	A qualimap output\n"

	echo -e "REQUIRES: The script requires java (v14.0.2), bamtools (v2.5.1), and qualimap (v2.3)"
	exit
fi



# Export path to Qualimap v2.3
export PATH=/cluster/work/pausch/cbortoluzzi/softwares/qualimap_v2.3:$PATH


# Load modules required to run the script
module load bamtools/2.5.1
module load openjdk/14.0.2


bam=$1
outbam=$(basename $bam)


mkdir -p bamstats
if [[ ! -f "bamstats/$outbam.bamtools.stats" ]];then
	echo -e "Run bamtools to obtain a set of general statistics for $bam\n\n"
	bamtools stats -in $bam > bamstats/$outbam.bamtools.stats
fi


if [[ ! -f "bamstats/$outbam/genome_results.txt" ]];then
	echo -e "Run qualimap to obtain a set of general statistics for $bam...this will take some time\n\n"
	qualimap bamqc --java-mem-size=32G -bam $bam -nt 16 -nw 500 -outdir bamstats/$outbam
fi

