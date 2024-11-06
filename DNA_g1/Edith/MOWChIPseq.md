#Introduction


#!/bin/bash

############################# MOWChIP Script #################################
# For proper script usage, place all fastq or fastq.gz files in a folder
# named Raw_Data within your directory of interest. Other folders will be
# generated within that same directory alongside the Raw_Data folder
#
# Please alter the input variables to match your parameters
# NOTE: If you are running this in a job queue, you may need to fetch the
# chromosome sizes manually.
##############################################################################

############################# Input Variables ################################
# Please fill in these variables with your specific parameters
##############################################################################

# This directory should contain the Raw_Data folder that holds the fastq files
PWD=Enter_Path_To_Directory

# Enter the name of the input file without the fastq extension 
# i.e. input-2.fastq is input-2
input_base=Input_File_Without_Extension 

# This is the path to the Bowtie genome files.
# See bowtie documentation for more details
BWT=Directory_Of_Bowtie_Genome

# Enter the name of the genome you are using (e.g. hg19, mm9)
# This is used by SICER and for naming files
genome=Name_of_genome_to_use

# Enter 1 if human, 0 if mouse. For other organisms, alter MACS2 section
human=1


########################## Unzipping and Trimming ############################
# Files are unzipped if necessary and low quality reads are removed
##############################################################################

cd $PWD/Raw_Data
gunzip *.gz
trim_galore *.fastq
cd ..

############################# EXTENSION SETUP ################################
# Sets up some common file extensions that will be used throughout the script
##############################################################################

FILES=$PWD/Raw_Data/*.fq
SAM=.sam
BAM=.bam
TRIM=_trimmed.fq
WIN=_window.bed
WINSORT=_window_sorted.bed
SORT=_sorted.bed
NOR=_normalized.bedGraph
BED=_Aligned.bed
EXT=_Extended.bed
FULLNOR=_Ext_Normalized.bedGraph
EXTWIN=_Extended_Window.bed
EXTSORT=_Extended_Sort.bed
EXTBW=_Ext_Normalized.bw
SICER=_SICER.bed

############################### FOLDER SETUP #################################
# Sets up the folders that will contain processed files
##############################################################################

mkdir Aligned_SAM #All aligned reads
mkdir Aligned_BAM #Unique aligned reads
mkdir Aligned_BED #All bed files (unique, extended, windowed)
mkdir Ext_Normalized_WIG #Normalized extended 100nt windows
mkdir Ext_Normalized_BW #Normalized extended 100nt windows. Use to visualize.
mkdir Normalized #Normalized 100nt windows
mkdir Input #Processed Input file
mkdir MACS2 #Narrow Peaks
mkdir SICER #Broad Peaks
mkdir Genome #Genome chromosome sizes. Used by script.

############################### MAKE WINDOWS #################################
# Gets the sizes of the chromosomes for genome of interest
# Makes 100nt windows of the genome of interest
# If using a job queue, fetchChromSizes might need to be performed manually
##############################################################################

fetchChromSizes $genome > $PWD/Genome/Genome-"$genome".bed
bedtools makewindows -g $PWD/Genome/Genome-"$genome".bed -w 100 > $PWD/Genome/GenomeWindowed-"$genome".bed

######################## PREPPING THE INPUT FILE #############################
#Aligns the input file to the genome
#Isolates unique bam reads
#Find the 100nt windows for correlation and sort it
#Extend the input 100nt on either side of the sequence then 100nt count windows
#Get the input_length for normalization
##############################################################################

bowtie -t -p 16 -S $BWT $PWD/Raw_Data/$input_base$TRIM -S $PWD/Input/Input.sam
samtools view -@ 16 -bq 1 $PWD/Input/Input.sam > $PWD/Input/Input.bam

bedtools coverage -counts -b $PWD/Input/Input.bam -a $PWD/Genome/GenomeWindowed-"$genome".bed > $PWD/Input/Window_Input.bed

bedtools bamtobed -i $PWD/Input/Input.bam > $PWD/Input/Input.bed
bedtools slop -i $PWD/Input/Input.bed -g $PWD/Genome/Genome-"$genome".bed -b 100 > $PWD/Input/Input_Extended.bed
bedtools coverage -counts -b $PWD/Input/Input_Extended.bed -a $PWD/Genome/GenomeWindowed-"$genome".bed > $PWD/Input/Window_Extended_Input.bed

input_length=$(samtools view -c "$PWD/Input/Input.bam")

####################### PROCESSING THE SAMPLE FILE ###########################
#Get the input length (for normalization)
#Loop through each sample. 
#Get the core name of the file
#Align the sample to the genome
#Obtain the unique bam file
##############################################################################

for fn in $FILES
do

#This pulls the extension off of the file name
f=`basename "$fn" _trimmed.fq` 

#This makes sure we aren't processing our input as a sample
if [ "$f" = "$input_base" ]; then
continue
fi

echo "Processing $f"

bowtie -t $BWT -p 16 -S $PWD/Raw_Data/$f$TRIM $PWD/Aligned_SAM/$f$SAM
samtools view -@ 16 -bq 1 $PWD/Aligned_SAM/$f$SAM > $PWD/Aligned_BAM/$f$BAM
bedtools bamtobed -i $PWD/Aligned_BAM/$f$BAM > $PWD/Aligned_BED/$f$BED

############################## Normalization  ################################
# This normalization can be used for correlation
# Get the length of the sample file for normalization
# Obtain 100nt windows for the sample and sort it
# Normalize the sample promotor file against input promotor file
##############################################################################

ChIP_length=$(samtools view -c $PWD/Aligned_BAM/$f$BAM)
bedtools coverage -counts -b $PWD/Aligned_BAM/$f$BAM -a $PWD/Genome/GenomeWindowed-"$genome".bed > $PWD/Aligned_BED/$f$WIN
sort -k1,1 -k2,2n -o $PWD/Aligned_BED/$f$WINSORT $PWD/Aligned_BED/$f$WIN
paste $PWD/Aligned_BED/$f$WINSORT $PWD/Input/Window_Input.bed | awk -v OFS="\t" '{print $1,$2,$3,$4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Normalized/$f$NOR

################################### MACS2 ####################################
# MACS2 for narrow peak calling (e.g. most transcription factors and some
# histone modifications)
#
# Checks if human flag raised to call appropriate form of MACS2
# To call other species, check MACS2 documentation and modify
# Uses a FDR cutoff of 0.05 to determine peaks
##############################################################################

if [ "$human" = "1" ]; then
	macs2 callpeak -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g hs -n $f -q 0.05 --outdir $PWD/MACS2
else
	macs2 callpeak -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -f BED -g mm -n $f -q 0.05 --outdir $PWD/MACS2
fi	

################################### SICER ####################################
# epic's implementation of SICER is for broad peak calling 
# Please check epic documentation for other acceptable genomes
##############################################################################

epic -t $PWD/Aligned_BED/$f$BED -c $PWD/Input/Input.bed -gn $genome --window-size 1000 --gaps-allowed 3 -o $PWD/SICER/$f$SICER

################################ Visualization ###############################
# Use the BigWig files in Ext_Normalized_BW for viewing in IGV
# Extends the reads by 100nt on either side before binning into 100nt windows
# Windows are sorted and normalized against the input
# They are sorted again and converted to BigWig format
##############################################################################

bedtools slop -i $PWD/Aligned_BED/$f$BED -g $PWD/Genome/Genome-"$genome".bed -b 100 > $PWD/Aligned_BED/$f$EXT
bedtools coverage -counts -b $PWD/Aligned_BED/$f$EXT -a $PWD/Genome/GenomeWindowed-"$genome".bed > $PWD/Aligned_BED/$f$EXTWIN
sort -k1,1 -k2,2n -o $PWD/Aligned_BED/$f$EXTSORT $PWD/Aligned_BED/$f$EXTWIN
paste $PWD/Aligned_BED/$f$EXTSORT $PWD/Input/Window_Extended_Input.bed | awk -v OFS="\t" '{print $1,$2,$3,$4/'$ChIP_length'*1000000-$8/'$input_length'*1000000}' > $PWD/Ext_Normalized_WIG/$f$FULLNOR
bedSort $PWD/Ext_Normalized_WIG/$f$FULLNOR $PWD/Ext_Normalized_WIG/$f$FULLNOR
bedGraphToBigWig $PWD/Ext_Normalized_WIG/$f$FULLNOR $PWD/Genome/Genome-"$genome".bed $PWD/Ext_Normalized_BW/$f$EXTBW

##############################################################################

done
