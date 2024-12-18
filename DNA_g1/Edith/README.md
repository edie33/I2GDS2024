<img src=https://github.com/edie33/I2GDS2024/blob/2212aa7db42493e201af9c4fe8e528271acde41b/DNA_g1/Edith/images/pipeline.png width=65% height=40%>


# MOWChIPseq for Histone Modification Profiling

## Introduction
This repository explains a general pipeline for data analysis with MOWChIP-Seq samples in a linux HPC environment. 
It was developed as part of the curriculum for Virginia Tech's ALS 5224 Intro to Genomic Data Science course. 

Epigenetic alterations, such as histone modifications, can help us understand how gene expression is regulated, how environmental factors can impact health, develop and evaluate  epigenetic therapeutic agents.

In order to identify significant histone modifications in samples and observe patterns across patients with and without treatment, we obtained and prepared samples with MOWChIP-seq. However, samples may be used from previously published data for running and testing the data analysis process. 

Contact: Edith Chen (edithchen@vt.edu)

# 1. Setup/Preparation

## 1.1 Create directories
For each project, it is important to have your own directory for the samples and scripts.

```
mkdir MOWChIP_practice
cd MOWChIP_practice

mkdir scripts
mkdir results
mkdir rawdata

cd rawdata
mkdir input
mkdir chip
```
The input samples are samples that has not been while the chip samples are samples with MOWChIP.

## 1.2 Create conda environment
A conda environment is a directory that contains a specific set of software packages, including libraries, dependencies, and Python versions.
This enables easy access to the specific versions of software packages needed for a project. For ChIP-seq analysis, I utilised the file ```spec_environment.txt```.

To create an environment:

```bash
conda create --name <my-env>
```
Replace <my-env> with the name of your environment.

For downloading with the file of the list of packages, replace <this file> with the name of the file.

```bash
conda create —name <env> python=3.7.12 --file <this file>
conda create —name chipseq python=3.7.12 --file `spec_environment.txt
```

For further information on managing conda environments, please visit https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.

## 1.3.1 Download Sample Data
The demo data set is from the Zhu et al., Nature Protocols 2019. can be found at: 

For information downloading the sample data from the paper
Each data sample used in papers will be accompanied with a GEO accession number.
Search for the GEO number/code on the GEO website: http://www.ncbi.nlm.nih.gov/geo/
For this practice, the GEO code is GSE123606.

Then a page for that Series will load with a ‘Samples’ section. Click the ‘More’ link if necessary to see all the samples in the entry. 
Obtain the SRR numbers for each sample and search for it on the European Nucleotide Archive (EBI) SRA page: http://www.ebi.ac.uk/ena/data/view/SRR_number.

On the following page, look for the fastq files (ftp) link and right-click on the link. A pop-up menu will appear and select 'Copy Link.:

Retrieve the file under designated folder using wget:
```
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR494/SRR494099/SRR494099.fastq.gz
```

For further information, please visit https://www.imm.ox.ac.uk/files/ccb/downloading_fastq_geo

## 1.3.2 Original Data

If you are using your own samples instead of those from papers, directly upload the data into the designated folder.
If demultiplexing is needed, use the following code:
```
zgrep  --no-group-separator -A 3 +<index sequence> /path-to-current-directory/ /path-to-designated-sample-directory/filename
```

Examples: 
```
zgrep --no-group-separator -A 3 +TAAGGCTC /projects/lu_lab/All_Raw_Sequencing_Data/20240613_Novogene_GL_XZ_EC/Undetermined_1.fq.gz>/projects/lu_lab/Gaoshan/scChIPs/bulk/H3K4me3/0.5Tn5_barcode_oligo_R1.fastq

U7=CAAGCTAGATCT +CGCTATGT
zgrep --no-group-separator -A 3 +$U7 /projects/lu_lab/All_Raw_Sequencing_Data/20240822_Novogene_GL_JN_XZ_EC/undetermined_1.fq.gz>/projects/lu_lab/edith/Practice/MOWChIP/rawdata/chip/ChIP1_0628_R1_5.fastq
```


# 2. Alignment
## 2.1 Index genome and download additional packages

1. Index genome of interest
2. Obtain blacklist of genes from genome of interest (i.e., repetitive regions)
3. fetchChromSizes for genome of interest (*genome.bed file)
4. bedtools makewindows to separate genome into windows of 100 bp (*genome_100.bed file)
5. Generate Promoters with promoter.sh (*Promoter.bed file) (for more info read promoter.sh script)

## 2.2 Align input samples
Use ```bowtie2``` to align the input files to the genome index.
This is FM-index based on Burrows-Wheeler transform.
It is memory-efficient but slower for longer reads 

```
bowtie2 -p 16 -x /home/gaoshanli/Data/mm10/Sequence/Bowtie2Index/genome -U $PWD/Raw_Data/$f$FQ -S $PWD/Aligned_SAM/$f$SAM 2>$PWD/Aligned_SAM/$f$LOG
```
The aligned reads are saved in `input.sam`, and any alignment errors are recorded in `alignSumm.log`.
For further information on bowtiw2, visit: https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

The ```preinput.sh``` file is a script that uses bowtie2 for alignment and performs downstream filtering and analysis.
After converting the SAM files into BAM files, the alignments are filtered and converted into BED format. This is followed by filtering with blacklist, extending the feature sizes, and completing coverage analysis.
```
sbatch preinput.sh
```
Use ```squeue``` to check the status of the processing on linux.

# 3. MACS2 for ChIP samples
MACS2 is a program for peak calling; specifically it is computational method to identify areas in the genome that have been enriched with aligned reads.


```
macs2 callpeak -t $PWD/BED/$f$BED -c /projects/intro2gds/DNAGroup/ChIP-seq/rawdata/edith/input/input.bed -f BED -g mm -n $f -q 0.05 --outdir $PWD/MACS
```
For further information please visit: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html

```
sbatch base_precorr.sh
```

# 4. IGB
The Integrated Genome Browser (IGB) can be used to explore & visually analyze vast genomic data and biologically-interesting patterns.
To download IGB, visit: https://www.bioviz.org/

<img src=https://github.com/edie33/I2GDS2024/blob/0b34aa4bcbbfbfb84e00789bd7d94a8efd51bc4e/DNA_g1/Edith/images/igb1.png width=65% height=40%>

Choose species and genome version. To choose a species and genome version. ...
Open data sets. Open data sets from remote data sources (Data Access tab) or by opening local files. ...
Zoom in. ...
Load data. ...
Configure tracks.
<img src=https://github.com/edie33/I2GDS2024/blob/0b34aa4bcbbfbfb84e00789bd7d94a8efd51bc4e/DNA_g1/Edith/images/igb2.png width=65% height=40%>

<img src=https://github.com/edie33/I2GDS2024/blob/172ac6b81c2696d4bea1f83c07414a12f0618856/DNA_g1/Edith/Images/mowchip_pearson.png width=65% height=40%>
**Figure Source:** [Nature Protocols, 2019](https://www.nature.com/articles/s41596-019-0223-x)
