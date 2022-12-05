## Setup

# Run from the server
mkdir micb405_final_project && cd micb405_final_project

# Due to problems with SRA toolkit on the ORCA1 server (which we can't do much about as students), try downloading from ENA using wget
# Hover over "generate fastq files" to see link; right click to copy it. Then rename to something more clear
# Check md5sum by going to column selection and checking fastq_md5

## Get all necessary data
# https://www.ebi.ac.uk/ena/browser/view/SRX1023255
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/008/SRR2016188/SRR2016188.fastq.gz
md5sum SRR2016188.fastq.gz # should be c5c51ff4a969c4513250c71181b120df
# Rename
mv SRR2016188.fastq.gz ./ctrl_SRX1023255.fastq.gz

# https://www.ebi.ac.uk/ena/browser/view/SRX1023254
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/007/SRR2016187/SRR2016187.fastq.gz
md5sum SRR2016187.fastq.gz # should be 27d45a96ada5a438ecca37eb33dfc377
mv SRR2016187.fastq.gz ./ctrl_SRX1023254.fastq.gz

# https://www.ebi.ac.uk/ena/browser/view/SRX1023258
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/001/SRR2016191/SRR2016191.fastq.gz
md5sum SRR2016191.fastq.gz # should be 1d04d664ebe44ce31663a8fa5c221155
mv SRR2016191.fastq.gz ./treat_SRX1023258.fastq.gz

# https://www.ebi.ac.uk/ena/browser/view/SRX1023257
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR201/000/SRR2016190/SRR2016190.fastq.gz
md5sum SRR2016190.fastq.gz # should be 90346caba85865cd3b58fd12a275d11f
mv SRR2016190.fastq.gz ./treat_SRX1023257.fastq.gz

# Generate FASTQC for the 4 RNA read files
fastqx *.fastq.gz

# Get files from the original study. Run this once.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/595/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna.gz

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/595/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.gtf.gz

# Check file integrity by checking the md5checksums.txt from https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/595/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5/
# Run once
md5sum GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.*.gz

# Unzip files
gunzip *.gz

## STAR alignment

# make directory
mkdir star_index

# Indexing ref genome
STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna --sjdbGTFfile GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.gtf --runThreadN 8

# Alignment (Run 4 times)
# https://www.biostars.org/p/331735/ single end read
# And of course, replace sequences_R1.fastq.gz with whichever of the 4 fastq.gz files you are aligning
STAR --genomeDir star_index --readFilesIn <sequences_R1.fastq.gz> --outFileNamePrefix ./star_alignments/<sequences_R1> --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --runThreadN 8

# Creating the project_essentials directory
mkdir project_essentials && cd project_essentials
mkdir star_tabs
mkdir qc_htmls

# Move needed files around
scp ../*.out.tab ./star_tabs # The tab files (which we will convert to tsv)
scp ../*.html ./qc_htmls # The fastQC html files

# Index the bam files
for bam in $(ls *.bam); do samtools index $bam; done

# Import to your user computer (use a new command line, from home computer)
scp -r <username>@orca1.bcgsc.ca:/home/<username>/micb405_final_project/project_essentials/* ./Documents/
