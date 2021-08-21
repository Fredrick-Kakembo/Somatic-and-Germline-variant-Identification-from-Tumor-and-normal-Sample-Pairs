<h1 align="center"> Identification of somatic and germline variants from tumor and normal sample pairs </h1>
<h1 align="center"><pre>Reproduced By Genomics Two A</pre></h1>


## `Introduction` <a name="introduction"></a>

Mutations (random single or multiple base changes) in DNA or RNA can have a beneficial (eg in evolution), neutral or harmful effect in an organism. Many diseases including mostly **cancers** (second leading cause of death) are as a result of harmful mutations in crucial genes eg Tumor suppressor genes, that cause cells to grow and divide uncontrollably, infiltrating and destroying normal body tissues. These mutations can be germline (inhrited) or somatic (acquired after birth), and a common kind of genetic mutation as a result of either is [Loss of Heterozygosity (LOH)](https://en.wikipedia.org/wiki/Loss_of_heterozygosity). LOH usually leads to loss of one normal copy or a group of genes, which is a common even in cancer development. Germline mutations can easily be identified by comparing a sample genome to a reference, however the story is quite different when it comes to somatic mutations as we need both a normal and tumor tissue DNA from the patient. 

<br>

> In this project, we aimed at reproducing a workflow that identifies germline and somatic variants, variants affected by LOH using both a health and tumor tissue, from which we would report variant sites and genes affected that could likely be the cause to the disease. Such insights can help us track the genetic events driving tumorigenesis in patients and might be useful in diagnosis, prognosis, developing and guiding therapeutics strategies.

<br>

Below is our Graphical abstract summarizing the key steps we took to achieve this.

![Graphical Abstract](Graphic_Abstract-Genomics-Two-A.png)
Figure 1: Graphical Abstract for Genomics-Two-A

<br>

We reproduced this tutorial both as a Galaxy Tutorial as well as Linux Pipeline.


### `Go To Section:`

1. [Introduction](#introduction)
2. [Section One: Linux Pipeline](#linux)
3. [Section Two: Galaxy Workflow](#galaxy)
5. [Contributors](#contributor)


# Section One: `Linux Pipeline` <a name="linux">.</a>

## Dataset Description

The datasets used in this analysis (reads from human chromosomes 5, 12 and 17), were obtained from a cancer patient’s tumor and normal tissue samples. The normal tissue coudn't be the only sample used because healthy tissue contains many variants and every individual inherits a unique pattern of many variants from their parents. The samples (paired end) were two in number.
	
## Data download
The datasets were downloaded from Zenodo using the wget command.
	
### Samples dataset

```
echo -e "\n Downloading data... \n"
	
mkdir -p raw_data 
cd raw_data
	
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-N_231335_r2_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r1_chr5_12_17.fastq.gz
wget https://zenodo.org/record/2582555/files/SLGFSK-T_231336_r2_chr5_12_17.fastq.gz	
```

### Reference sequence

```
echo -e "\n Downloading reference sequence... \n"
	
wget https://zenodo.org/record/2582555/files/hg19.chr5_12_17.fa.gz

#unzip reference
unzip hg19.chr5_12_17.fa.gz
```

## Pre-processing and Trimming

### i) Quality check.

The reads quality were examined using fastqc and an aggregate report generated with multiqc.

#### Description

FastQC aims to provide a way to do quality control checks on sequence data. Within the `fastq` file is quality information that refers to the accuracy of each base call. This helps to determine any irregularies or features that make affect your results such as adapter contamination

#### Installation

```{bash}
conda install -c bioconda fastqc multiqc --yes
```
#### Command
	
```{bash}
echo -e "\n Data Preprocessing... \n"

mkdir -p Fastqc_Reports  #create directory for the fastqc output
```
```
#Qc on reads
for sample in `cat list.txt`
do
	fastqc raw_data/${sample}*.fastq.gz -o Fastqc_Reports
done

multiqc Fastqc_Reports -o Fastqc_Reports	
```

The multiqc report can be examined from [here](multiqc_report_linux.html). 
From the report, the reads quality are great, a few adapters are however observed.

### ii) Removing low quality sequences using Trimmomatic
<http://www.usadellab.org/cms/?page=trimmomatic>

#### Description

`Trimmomatic` is a wrapper script that automate quality and adapter trimming. After analyzing data quality, the next step is to remove sequences that do not meet quality standards. 	

#### Installation
```	
conda install -c bioconda trimmomatic --yes
```	

#### Command
```
mkdir -p trimmed_reads

for sample in `cat list.txt`
do
       trimmomatic PE -threads 8 raw_data/${sample}_r1_chr5_12_17.fastq.gz raw_data/${sample}_r2_chr5_12_17.fastq.gz \
               trimmed_reads/${sample}_r1_paired.fq.gz trimmed_reads/${sample}_r1_unpaired.fq.gz \
               trimmed_reads/${sample}_r2_paired.fq.gz trimmed_reads/${sample}_r2_unpaired.fq.gz \
               ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:keepBothReads \
               LEADING:3 TRAILING:10 MINLEN:25
       
       fastqc  trimmed_reads/${sample}_r1_paired.fq.gz  trimmed_reads/${sample}_r2_paired.fq.gz \
                 -o trimmed_reads/Fastqc_results
done 

multiqc  trimmed_reads/Fastqc_results  -o trimmed_reads/Fastqc_results

```
	
The parameters shown below were used during trimming:
* PE - paired end
* threads - number of cores assigned to the task, 
* LEADING - remove lead bases with low quality of 3
* TRAILING - remove trailing bases with low quality of 10
* MINLEN - remove reads below 25 bases long
* ILLUMINACLIP - used to remove adapters
	* _Trused3-PE - adapter_, _2 - Maximum mismatch count_, _30 - Accuracy of the match between the two ‘adapter ligated’ reads for PE palindrome read alignment_, _10 - Accuracy of the match between any adapter against a read_, _8 - Minimum length of adapter that needs to be detected (PE specific/ palindrome mode_
	
	
The post trimming multiqc report can be found [here](post_trim_multiqc_report_linux.html). It is evident from the report that the quality of the reads improved having per base quality scores above 35 and no adapters observed. After trimming an average of 0.73% normal reads and 1.24% tumor reads were lost.

**NB: To view the multiqc html reports download the files and view them from your browser.**

	
## Mapped read postprocessing
	
### Description	
Mapping of sample sequences against the reference genome is conducted with an aim of determining the most likey source of the observed sequencing reads.
`BWA-MEM` was used for alignment. The results of mapping is a sequence alignment map (SAM) format. The file has a single unified format for storing read alignments to a reference genome.

### Installation 
```
conda install -y -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda bamtools

```

### Command
	
#### Read mapping
	
In order to align the data, we need a reference to align against.  First, a directory is created for the reference and then copied. The reference is  indexed to be able to align the data.This is done using the command;
	
``` 
bwa index hg19.chr5_12_17.fa 
```
	
This produces 5 files in the reference directory that BWA uses during the alignment phase. The 5 files have different extensions named amb,ann,bwt pac and sa. Alignment can be done using the command;
	
```
Bwa mem
```
	
Note that  bwa is given a location , which is the path to the reference. Now, the two paired-end files are aligned and the alignment output (in SAM format) directed to a file. 24 threads (processors) were used to speed up this process and a read group (i.e sample ID) information was added to the alignment:
	
 ```
 bwa mem -t 24 -R '@RG\tID:231335\tSM:Normal' ./reference/hg19.chr5_12_17.fa.gz SLGFSK-N_231335_r1_chr5_12_17.fastq.gz  SLGFSK-N_231335_r2_chr5_12_17.fastq.gz >SLGFSK-N_231335_paired.sam
	
 bwa mem -t 24 -R '@RG\tID:231336\tSM:Tumor' ./reference/hg19.chr5_12_17.fa.gz SLGFSK-T_231336_r1_chr5_12_17.fastq.gz SLGFSK-T_231336_r2_chr5_12_17.fastq.gz >SLGFSK-T_231336_paired.sam
 
 ```
	

	
#### Conversion of the SAM file to BAM file, sorting and indexing
A Binary Alignment Map (BAM) format is an equivalent to sam but its developed for fast processing and indexing. It stores every read base, base quality and uses a single conventional technique for all types of data.
The produced BAM files were sorted by read name and indexing was done for faster or rapid  retrieval. At the end of the every BAM file,  a special end of file (EOF) marker is usually written, the samtools index command also checks for this and produces an error message if its not found.
	
```
for sample in `cat list.txt`
do
        Convert SAM to BAM and sort it 
        samtools view -@ 20 -S -b Mapping/${sample}.sam | samtools sort -@ 32 > Mapping/${sample}.sorted.bam
        
        Index BAM file
        samtools index Mapping/${sample}.sorted.bam
done
```	

#### Mapped reads filtering
	
```
for sample in `cat list.txt`
do
	#Filter BAM files
        samtools view -q 1 -f 0x2 -F 0x8 -b Mapping/${sample}.sorted.bam > Mapping/${sample}.filtered1.bam
done
```

To view the output of the results use :
```
samtools flagstat <bam file>
```

#### Duplicates removal
During library construction sometimes there's introduction of PCR (Polymerase Chain Reaction) duplicates, these duplicates usually can result in false SNPs (Single Nucleotide Polymorphisms), whereby the can manifest themselves as high read depth support. A low number of duplicates (<5%) in good libraries is considered standard.

```
#use the command <markdup>
for sample in `cat list.txt`
do
	samtools collate -o Mapping/${sample}.namecollate.bam Mapping/${sample}.filtered1.bam
        samtools fixmate -m Mapping/${sample}.namecollate.bam Mapping/${sample}.fixmate.bam
        samtools sort -@ 32 -o Mapping/${sample}.positionsort.bam Mapping/${sample}.fixmate.bam
        samtools markdup -@32 -r Mapping/${sample}.positionsort.bam Mapping/${sample}.clean.bam
done
	
#or <rmdup>
samtools rmdup SLGFSK35.sorted.bam  SLGFSK35.rdup and samtools rmdup SLGFSK36.sorted.bam  SLGFSK36.rdup.
```

#### Left Align BAM
```
for sample in `cat list.txt`
do
        #-c -> compressed, -m -> max-iterations
        cat Mapping/${sample}.clean.bam  | bamleftalign -f hg19.chr5_12_17.fa -m 5 -c > Mapping/${sample}.leftAlign.bam
```

#### Recalibrate read mapping qualities
```
for sample in `cat list.txt`
do
        samtools calmd -@ 32 -b Mapping/${sample}.leftAlign.bam hg19.chr5_12_17.fa > Mapping/${sample}.recalibrate.bam
done
```
	
#### Refilter read mapping qualities

```
for sample in `cat list.txt`
do
        bamtools filter -in Mapping/${sample}.recalibrate.bam -mapQuality "<=254" > Mapping/${sample}.refilter.bam
done
```

## Variant calling and classification
<http://varscan.sourceforge.net/somatic-calling.html>
	
### Description	
To be able to identify variants from the mapped samples, the tool `VarScan somatic` was used. 
The command expects both a normal and tumor sample in `Samtools pileup` format and outputs an indel file and snp file.
The command reports germline, somatic, and LOH events at positions where both normal and tumor samples have sufficient coverage 

### Installation 
```
wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar		
```

### Command
#### Convert data to pileup

```
mkdir Variants

for sample in `cat list.txt`
do
        samtools mpileup -f hg19.chr5_12_17.fa Mapping/${sample}.refilter.bam --min-MQ 1 --min-BQ 28 \
                > Variants/${sample}.pileup
done
```

#### Call variants
```
java -jar VarScan.v2.3.9.jar somatic Variants/SLGFSK-N_231335.pileup \
        Variants/SLGFSK-T_231336.pileup Variants/SLGFSK \
        --normal-purity 1  --tumor-purity 0.5 --output-vcf 1 
```

#### Merge vcf
VarScan generates 2 outputs (indel.vcf and snp.vcf), merge the two into one vcf file using `bcftools.`
```
#merge vcf
bgzip Variants/SLGFSK.snp.vcf > Variants/SLGFSK.snp.vcf.gz
bgzip Variants/SLGFSK.indel.vcf > Variants/SLGFSK.indel.vcf.gz
tabix Variants/SLGFSK.snp.vcf.gz
tabix Variants/SLGFSK.indel.vcf.gz
bcftools merge Variants/SLGFSK.snp.vcf.gz Variants/SLGFSK.indel.vcf.gz > Variants/SLGFSK.vcf
```
		
## Variant Annotation
### Functional Annotation using `SnpEff`
<https://pcingola.github.io/SnpEff/examples/>
		
#### Description
`SnpEff` is a variant annotator and functional effect predictor. The output is appended to the vcf file with the field `ANN`. A snpEff database is required prior to performing annotation. In case the organism of interest is not present in the snpEff database, you can build the database using the snpEff command. If the organism is present in the database, download it using the snpEff command.

#### Installation 
```
#download jar file
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

# Unzip file
unzip snpEff_latest_core.zip
		
#download snpEff database
java -jar snpEff.jar download hg19		
```

#### Command
```
#annotate variants
java -Xmx8g -jar snpEff/snpEff.jar hg19 Variants/SLGFSK.vcf > Variants/SLGFSK.ann.vcf
```		
		
### Clinical Annotation using `gemini`
<https://gemini.readthedocs.io/en/latest/content/preprocessing.html>
		
#### Description
		
		
#### Installation
```
wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py
python gemini_install.py /usr/local /usr/local/share/gemini
```
		
#### Command
```
gemini load -v Variants/SLGFSK.ann.vcf -t snpEff Annotation/gemini.db
		
```
	
	
# Section Two:  `GALAXY WORKFLOW` <a name="galaxy">.</a>

<Lets add the galaxy sections here>

## 1. Data Preparation:

The sequencing reads that were used for analysis were obtained from a cancer patient's normal and tumor tissues.
There were a total of four samples. A forward reads sample and a reverse reads sample was obtained for both the normal and tumor tissue. A human reference genome, hg19 version was also used for analysis.

The first step was to create a new workflow we named Genomics_2_A in the galaxy window. We then imported the 4 fastq files and the reference using the `Upload` button then `Paste/Fetch Data` and pasting the corresponding links to the data, selecting datatype as `fastqsanger.gz` for the fsatq samples and `fasta` for the reference. Once done uploading, the attributes of the samples were edited using the pencil mark to **tumor** and **normal** correspondily for easier identification.

## 2.  Quality Control & Check:

•	FastQC:  is a quality control tool for high throughput sequence data that gives a summary report about the sequence.

•	MultiQC: A modular tool to aggregate results from bioinformatics analyses across many samples into a single report.


The MultiQC Output & Report:

For the quality control analysis, the following figures show that the chosen dataset is of high quality for both Normal R1 & R2 datasets and both Tumor R1 & R2 datasets even the tumor ones are of poorer quality than the normal ones:
o   All nucleotides have high-quality scores, (the forward and reverse reads of normal and tumor patient’s tissues), as they all are present in high/good quality region.

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/239887479_10225646820098581_7852584000625137608_n.jpg?_nc_cat=103&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeEGg5POsRQ22z30bSKJKOHMNyXKaVeIx-Q3JcppV4jH5H7955fLdsp5_S5xl6vVh4A&_nc_ohc=qsVI3WWfM6MAX-9V1fz&_nc_ht=scontent.fcai20-3.fna&oh=0f87104f5f97008322df582eea23119f&oe=61249517)

 o	Good quality score distribution as the mean is high with a sharp, distinct peak.

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/240386199_10225646820338587_8217938469652979046_n.jpg?_nc_cat=107&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeHfYmVb7PRpPMCZZSub8LeiagTNa3_Pz8lqBM1rf8_PyQOIPyzQKNJrrF5grXk5bpU&_nc_ohc=G1o6LS9X_hYAX-478Tr&_nc_ht=scontent.fcai20-3.fna&oh=e2a9dbde12186e84f979a5391d58fa84&oe=6124BBAD)

o	The actual mean of the GC% content is lower than the theoretical one and that is a non-normal distribution which may indicate some contamination; however, this peculiar bimodal distribution is considered to be a hallmark of the captured method as using Agilent’s SureSelect V5 technology for exome enrichment.

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/240332287_10225646820538592_1031475029124752916_n.jpg?_nc_cat=103&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeF1xpsEQ-0f3Eb3Mqxtc-Npz8gRJUkIxxLPyBElSQjHEp63raX5oUlFIBOhFqte0o4&_nc_ohc=3X7mee3XrsMAX8VAWOt&_nc_ht=scontent.fcai20-3.fna&oh=dadab37eda7bc7095539a18260709433&oe=61251953)

o	The N content is zero, thus not indicating any bad base detection.

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/239932544_10225646821178608_1804771247842536342_n.jpg?_nc_cat=111&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeFJ5LZYvlO_zgrhPl4Qvqy0RC3ghuo8qpxELeCG6jyqnLVOoaNAXDmrf5lC-FTyqoA&_nc_ohc=TsZziKppm1gAX8l1FUW&tn=D8rHGaIJCwa8Og3I&_nc_ht=scontent.fcai20-3.fna&oh=3753222c1f2642ed1b8f013d8ae71218&oe=6123FF0F)

o	Sequence duplication level is good as well with almost no PCR biases when library was prepared as there is no over-amplified fragment.

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/240452151_10225646820978603_4722043033758089712_n.jpg?_nc_cat=103&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeFO1KjcfjwGSAM-xaMakh3a9bVbKhLaQoT1tVsqEtpChIiRWFZIA79Wijx9qGdCRHA&_nc_ohc=EVOa3QcTTFUAX-an6Ea&tn=D8rHGaIJCwa8Og3I&_nc_ht=scontent.fcai20-3.fna&oh=15a70d18a95a3ab36bec20f483ece380&oe=612406E1)

o	There is a little presence of adaptors in the sequence.

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/240066890_10225646821418614_9064760851692879195_n.jpg?_nc_cat=107&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeHvZIsBRpPCCRbJELOa_nlbsKsWUcxpgYmwqxZRzGmBid0AqMMgD3B6gqZWkCywp5c&_nc_ohc=fjLSS02SO9wAX9MNlaP&_nc_ht=scontent.fcai20-3.fna&oh=4fa0d302f28428aa9f67700fa89854d6&oe=612565FC)



## Read Trimming and Filtering

The next process after quality control is Trimming and Filtering. This process further helps to trim and filter raw reads to give a more improved quality.

° Tool - TRIMMOMATIC (Galaxy Version 0. 38. 0) is a fast, multithreaded command line tool that can be used to trim and crop low quality reads and remove the present adaptors in the reads to improve the quality of these processing datasets.

° Input data - The forward read FASTQ file (r1) and the reverse read FASTQ file (r2) of the normal tissue were run concurrently as paired-end by performing initial ILLUMINACLIP step.

° Parameters - 
![Parameters 1](https://github.com/Fredrick-Kakembo/Somatic-and-Germline-variant-Identification-from-Tumor-and-normal-Sample-Pairs/blob/Trimmomatic-Parameters/IMG_20210820_230836.jpg)
![Parameters 2](https://github.com/Fredrick-Kakembo/Somatic-and-Germline-variant-Identification-from-Tumor-and-normal-Sample-Pairs/blob/Trimmomatic-Parameters/IMG_20210820_231013.jpg)

° Output data - Trimmed forward and reverse reads for each normal and tumor tissue. 
            - Orphaned forward and reverse reads for each normal and tumor tissue which the corresponding mate got dropped because of insufficient length after trimming. These datasets are empty and therefore deleted. 



° The MultiQC Output & Report after Trimming:

The quality reads did not change much as the datasets were already of high quality although a small fraction of adapter is successfully removed as shown in the following report:

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/239936305_10225646826698746_7439285328275684759_n.jpg?_nc_cat=104&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeEHowE9hJvIdP00B3bXEG5QGNS2Mrs4sTEY1LYyuzixMaJmAJ4zP9daBQiwILofdJ8&_nc_ohc=tUKGyhg8GkgAX8GLGoh&_nc_ht=scontent.fcai20-3.fna&oh=9a6524fc18a614b7d92f3cd04231ef24&oe=6124AFE5)

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/240397072_10225646826898751_7144027909249758575_n.jpg?_nc_cat=100&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeGiMN88ld5QT1-kVw42WldImmCThDU_8KmaYJOENT_wqY4BeSIlAskfw6yiEsQrkqw&_nc_ohc=GsLoATZazrIAX9uuOWd&_nc_oc=AQkvTp0FpyMonETpOoyRdYs3gDhNhZRgJJY5mRPVZ5ZKlSDkl5-VHa_qi5bhfVcQDxM&tn=D8rHGaIJCwa8Og3I&_nc_ht=scontent.fcai20-3.fna&oh=215b57e346986a002f0d71d5ac061d18&oe=6124CD76)

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/240390881_10225646827178758_2588812202487255242_n.jpg?_nc_cat=110&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeH6iYY4sXoq-PtVR8WdrtGfKvLEkz25rGwq8sSTPbmsbI0W9PNZYH6qo6HWtw4s8wU&_nc_ohc=zSrgQ62Q78EAX-Y-u1r&_nc_ht=scontent.fcai20-3.fna&oh=4cf0962b0311512e5bc0348ebc787e5a&oe=61254AF6)

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/240169098_10225646827578768_8014971172963904116_n.jpg?_nc_cat=104&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeGFQKCezKwD92wndqrvfaeDG5H1YPTmNIwbkfVg9OY0jPmwsfygta-ncYjVKHNB7hY&_nc_ohc=D2HaXxam1p4AX9Fajg2&_nc_ht=scontent.fcai20-3.fna&oh=3a60075320e6cdaffa5a4d466a2aa254&oe=612423E1)

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/239987436_10225646827778773_1891502870608261843_n.jpg?_nc_cat=105&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeGnbdwBY8IClJwm7Pt5IMyEgm4s8KvpWKSCbizwq-lYpEBMKxp3jTzE1uCc4A4hPZU&_nc_ohc=dNsRndgpP9cAX_E12j8&_nc_ht=scontent.fcai20-3.fna&oh=f8ef18b98ef5f7e9300616fe67585935&oe=6124D5BA)

![Quality Report](https://scontent.fcai20-3.fna.fbcdn.net/v/t39.30808-6/239939775_10225646828178783_4978182995878171958_n.jpg?_nc_cat=102&ccb=1-5&_nc_sid=730e14&_nc_eui2=AeEbSPFnKNZSQ1lwLKiEAJN7PU9sMEthmiU9T2wwS2GaJaw4KG3CdEHZIDmkt9PEA10&_nc_ohc=LCgTxenAqUQAX8ba8w4&_nc_ht=scontent.fcai20-3.fna&oh=fe03893345cfc2a4f0382c2c40544e91&oe=6125378E)

	
## Read Mapping
Once the sequence reads have been filtered and trimmed, read mapping or alignment is the next step in the bioinformatic pipeline. Read mapping is the process of aligning a set of reads to a reference genome to determine their specific genomic location. Several tools are available for read mapping but BWA-MEM (Galaxy Version 0.7.17.2) was used for our analysis as it is faster, accurate and supports paired-end reads. Read mapping was performed independently for the normal and tumor tissue using thesame parameters except otherwise stated.

Paramaters
- Locally cached human hg19 reference genome, Paired end reads, Forward and reverse trimmed reads (output of trimmomatic), Set read groups (SAM/BAM specification), auto-assign (no), Read Group Identifier (231335 for normal tissue and 231336 for tumor tissue), Read group sample name (normal for normal tissue and tumor for tumor tissue),
- For parameters not listed, default setting was used.

## Mapped reads postprocessing
	
### Mapped reads filtering
The tool used was a BAM tools filter called: ![](https://i.imgur.com/xkurc1C.png) available on Galaxy. It produces newly filtered BAM datasets and only retains reads mapped to the reference successfully and have a minimal mapping quality of 1 and for which the mate read has also been mapped.
The quality of the output data is controlled by a series ofconditions and filters.

The BAM tools filter was run with these parameters:
The BAM datasets we filtered were:
1. The output of Map with BWA (mapped reads in BAM format)![](https://i.imgur.com/AWEyljF.png)
2. The output of Map with BWA (mapped reads in BAM format)![](https://i.imgur.com/VmShNr5.png)

The quality of the output data is controlled by a series of conditions and filters.
The Conditions set were as below:
The first filter involved selecting a BAM property to filter on which was mapQuality+ the filter on read mapping quality (on a phred scale):>=1 
The mapping quality scale quantifies the probability that a read was misplaced.

The second filter involved selecting another BAM property to filter,for which we selected:isMapped(for mapped reads)+Selected mapped reads>Yes

The third filter involved selecting yet another BAM property to filter for which we selected isMateMapped (for paired-end reads with long inserts)+a confirmation to select mapped reads>yes

The last condition set involving opting to set rules for which we selected >No ![](https://i.imgur.com/p9NcqM1.png)
Then we ran the job. This was done for both the normal and tumor tissue data thus resulting in two datasets in the output results

	
## Duplicate Reads Removal
Tool: ![](https://i.imgur.com/OPq6wgU.png)

#### Significance
RmDup is a tool that identifies PCR duplicates by identifying pairs of reads where multiple reads align to the same exact start position in the genome. PCR duplicates arise from multiple PCR products from the same template binding on the flow cell.These are usually removed because they can lead to false positives<br>The read pair with the highest mapping quality score is kept and the other pairs are discarded.<br>It is important to note that this tool does not work for unpaired reads(in paired end mode) or reads that would pair where each maps to different chromosomes.<br>We used filtered reads datasets(BAM file) from the normal and the tumor tissue data- *outputs of Filter BAM datasets on a variety of attributes* .
We run RmDup on the following parameters:![](https://i.imgur.com/b2IeqaC.png)and ![](https://i.imgur.com/3m1NMRc.png)The result was two new datasets in BAM format.The duplicate rate for both sets was well below 10% which is considered good.The tool standard error reflected the results below of unmatched pairs on chr5 and chr12 that otherwise were not included in the output data.<br>
![](https://i.imgur.com/PCxnoWE.png)
	

### Left-align reads around indels

The first Step in this is running the BamLeftAlign tool from the Tools set available on Galaxy. Then we have chosen the source for the reference genome as Locally cached and selected the filtered and dedicated reads datasets from the normal and the tumor tissue data which were the outputs of RmDup. Then we used Human: hg19 aa the genome reference  and set the maximum number of iterations as 5, keeping all other settings as default and finally this will generate two new datasets, that is,one for each of the normal and tumor data.

### Recalibrate read mapping qualities

The next step after Left aligning the reads around indels is  Recalibrating the read mapping qualities.

•RECALIBRATE READ QUALITY SCORES :
The first Step in Recalibrating read mapping qualities is running CalMD tool from Galaxy tool set. Firstly we have selected the left-aligned datasets from the normal and the tumor tissue data; the outputs of BamLeftAlign tool as the input for the BAM file to recalculate. Them we chose the source of reference genome as Use a built in genome as the required hg 19 reference genome was already in built in the Galaxy version we were using. We chose Advanced options as the choice for Additional options and we selected 50 as the Coefficient to cap the mapping quality of poorly mapped reads. And finally this step would produce two new datasets, that is one for each of the normal and tumor data.

### Refilter reads based on mapping quality

Eliminating reads with undefined mapping quality
We ran Filter BAM datasets on a variety of attributes tool using some parameters.          The  recalibrated datasets from the normal and the tumor tissue data which were the outputs of CalMD were selected as the BAM datasets to filter. Then we applied certain conditions as the options , in Filter, we selecte the MapQuality as the BAM property to Filter. Then set the value of less than or equal to 254 (<=254) as the Filter on read mapping quality (phred scale).
		
## Variant Calling and Classification using VarScan Somatic

The purpose of this step is to use the tool "VarScan Somatic" to detect variant alleles in tumor or normal sample pairs, call sample genotypes at variant sites, as well as classify variants into germline, somatic and LOH event variants using solid classical statistics even in the presence of non-pure samples like those obtained from biopsies.

After generating high quality set of mapped read pairs, we ran the "VarScan Somatic" with some parameters using the Human: h19 genome as our reference genome. 

The mapped and filtered CalMD outputs of the normal and tumor tissue datasets were aligned to be read, and the estimated purity content for the normal and tumor samples were set to 1 and 0.5 respectively. 

We then customized the settings for variants calling and set the "Minimum Base Quality" to 28. This was done to increase the base quality of our sequence data without throwing away a significant portion of the data. 

The "Minimum Mapping Quality" was assigned to 1, in order to filter our sequence reads with a mapping quality of at least one as CalMD might have lowered some mapping qualities to zero.

Then, we set the "Settings for Posterior Variant Filtering"as default values while leaving other settings to their default values as wellbefore executing. The result of this was an output file in VCF format.

## Variant annotation and reporting

The next step after refiltering reads based on mapping quality was to add annotation to the variants.
This was done by importing variant annotations datasets from four different sources into galaxy workspace.
Also gene-level annotation files was also imported from Zenodo, Since galaxy workspace has SnpEff functional genomic annotations installed in it database, Homo sapiens: hg19 was accessible. Alternatively we can use SnpEff Download tool to download genome annotation database hg19.

## Functional annotations to the called variants

The next step after importing variant and gene-level annotations was to utilize SnpEff tool to add functional annotations to the variants in comparison with the reference genome (Homo sapiens hg19).
SnpEff is used to annotate and predict the effects of genetic variants on genes and proteins (such as amino acid changes).
To do this, we ran the SnpEff eff Tool with the following parameters:
* Sequence changes (SNPs, MNPs, InDels) : the output of VarScan somatic tool
* Input format: VCF
* Output format: VCF (only if input is VCF)
* Genome source: Locally installed reference genome
* Genome: Homo sapiens: hg19 (or a similarly named option)

## Adding genetic and clinical evidence_based annotation : Creating a GEMINI database for variants

The next step after adding functional annotation to the called variants was to add "genetic and clinical based annotations"
The processes of the step will help observe more information on the variants like: the clinical and genetic aspects, prevalence in the population and the frequency of occurency.
Firstly, the output from function annotation was loaded into the GEMINI database so as to create a Gemini database where further annotations can be effectively carried
out.
The following optional contents of the variant were loaded as well: GERP scores- these are scores gotten from the constrained elements in multiple alignments by
quantifying the substitution deficits, CAAD scores, [N:B-high CAAD and GERP scores were observed in all pathogenicity components], gene tables, sample genotypes and
variant INFO fields.
These loaded variants in the gemini database are then annoatated by GEMINI annotate; this tool adds more explicit information on the output of VarScan that could not be indentified by Gemini load.

## Making variant call statistics accessible

Hence, we used Gemini annotate to extract three values : Somatic Status(SS), Germline p-value (GPV) and Somatic p-value(SPV) from the info generated by VarScan and added them to the Gemini database.
In order to crosscheck if all information extracted by the GEMINI database are in relation to variants observed in the population, we decided to annote by adding more information from the Single Nucleotide Polymorphism Database(dbSNP), also from Cancer Hotspots, links to CIViC database as well as more information from the Cancer Genome Interpreter (CGI)

For dbSNP: the last output from Gemini annotate was annotated with the imported dbSNP using the Gemini annotate tool. This process extrated dbSNP SNP Allele Origin (SAO) and adds it as "rs_ss" column to the existing database.

For Cancerhotspots: the last ouput generated from annotating dbSNP information was then annotated using GEMINI annotate tool, using the imported cancer hotspots as annotation source to extract "q-values" of overlapping cancerhotspots and add them as "hs_qvalue" column to the existing database.

Links to CIVic: the output of the last Gemini annoate for cancerhotspots was annotated using the imported CIViC bed as annotation source. We extracted 4 elements from this source and again added them as a list of "overlapping_civic_urls" to the existing Gemini database.

For the Cancer Genome Interpreter: the last output with the extracted infomation linking to CIViC was further annotated using the tool Gemini annotate with the imported CGI variants as an annotation source. the information extracted was recorded in the Gemini database as "in_cgidb" being used as the column name.


 ## Reporting Selected Subsets of Variants with GEMINI Query

GEMINI query syntax is built on the SQLite dialect of SQL. This query language enables users express different ideas when exploring variant datasets, for this analysis four (4) stepwise GEMINI queries were carried out.

1.	A query to obtain the report of bona fide somatic variants

“GEMINI database”: the fully annotated database created in the last GEMINI annotate step i.e., Cancer Genome Interpreter(CGI)

“Build GEMINI query using”: *Basic variant query constructor*

“Insert Genotype filter expression”: ```gt_alt_freqs.NORMAL <= 0.05 AND gt_alt_freqs.TUMOR >= 0.10```

This genotype filter aims to read only variants that are supported by less than 5% of the normal sample, but more than 10% of the tumor sample reads collectively

 “Additional constraints expressed in SQL syntax”: ```somatic_status = 2```

This somatic status called by VarScan somatic is one of the information stored in the GEMINI database.

By default, the report of this run would be output in tabular format; and a column header is added to it output.

The following columns were selected

* “chrom”
* “start”
* “ref”
* “alt”

* “Additional columns (comma-separated)”: ```gene, aa_change, rs_ids, hs_qvalue, cosmic_ids```
 These columns are gotten from the variants table of the GEMINI database.

2. This second step has the same settings as the above step except for:

* “Additional constraints expressed in SQL syntax”: ```somatic_status = 2 AND somatic_p <= 0.05 AND filter IS NULL```

3. Run GEMINI query with same settings as step two, excepting:

* In “Output format options”

“Additional columns (comma-separated)”: ```type, gt_alt_freqs.TUMOR, gt_alt_freqs.NORMAL, ifnull(nullif(round(max_aaf_all,2),-1.0),0) AS MAF, gene, impact_so, aa_change, ifnull(round(cadd_scaled,2),'.') AS cadd_scaled, round(gerp_bp_score,2) AS gerp_bp, ifnull(round(gerp_element_pval,2),'.') AS gerp_element_pval, ifnull(round(hs_qvalue,2), '.') AS hs_qvalue, in_omim, ifnull(clinvar_sig,'.') AS clinvar_sig, ifnull(clinvar_disease_name,'.') AS clinvar_disease_name, ifnull(rs_ids,'.') AS dbsnp_ids, rs_ss, ifnull(cosmic_ids,'.') AS cosmic_ids, ifnull(overlapping_civic_url,'.') AS overlapping_civic_url, in_cgidb```

## Generating Reports of Genes Affected by Variants

In this step, gene-centred report is generated based on the same somatic variants we selected above.
As in the previous step we run GEMINI query but in advanced mode

•	“Build GEMINI query using”: *Advanced query constructor*

•	“The query to be issued to the database”: ```SELECT v.gene, v.chrom, g.synonym, g.hgnc_id, g.entrez_id, g.rvis_pct, v.clinvar_gene_phenotype FROM variants v, gene_detailed g WHERE v.chrom = g.chrom AND v.gene = g.gene AND v.somatic_status = 2 AND v.somatic_p <= 0.05 AND v.filter IS NULL GROUP BY g.gene```

However the “Genotype filter expression”: ```gt_alt_freqs.NORMAL <= 0.05 AND gt_alt_freqs.TUMOR >= 0.10``` remains the same








## Adding additional Annotation to the Gene-Centered Report
The aim of including extra annotations to the GEMINI-generated gene report (that is, the output of the last GEMINI query) is to make interpreting the final output easier. While GEMINI-annotate allowed us to add specific columns to the table of the database we created, it does not allow us to include additional annotations into the tabular gene report.

By simply using the Join two files tools on Galaxy, this task was  achieved. After which, irrelevant columns were removed by specifying the columns that are needed. Three step wise process were involved here: One, we pulled the annotations found in Uniprot cancer genes dataset; second, we used the output of the last Join operation, annotated the newly formed gene-centered report with the CGI biomarkers datasets; and three, we used the output of the second Join operation, add the Gene Summaries dataset. Lastly, we ran Column arrange by header name to rearrange the fully-annotated gene-centered report and eliminate unspecified columns.

The last output of the Join operation was selected in the “file to arrange” section. The columns to be specified by name are: gene, chrom, synonym, hgnc_id, entrez_id, rvis_pct, is_TS, in_cgi_biomarkers, clinvar_gene_phenotype, gene_civic_url, and description. The result gotten was a [tabular gene report](https://github.com/Fredrick-Kakembo/Somatic-and-Germline-variant-Identification-from-Tumor-and-normal-Sample-Pairs/blob/main/Galaxy54-%5BColumn_arrange_on_data_53%5D%20(1).tabular), which was easy to understand and interpret.

## Conclusion
The Tumor/Normal data analysis workflow, consisting of Alignment and VarScan Somatic variant calling, demonstrates excellent performance for the detection of somatic variants. <br>
Somatic variant calling not only calls variants but also distinguishes Cancer-specific variants (Somatic mutations) found only in tumor tissues from germline mutations that are shared by tumor and healthy tissue, and loss-of-heterozygosity events; that is, the absence of one of two alleles found at a biallelic site of healthy tissues, in tumor tissues. This makes it a more optimal approach for applications requiring high precision such as novel mutation detection and mutation signature analysis. <br>
The interpretation of any list of variants (somatic, germline or LOH) almost always depends crucially on rich genetic and cancer-specific variants and gene annotations which can be inferred from the analysis.



---
##  Contribution of team members according to the environment used<a name="contributor">:</a>

1. Galaxy Workflow:
Each person that chose to work on Galaxy ran the tutorial from start to end and Links to their workflow added besides their names. After each member had to document a chosen section within the Galaxy workflow section above.    
In additon, We all came up with one comprehensive workflow for the entire team that can be [Found Here](https://usegalaxy.eu/u/rachael-eo/w/workflow-constructed-from-history-team-genomicstwoa)

- @Rachael - Documented genetic and clinical evidence-based annotations under Annotation with Gemini. [Link to her whole galaxy workflow](https://usegalaxy.eu/u/rachael-eo/w/workflow-constructed-from-history-genomics-twoarachael-1)
- @Mercy - Documented Variant Calling and Classification Using VarScan Somatic [Link to her whole galaxy workflow](https://usegalaxy.eu/u/mercyoni/w/workflow-constructed-from-history-genomic-two-a-mercy)
- @Neesah - Written about Variant Functional annotation using snpEff eff [Link to her complete galaxy workflow](https://usegalaxy.eu/u/nerdy_neesah1./w/workflow-constructed-for-identification-of-somatic-and-germline-variants-from-normal-and-tumor-sample-pairs-tutorial)
- @Orinda - Documented about Filtering of Mapped reads and Duplicate reads removal. [Link to her galaxy workflow](https://usegalaxy.eu/workflow/display_by_id?id=52354430d02f285c)
- @Heshica - Documented about Postprocessing of Mapped Read ie Left-align reads around indels , Recalibrate read mapping qualities and Refilter reads based on mapping quality)[Link to her full Galaxy Workflow](https://usegalaxy.eu/u/heshica_battina_chowdary/w/normal-and-tumor-analysisheshica-genomics-2a)
- @VioletNwoke - Documented about Read mapping using BWA-MEM [Link to galaxy workflow](https://usegalaxy.eu/u/violet/w/workflow-constructed-from-history-hackbiogenomicstwoaviolet-4)
- @AmaraA - Documented about 
- @Amarachukwu -Reporting Selected Subsets of Variants and Generating Reports of Genes Affected by Variants(GEMINI Query) [Link to Galaxy workflow](https://usegalaxy.eu/u/amara_chike/w/somatic-variant-tutorial-genomics-2-a-1) 
- @Mallika [Link to Galaxy Workflow](https://usegalaxy.eu/u/mallika_g/w/variant-analysis-mallika)
- @Olamide - Read Trimming and Filtering [Link to Galaxy Workflow](https://usegalaxy.eu/u/olamide21/w/identification-of-somatic-and-germline-variants-from-tumor-and-normal-sample-pairs) 
- @NadaaHussienn - Quality Control and Check [Link to Galaxy Workflow](https://usegalaxy.eu/u/nadahussien/w/workflow-constructed-from-history-identification-of-somatic-and-germline-variants-from-tumor-and-normal-sample-pairs-3)
- @Christabel- Conclusion [link to galaxy workflow](https://usegalaxy.eu/u/christabelmn1/w/somatic-and-germline-variants-and-gene-mutation-2)
- @Marvellous - Adding additional Annotation to the Gene-Centered Report [Galaxy Workflow](https://usegalaxy.eu/u/marvellous_oyebanjo/w/workflow-constructed-from-history-identification-of-somatic-and-germline-variants-from-tumor-and-normal-sample-pairs3) 
- @juwon - Introduction


2. Linux Workflow
- @Praise
- @Fredrick
- @RuthMoraa
- @Kauthar - Preprocessing(pre/post trim qc) and read trimming.
- @Gladys
- @Nanje
