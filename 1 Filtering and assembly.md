# 1. FILTER RAW READS BY QUALITY
---
## Visualising raw reads

*FastQC* is an extremely popular tool for checking your sequencing libraries, as the visual interface makes it easy to identify issues like:

1. Adapter/barcode sequences
2. Low quality regions of sequence
3. Quality drop-off towards the end of read-pair sequence
   
### Counting number of reads

Before using our first software, a recommended step is to see how many reads we have. It is not weird that we encounter errors when downloading our fasta files. As we are working with Illumina paired-end reads, both files (forward and reverse) should have same number of reads. Most sequencing reads are in FASTQ format. Each sequence in FASTQ format in represented in four consecutive lines:

```bash 
@HWUSI-EAS1752R:23:FC62KHPAAXX:6:3:3542:1008 1:N:0:GCCAAT
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
````
- First line, is the **sequence name** beginning with "@". The name of the sequence usually contains information of the instrument, the number of run, the lane, etc.
- Second line is the **sequence itself**.
- Third line is for **comments** (usually empty).
- Last line is the **quality** of the sequence code using Phred scale (currently most sequencing platforms used Phred+33 scale).

So, the easy way of counting the number of reads in a file is using wc linux command (word count):

```bash
wc -l A_1_R1.fastq
wc -l A_1_R2.fastq
````
-l option is applied to count the number of lines. However, we have to divide by 4 to get the number of reads. To avoid this, we can apply some piping:
```bash
wc -l A_1_R1.fastq | awk '{print $1/4}' #37560016
wc -l A_1_R2.fastq | awk '{print $1/4}' #37560016
````

### Check sequence quality with *FastQC*

The first step is to install *FastQC*
```bash
conda create env -n fastqc
conda activate fastqc
conda install -c bioconda fastqc -y
````
Now, run the program to check our reads:
```bash
mkdir A_quality_reads
fastqc A_1_R1.fastq.gz -o A_quality_reads/
fastqc A_1_R2.fastq.gz -o A_quality_reads/
````
>NOTE: In order to repeat this process with all the samples (12) we can run the command manually or use a loop. We will use them later on.
```bash 
for item in [LIST]
do 
    [COMMANDS]
done
````
Open the html files to get information about the number of sequences, the length distribution, the %GCs content, the average quality, etc.
#### R1 Quality plot
![r1_quality](https://user-images.githubusercontent.com/13121779/162767027-92e2adeb-bec8-4571-8257-3196cd7de944.png) 
#### R2 Quality plot
![r2_quality](https://user-images.githubusercontent.com/13121779/162767078-d14e19d6-40a9-498a-828d-f3b65ebad31e.png)

## Read trimming and adapter removal with *trimmomatic*

There are multitude of programs which can be used to quality trim sequence data and remove adapter sequence. We will use *trimmomatic*.

*Trimmomatic* installation (inside fastqc enviroment)
```bash
conda install -c bioconda trimmomatic -y
````
Running *trimmomatic* We will create a loop file saved as run_trimmomatic.sh. Don't forget to make it executable with chmod +x.
```bash
#!/bin/bash

# Loop from 1 to 12
for i in {1..12}
do
    # Construct file names
    R1="A_${i}_R1.fastq.gz"
    R2="A_${i}_R2.fastq.gz"
    R1_paired="A_${i}_R1_qf_paired.fastq.gz"
    R1_unpaired="A_${i}_R1_qf_unpaired.fastq.gz"
    R2_paired="A_${i}_R2_qf_paired.fastq.gz"
    R2_unpaired="A_${i}_R2_qf_unpaired.fastq.gz"

    # Run Trimmomatic
    trimmomatic PE -phred33 $R1 $R2 $R1_paired $R1_unpaired $R2_paired $R2_unpaired SLIDINGWINDOW:4:30 MINLEN:80
done
conda deactivate
````

Let's take a look to some parameters:

- SLIDINGWINDOW: Performs a sliding window trimming approach. It starts scanning at the 5‟ end and clips the read once the average quality within the window falls below a threshold.

- MINLEN: Drop the read if it is below a specified length.

> Now try to check the number of high-quality paired-end reads selected in this pre-processing step and compare their quality profile with that of the input

#### Considerations when working with *trimmomatic*
The basic format for a *trimmomatic* command is 
```bash 
trimmomatic PE <keyword flags> <sequence input> <sequence output> <trimming parameters>
````
The trimming parameters are processed in the order you specify them. This is a deliberate behaviour, but can have some unexpected consequences for new users. For example.
```bash 
trimmomatic PE <keyword flags> <sequence input> <sequence output> SLIDINGWINDOW:4:30 MINLEN:80

trimmomatic PE <keyword flags> <sequence input> <sequence output> MINLEN:80 SLIDINGWINDOW:4:30
````
In the first run, we would not expect any sequence shorter than 80 base pairs to exist in the output files. However, we might encounter them in the second command. This is because in the second command we remove sequences shorter than 80 base pairs, **then** perform quality trimming. If a sequence is trimmed to a length shorter than 80 base pairs **after** trimming, the MINLEN filtering does not execute a second time. In the first instance, we do not perform trimming **before** size selection, so any reads that start longer than 80 base pairs, but are trimmed to under 80 base pairs during quality trimming will be caught in the MINLEN run.

## Removal of human genomes (Bowtie 2)
Human contamination is very frequent, as humans prepare and manipulate the samples before being sequenced. In some cases, the genome of other organisms, such as mouse or monkey, should be used instead of human, depending on the experimental set. To perform short read alignment we are going to use Bowtie2. There are other popular aligners such as BWA, but Bowtie2 parameters are easy to understand and modify.

Before preforming the alignment, Bowtie2 requires the preparation of an index containing the reference sequences to align against. 

```bash
conda create -n bowtie2 -c bioconda bowtie2 -y
```
**Decontamination human reads**

The file containing the human reads can be downloaded from [Zenodo](https://zenodo.org/records/1208052). Then, we build the index:

```bash
wget https://zenodo.org/records/1208052/files/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz?download=1
bowtie2-build hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz human_cds
````
**Aligning reads against human reads**

```bash
for i in {1..12}; do
  sample="A_${i}"
  bowtie2 -x human_cds \
    -1 ${sample}_R1_qf_paired.fastq \
    -2 ${sample}_R2_qf_paired.fastq \
    --un-conc ${sample}_qf_paired_nohuman_R%.fastq \
    -S ${sample}_tmp.sam
done
```
Now we can count clean reads

```bash
wc -l *_qf_paired_nohuman* | awk '{print $1/4}'
````


# 2. Assembling reads -> contigs
---
We are going to use [SPAdes](https://github.com/ablab/spades) one of the most popular de Bruijn graph de novo assemblers. 

This assembler has specific assembly protocols for different experimental set (--sc, single cell data; --meta, metagenomics data; --isolate, a single genome sequencing or --rna, for RNA-seq) and sequencing platforms (Illumina, IonTorrent, Nanopore or PacBio). The program does not only make de novo assembly but also run an error correction tool for Illumina sequences. If we run the program with the flag --careful, SPAdes will also run MismatchCorrector, a post processing tool which uses BWA aligner. However, this flag is only recommended for assembly of small genomes and is not compatible with some other options such as single cell (--sc) or metagenomics (--meta). So we won't be running that flag today.

Another important option is -k flag with which we can define the kmer length to use in the assembly. If we do not provide a list of kmers, the program will select automatically three of them based on average reads length and perfoms a combined assembly. Nowadays, with this modern assembler it is better to allow the program to choose the best kmers, but at the beginning of the genomic era researches had to choose them by making multiple assemblies. As a general rule, kmers have to be smaller than average read length and must be odd and less than 128 (althougth other assemblers can use longer kmer sizes).

Awkwardly, while SPAdes accepts multiple input libreries in a single assembly, this behaviour does not work with the --meta flag enabled, which would be the desire flag to use in our samples. **This is why we are going to assemble each sample sepparetly with the --careful flag**. We also could use --isolate flag, but literature shows that --careful flag gives better results.

If we wanted to use --meta flag, we would need to concatenate R1 files and R2 files separetly.

But first, we need to install SPAdes:
```bash
conda create env -n spades
conda activate spades
conda install -c bioconda spades -y
#############################################################################################
#                                                                                           #
#   Note: SPAdes installed through bioconda on MacOS may be somewhat slower than the SPAdes #
#   binaries distributed by the authors at                                                  #
#                                                                                           #
#   http://cab.spbu.ru/files/release3.15.2/SPAdes-3.15.2-Darwin.tar.gz
#                                                                                           #
#   due to unavailability of parallel libstdc++ for the Clang compiler used by bioconda on  #
#   MacOS; see https://github.com/ablab/spades/issues/194#issuecomment-523175204            #
#                                                                                           #
#############################################################################################
````
#### Run SPAdes (.sh)
```bash
#!/bin/bash

# List of indices you want to iterate over
for i in {1..12} # Change this range as needed
do
    echo "Running SPAdes for sample $i"
    
    # Construct the input file names
    R1="A_${i}_R1_qf_paired.fastq.gz"
    R2="A_${i}_R2_qf_paired.fastq.gz"
    
    # Output directory
    output_dir="Assembly_careful/A_${i}_careful"
    
    # Run SPAdes
    spades.py --careful -t 2 -1 $R1 -2 $R2 -o $output_dir
done
````
Spades is also really useful when working with isolate genomes. We can add the flag --isolate and then compare what assembly turned "better" :
```bash 
mkdir Assembly_isolate
# List of indices you want to iterate over
for i in {1..12} # Change this range as needed
do
    echo "Running SPAdes for sample $i"
    
    # Construct the input file names
    R1="A_${i}_R1_qf_paired.fastq.gz"
    R2="A_${i}_R2_qf_paired.fastq.gz"
    
    # Output directory
    output_dir="/Assembly_isolate/A_${i}_isolate"
    
    # Run SPAdes
    spades.py --isolate -t 2 -1 $R1 -2 $R2 -o $output_dir
done
```

### Check quality of assemblies:

We can evaluate the assemlies using longMeta-summary (Varliero et al., 2021). [See installation](https://github.com/gvMicroarctic/LongMeta). But before that, it is recommended to see how many contigs or scaffolds we have. We can count sequence numbers using grep.
It is recommended to have all assemblies (fasta files and scaffolds) in the same folder. To reduce redundancy and save a bit of disk space, we are going to use symbolic links:
```bash
mkdir Quast # This folder will be useful for later
ln -rs ./A_1_careful-contigs.fasta ./Quast/A_1_careful-contigs.fasta
ln -rs ./A_1_careful-scaffolds.fasta ./Quast/A_1_careful-scaffolds.fasta
ln -rs ./A_1_isolate-contigs.fasta ./Quast/A_1_isolate-contigs.fasta
ln -rs ./A_1_isolate-scaffolds.fasta ./Quast/A_1_isolate-scaffolds.fasta
# Repeat with all assemblies or use a for loop ;)
```
We can take a look to the number of contigs/scaffold we have obtained from the assemblies:
```bash
grep -c '>' ./Quast/*.fasta
````
### longMeta-summary (Varliero et al., 2021):

```bash
for file in ./Quast/*.fasta; do longMeta-summary --assembly-input "$file" >> assemblies_stats.txt; done
```
The downside of longMeta-summary is that it doesn't tell us the % of contigs above a certain length. **Bbmap's stats.sh** can do that. So if we wanted that metric, we can run stats.sh on our assemblies:
```bash
for file in /Quast/*.fasta; do bash ./stats.sh in= "$file" >> assemblies_results_bbmap.txt; done
```
Ok, we have our stats on every contig and scaffold file with both spades methods but, how do we chose one assembly? There is actually no right answer. There are some general things we can look at, like N50 or largest contig, or fraction of reads that successfully recruit to our assembly. But for one, these don't really have any solid context to know if they're "good" or not unsless we're comparing multiple assemblies of the same data; and two, we can have metagenomic assemblies with "worse" overall summary statistics, but that might enable us to recover more high-quality bins than an assembly with "better" summary statistics. Having a reference genome, can make things a lot easier. 

*From now on, sample names will be different, since we want to work with metagenomes from a bunch of Tolypothrix isolates*

## Comparing assemblies - QUAST
A better and more complete way of making assebmblies comparision is by using dedicated tools such as [QUAST](https://www.google.com/url?sa=t&source=web&rct=j&opi=89978449&url=https://github.com/ablab/quast&ved=2ahUKEwiX0-q-vPeKAxWZSfEDHXQkG48QFnoECAwQAQ&usg=AOvVaw3ngoxAk8Wa13hhRzdGn08y).
### Installing QUAST
We must create a new enviroment
```bash
# conda deactivate first
conda create -n quast -c bioconda quast -y
conda activate quast
````
*Note: We don't have a reference genome of our isolates. Either way, I am going to explain how it would be done in case we had one*

We can provide QUAST with all of our assemblies, a fasta file of our reference genome, and a .gff (**g**eneral **f**eature **f**ormat) file of our reference genome. I downloaded two reference files for our *Tolypothrix sp.* [PCC 7712](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_025860405.1/) from NCBI. Now we can run QUAST:
```bash
cd Quast
quast -r Tolyp_ref.fna -g Tolyp_ref.gff *.fasta
conda deactivate
````
We can see the report.html to see how well we did. Genome fraction (%) values are low. This means our assemblies have other bacterial genomes in it. Since we know this, we will proceed with careful spades assemblies for both samples. 
We could trim our samples even more and reduce our coverage. But for now, we are doing good. 
