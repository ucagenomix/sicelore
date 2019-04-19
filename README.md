# SiCeLoRe

Acronyme for [Si]ngle [Ce]ll [Lo]ng [Re]ad is a suite of tools dedicated to the bioinformatics processing, analysis and exploration of highly 
multiplexed single cell droplets-based approaches experiments sequenced with long reads. Typically starting with a possorted_genome.bam file from 
a 10xGenomics single cell cellranger pipeline, the workflow integrate several sequential steps for standard analysis and processing of sequencing 
runs of the same unfragmented library on Nanopore/PacBio device.

[![GitHub license]()]((https://github.com/ucagenomix/Sicelore/blob/master/LICENCE.md))
[![Twitter Follow](https://img.shields.io/twitter/follow/kevinlebrigand.svg?style=social&logo=twitter)](https://twitter.com/kevinlebrigand)

## Installation

*sicelore* requires Java 8, <a href="https://github.com/lh3/minimap2">minimap2</a>, <a href="http://mccarrolllab.com/download/1276/">Drop-seq tools v1.13</a>, <a href="https://github.com/isovic/racon">racon</a>

## Features

* [10xGenomics bam file preprocessing](#preprocessing)
* [Long reads scanning and filtering](#scanning)
* [minimap2 mapping of long reads](#minimap2-mapping)
* [Gene name (GE) SAMrecords tagging (dropseq.jar)](#gene-name-tagging)
* [Read sequence (US) and quality values (UQ) SAMrecords tagging](#TagReadWithSequence)
* [Transfert of 10xGenomics attributes cell barocde (BC) and UMI (U8) to long reads](#10x-attributes-association)
* [Isoforms expression quantification](#IsoformMatrix)
* [Consensus molecule sequence calculation](#ComputeConsensus)
* [Molecule sequence genome remapping, production of a molecular BAM file](#molecule-mapping)
* [SNPs calling, fusion transcript detection]()


## preprocessing
Short reads 10xGenomics possorted_genome_bam.bam file parsing and production of the parsedForNanopore.obj java object 
required during illumina cellBC/UMI (BC/U8) transfert to long reads.

```bash

```
## scanning
Long read scanning operation to filter out low quality reads in which we won't be abble to find the Illlumina barcodes. 
Reads where we detect a polyA tail are kept and eventually reversed/complemented to fit for structure (TSO)-(cDNA)-(polyA)-(BC)-(U8)-(ADAPTOR)
so that the next steps can be considered as stranded.

```bash

```

## minimap2 mapping
For time calculation purpose the next steps can be parrallelized and distribute independently accross compuet node of a calcul cluster.

```bash
# fastq splitting in 24 chunks for time optimization
fastp -i PROMxxxxxx.fastq -Q -A --thread 20 --split_prefix_digits=4 --out1=sub.fastq --split=24

# parallel minimap2 mapping using 20 threads
minimap2 -a -x splice -t 20 -N 100 $BUILD.mmi 0001.sub.fastq > 0001.sub.sam
"/bin/awk '{ if($3 !="*") print $0 }' 0001.sub.sam > 0001.sub.match.sam
samtools view -Sb 0001.sub.match.sam -o 0001.sub.unsorted.bam
samtools sort 0001.sub.unsorted.bam -o 0001.sub.bam
samtools index 0001.sub.bam
```

## gene name tagging
SAMrecords gene name (GE) tagging using dropseq.jar from McCarrol lab.

```bash
cd ~/Drop-seq_tools-1.12/jar/
java -jar -Xmx12g dropseq.jar TagReadWithGeneExon I=0001.sub.bam O=0001.sub.GE.bam ANNOTATIONS_FILE=~/cellranger_references/refdata-cellranger-mm10-1.2.0/genes/genes.gtf TAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index 0001.sub.GE.bam
```

## TagReadWithSequence
SAMrecords read sequence (US) and read quality (UQ) tagging that will be used for BC/U8 assignment and molecule consensus calculation.

```bash
java -jar -Xmx12g sicelor.jar AddBamReadSequenceTag I=0001.sub.GE.bam O=0001.sub.GEUS.bam FASTQ=nanopore.fastq
samtools index 0001.sub.GEUS.bam
```

## 10x attributes association
Transfert of 10xGenomics attributes cellBC (BC) and UMI (UB) to Nanopore long reads

```bash
java -jar -Xmx40g IlluminaOxfordMergerNew.jar -i 0001.sub.GEUS.bam -o 0001.sub.GEUS10xAttributes.bam -k parsedForNanopore.obj -p CTTCCGATCT -a 140 -s GTACATGG  --maxUMIfalseMatchPercent 6 --maxBCfalseMatchPercent 5 -l minimap.GEUS10xAttributes.log
```

## IsoformMatrix
The IsoformMatrix pipeline takes as input: 
(i) CSV: .csv file listing, one per line (no '-1'), the cell barcodes that need to be quantified (cf. brain951.barcodes.tsv in barcodes dir); 
(ii) REFFLAT: .txt file describing the gene model for the organism concerned (cf. refFlat_gencode.vM18.txt in refFlat dir);
(iii) INPUT: global .bam file containing IG/BC/U8 short reads transfered barcodes. SAMrecords whithout those 3 required fields are discards 
from further analysis. SAMrecords starting or ending with Hard or Soft clipping sequence longer than 150 nt are also discarded in order 
to eliminate the chimeric reads that can arise during Nanopore library preparation. GeneId and TranscriptId(isoform) annotation is then 
set for each molecule detected by comparing all the reads of the molecule versus the lists of the isoforms belonging to the genes 
corresponding to the SAM records IG Tags. This attribution requires that the SAMrecord recapitulates the exon-2-exon junctions layout 
of the model transcript with a latitude given on the mapping precision equal to DELTA (default=5). When more than one read is
available for the molecule, all are processed and the TranscriptId is set to the most often identified transcripts.
A SOFT (=true) mode can be appliyed to IsoformMatrix pipeline to quantify truncated molecules not spanning the entire set 
of exon-2-exon junction of the transcript model due to protocol issue and internal polyA-like priming which can arise in 
single cell leading to truncated molecules profiling. In SOFT mode, you need to provide a REFFLAT gene model file listing 
sub-sampling of the entire set of the transcripts exons, in order to quantify for instance the inclusion/exclusion of one specific exon 
that we know is specific of an isoform. At the end of the pipeline, several files (names starting with PREFIX paramter) are produced including 
Isoforms and Genes metrics and matrices used for subsequent statistical analysis.

```bash
# chunks bam file aggregation
java -jar -Xmx44g ~/picard-tools-1.119/MergeSamFiles.jar INPUT=0001.sub.GEUS10xAttributes.bam INPUT=0002.sub.GEUS10xAttributes.bam INPUT=...  ASSUME_SORTED=true USE_THREADING=true TMP_DIR=/scratch/tmp/ MAX_RECORDS_IN_RAM=100000000 OUTPUT=GEUS10xAttributes.umifound.bam VALIDATION_STRINGENCY=SILENT
samtools index GEUS10xAttributes.umifound.bam

# IsoformMatrix pipeline
java -jar -Xmx44g sicelor.jar IsoformMatrix I=GEUS10xAttributes.umifound.bam REFFLAT=refFlat_gencode.vM18.txt CSV=10xgenomics.barcodes.csv DELTA=5 SOFT=false PREFIX=sicelore
```

## ComputeConsensus
ComputeConsensus pipeline allows to compute the consensus sequence for each of the molecule identified in the INPUT .bam file. 
The first steps are very similar to the IsoformMatrix pipeline except that the cDNA sequence is loaded into LongreadRecord Objects 
so that it can be used for consensus sequence computation. The default number of threads is set to 20 but can be change using the “T” argument. 
Briefly, each molecule are processed as follow depending on the number of reads the molecule behave: (i)case of a 1-read molecule, the consensus sequence is set to the 
read sequence; (ii) case of a 2-reads molecule, the consensus sequence is set as the sequence of the best quality read according 
to the minimal “dv” minimap2 SAMrecord tag value; (iii) Case of a multi-reads molecule (i.e. > 2), a consensus sequence is set 
using Biojava API multiple alignment routines using a maximum of 10 reads selected as the ones minimizing the “dv” value. 
The consensus sequence is then “racon” polished using the whole set of reads for the molecule. All the computed molecules 
consensus sequence are finally grouped as a “OUTPUT” Fasta formatted file so that they can be mapped-back to the genome using minimap2 or gmap 
mapper for subsequent analysis such as SNPs or editing events calling. For time calculation optimization, this step should be parrallelized, 
for instance on a per chromosome basis.

```bash

# splitting by chromosomes
@jobs = ('1','2','X','MT','3','4',....);
for($i=0; $i<@jobs; $i++){
    samtools view -Sb GEUS10xAttributes.umifound.bam ".$jobs[$i]." -o GEUS10xAttributes.umifound.chr".$jobs[$i].".bam
    samtools index GEUS10xAttributes.umifound.chr".$jobs[$i].".bam
}

# ComputeConsensus pipeline for chr1
java -jar -Xmx22g sicelor.jar ComputeConsensus I=GEUS10xAttributes.umifound.chr1.bam O=molecules.fa REFFLAT=refFlat_gencode.vM18.txt T=20 DELTA=5 SOFT=false TMPDIR=/tmp/ RACONPATH=/share/apps/local/racon/bin/ MINIMAP2PATH=/share/apps/local/minimap2/
```

## molecule mapping
Molecule consensus sequence can be then mapped to the genome to produce a molecule .bam file from which subsequent analysis can be performed
such as editing or SNPs calling. Gene (IG), cellBC (BC) and UMIs (U8) tags can be added to records using AddBamMoleculeTags pipeline. 
This file can then be splitted by clusters using SplitBamPerCluster pipeline gicing a .csv file of 'cellBC,clusterName', one cell per line
to produce per cluster bam files. 

```bash
minimap2 -a -x splice -t -N 100 $BUILD.mmi molecules.fa > molecules.sam
"/bin/awk '{ if($3 !="*") print $0 }' molecules.sam > molecules.match.sam
samtools view -Sb molecules.match.sam -o molecules.unsorted.bam
samtools sort molecules.unsorted.bam -o molecules.bam
samtools index molecules.bam

java -jar -Xmx22g sicelor.jar AddBamMoleculeTags I=molecule.bam O=molecule.tags.bam
java -jar -Xmx22g sicelor.jar SplitBamPerCluster I=molecule.tags.bam O=./bams_per_cluster_directory/ CSV=./R/data/951c/illumina/barcodes.tsv CELL_FLAG=BC

```
