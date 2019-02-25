# Sicelore

Acronyme for [Si]ngle [Ce]ll [Lo]ng [Re]ad is a suite of tools dedicated 
to the bioinformatics processing, analysis and exploration of highly 
multiplexed single cell droplets-based approach experiments sequenced 
with Oxford Nanopore Technology long reads. Typically starting with a 
possorted_genome.bam file from a 10xGenomics single cell cellranger pipeline, 
the workflow integrate several sequential steps for standard analysis and 
processing of sequencing runs of the same unfragmented library on Nanopore device.

[![GitHub license]()]((https://github.com/hyeshik/poreplex/blob/master/LICENSE.txt))
[![Twitter Follow](https://img.shields.io/twitter/follow/kevinlebrigand.svg?style=social&logo=twitter)](https://twitter.com/kevinlebrigand)

## Installation

*sicelore* requires Java 8, <a href="https://github.com/lh3/minimap2">minimap2</a>, <a href="http://mccarrolllab.com/download/1276/">Drop-seq tools v1.13</a>, <a href="https://github.com/isovic/racon">racon</a>

## Features

* [10xGenomics possorted_genome_bam.bam preprocessing](#preprocessing)
* [Long reads scanning](#scanning)
* [minimap2 mapping of Nanopore reads](#minimap2-mapping)
* [Long reads SAM records Gene name tagging (GE) with dropseq.jar](#gene-name-tagging)
* [Long reads SAM records read sequence tagging (US)](#TagReadWithSequence)
* [Short reads 10xGenomics attributes cell barocde (BC) and UMI (U8) association to Long reads](#10x-attributes-association)
* [Long reads Isoforms Expression assessment](#IsoformExpressionMatrix)
* [Molecule Consensus sequence computation for long read error correction](#MoleculeConsensus)
* [Mapping of molecule sequence to obtain molecular BAM file](#molecule-mapping)
* [SNPs calling, fusion transcript detection]()


## preprocessing
Parsing of 10xGenomics possorted_genome_bam.bam file for production of parsedForNanopore.obj java object required during illumina barcodes (BC/UMI) transfert to Nanopore reads.

```bash

```
## scanning
Long read scan operation to filter out low quality reads in which we won't be abble to find the Illlumina barcodes. The reads where we detect a polyA tail are kept and eventually reversed/complemetned to fit for structure (TSO)-(cDNA)-(polyA)-(BC)-(UMI)-(ADAPTOR)
so the next steps could considered as stranded.

```bash

```

## minimap2 mapping
Reads are first splitted in 24 chunks for time optimization and distribution accross the calcul custer.
Parallel Nanopore long read genome mapping  is then performed using minimap2.


```bash
# fastq splitting
fastp -i PROMxxxxxx.fastq -Q -A --thread 20 --split_prefix_digits=4 --out1=sub.fastq --split=24

# then parallel mapping
minimap2 -a -x splice -t 20 -N 100 $BUILD.mmi 0001.sub.fastq > 0001.sub.sam
"/bin/awk '{ if($3 !="*") print $0 }' 0001.sub.sam > 0001.sub.match.sam
samtools view -Sb 0001.sub.match.sam -o 0001.sub.unsorted.bam
samtools sort 0001.sub.unsorted.bam -o 0001.sub.bam
samtools index 0001.sub.bam
```

## gene name tagging
Assignment of GE tag to Nanopore reads using dropseq.jar from McCarrol lab.
Optimization for taking into account long read particularities compared to short read for wich this pipeline has been implemented.

```bash
cd ~/Drop-seq_tools-1.12/jar/
java -jar -Xmx12g dropseq.jar TagReadWithGeneExon I=0001.sub.bam O=0001.sub.GE.bam ANNOTATIONS_FILE=~/cellranger_references/refdata-cellranger-mm10-1.2.0/genes/genes.gtf TAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index minimap.GE.bam
```

## TagReadWithSequence
Sam flag US assigment of Nanopore raw read sequence

```bash
java -jar -Xmx12g sicelor.jar TagReadWithSequence I=0001.sub.GE.bam O=0001.sub.GEUS.bam FASTQ=nanopore.fastq
samtools index minimap.GEUS.bam
```

## 10x attributes association
Assignation of Illumina cell barcode and UMI to Nanopore mapped reads

```bash
java -jar -Xmx25g IlluminaOxfordMergerNew.jar -i 0001.sub.GEUS.bam -o 0001.sub.GEUS10xAttributes.bam -k parsedForNanopore.obj -p CTTCCGATCT -a 140 -s GTACATGG  --maxUMIfalseMatchPercent 6 --maxBCfalseMatchPercent 5 -l minimap.GEUS10xAttributes.log
```

## IsoformExpressionMatrix

```bash
java -jar -Xmx44g Sicelor-1.0-SNAPSHOT-jar-with-dependencies.jar IsoformExpressionMatrix I=0001.sub.GEUS10xAttributes.umifound.bam REFFLAT=refFlat_gencode.vM18.txt CSV=10xgenomics.barcodes.csv MATRIX=MatrixIsoforms.txt DELTA=10 METRICS=MetricsIsoforms.txt
```

## MoleculeConsensus

```bash
java -jar -Xmx22g sicelor.jar MoleculeConsensus
```

## molecule mapping

```bash
minimap2 -ax splice -t ... $BUILD.mmi molecule_consensus.fasta > molecule.sam
```
