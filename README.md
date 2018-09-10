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
* [minimap2 mapping of Nanopore reads](#minimap2-mapping)
* [Long reads SAM records Gene name tagging (GE) with dropseq.jar](#gene-name-tagging)
* [Long reads SAM records read sequence tagging (US)](#TagReadWithSequence)
* [Short reads 10xGenomics attributes cell barocde (BC) and UMI (U8) association to Long reads](#10x-attributes-association)
* [Long reads Isoforms Expression assessment](#IsoformExpressionMatrix)
* [Molecule Consensus sequence computation for long read error correction](#MoleculeConsensus)
* [Mapping of molecule sequence to obtain molecular BAM file](#molecule-mapping)
* [SNPs calling, fusion transcript detection]()


## preprocessing
Parsing of 10xGenomics possorted_genome_bam.bam file for production of parsedForNanopore.obj java object required for illumina cell barcode and UMI transfert to Nanopore reads

```bash

```

## minimap2 mapping
Nanopore reads genome mapping step

```bash
minimap2 -ax splice -t ... $BUILD.mmi nanopore_reads.fastq > minimap.sam
"/bin/awk '{ if($3 !="*") print $0 }' minimap.sam > minimap.match.sam
samtools view -Sb minimap.match.sam -o minimap.unsorted.bam
samtools sort minimap.unsorted.bam -o minimap.bam
samtools index minimap.bam
```

## gene name tagging
Assignment of GE tag to Nanopore reads using dropseq.jar from Maccarrol lab

```bash
cd ~/Drop-seq_tools-1.12/jar/
java -jar -Xmx12g dropseq.jar TagReadWithGeneExon I=minimap.bam O=minimap.GE.bam ANNOTATIONS_FILE=~/cellranger_references/refdata-cellranger-mm10-1.2.0/genes/genes.gtf TAG=GE ALLOW_MULTI_GENE_READS=TRUE USE_STRAND_INFO=FALSE
samtools index minimap.GE.bam
```

## TagReadWithSequence
Sam flag US assigment of Nanopore raw read sequence

```bash
java -jar -Xmx12g sicelor.jar TagReadWithSequence I=minimap.GE.bam O=minimap.GEUS.bam FASTQ=nanopore.fastq
samtools index minimap.GEUS.bam
```

## 10x attributes association
Assignation of Illumina cell barcode and UMI to Nanopore mapped reads

```bash
java -jar -Xmx25g IlluminaOxfordMergerNew.jar -i minimap.GEUS.bam -o minimap.GEUS10xAttributes.bam -k parsedForNanopore.obj -p CTTCCGATCT -a 140 -s GTACATGG  --maxUMIfalseMatchPercent 6 --maxBCfalseMatchPercent 5 -l minimap.GEUS10xAttributes.log
```

## IsoformExpressionMatrix

```bash
java -jar -Xmx44g Sicelor-1.0-SNAPSHOT-jar-with-dependencies.jar IsoformExpressionMatrix I=minimap.GEUS10xAttributes.umifound.bam REFFLAT=refFlat_gencode.vM18.txt CSV=10xgenomics.barcodes.csv MATRIX=MatrixIsoforms.txt DELTA=10 METRICS=MetricsIsoforms.txt
```

## MoleculeConsensus

```bash
java -jar -Xmx22g sicelor.jar MoleculeConsensus
```

## molecule mapping

```bash
minimap2 -ax splice -t ... $BUILD.mmi molecule_consensus.fasta > molecule.sam
```
