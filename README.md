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

*sicelore* requires Java 8, <a href="https://github.com/lh3/minimap2">minimap2</a>, <a href="http://mccarrolllab.com/download/1276/">Drop-seq tools v1.13</a>

## Features

* [10xGenomics possorted_genome_bam.bam preprocessing](#preprocessing)
* [minimap2 mapping of Nanopore reads](#minimap2-mapping)
* [Long reads SAM records Gene name tagging (GE) with dropseq.jar](#gene-name-tagging)
* [Long reads SAM records read sequence tagging (US)](#TagReadWithSequence)
* [Short reads 10xGenomics attributes cell barocde (BC) and UMI (U8) association to Long reads](#10x-attributes-association)
* [Filtering of BC/U8 associated long reads records](#FilterGetBcUmiReads)
* [Long reads Isoforms Expression assessment](#IsoformExpressionMatrix)
* [Molecule Consensus sequence computation for long read error correction](#MoleculeConsensus)
* [Mapping of molecule sequence to obtain molecular BAM file](#molecule-mapping)
* [SNPs calling, fusion transcript detection]()


## preprocessing
```bash

```

## minimap2 mapping
```bash
minimap2 -ax splice -t ... $BUILD.mmi nanopore_reads.fastq > minimap.sam
```

## gene name tagging
```bash

```

## TagReadWithSequence
```bash
java -jar -Xmx22g sicelor.jar TagReadWithSequence
```

## 10x attributes association
```bash

```

## FilterGetBcUmiReads
```bash
java -jar -Xmx22g sicelor.jar FilterGetBcUmiReads
```

## IsoformExpressionMatrix
```bash
java -jar -Xmx22g sicelor.jar IsoformExpressionMatrix
```

## MoleculeConsensus
```bash
java -jar -Xmx22g sicelor.jar MoleculeConsensus
```

## molecule mapping
```bash
minimap2 -ax splice -t ... $BUILD.mmi molecule_consensus.fasta > molecule.sam
```
