# Sicelore
Acronyme for [Si]ngle [Ce]ll [Lo]ng [Re]ad is a suite of tools dedicated to the bioinformatics processing, analysis 
and exploration of highly multiplexed single cell droplets-based approach experiments sequenced with Oxford Nanopore 
Technology long reads.

[![GitHub license]()]((https://github.com/hyeshik/poreplex/blob/master/LICENSE.txt))
[![Twitter Follow](https://img.shields.io/twitter/follow/kevinlebrigand.svg?style=social&logo=twitter)](https://twitter.com/kevinlebrigand)

## Summary
Typically starting with a possorted_genome.bam file from a 10xGenomics single cell cellranger pipeline, the workflow integrate several sequential steps for standard analysis and processing of sequencing runs of the same unfragmented library on Nanopore device.

## Installation
*sicelore* requires Java 8, <a href="https://github.com/lh3/minimap2">https://github.com/lh3/minimap2</a>, <a href="http://mccarrolllab.com/download/1276/">Drop-seq tools v1.13</a>

```bash
pip 
```

# Features
* [10xGenomics possorted_genome_bam.bam preprocessing]()
* [minimap2 mapping of Nanopore reads]()
>minimap2 -ax splice -t ... $BUILD.mmi nanopore_reads.fastq > minimap.sam
* [Long reads SAM records Gene name tagging (GE) with dropseq.jar]()
* [Long reads SAM records read sequence tagging (US)]()
>java -jar -Xmx22g sicelor.jar TagReadWithSequence
* [Short reads 10xGenomics attributes cell barocde (BC) and UMI (U8) association to Long reads]()
* [Filtering of BC/U8 associated long reads records]()
>java -jar -Xmx22g sicelor.jar FilterGetBcUmiReads
* [Long reads Isoforms Expression assessment]()
>java -jar -Xmx22g sicelor.jar IsoformExpressionMatrix
* [Molecule Consensus sequence computation for long read error correction]()
>java -jar -Xmx22g sicelor.jar MoleculeConsensus
* [Mapping of molecule sequence to obtain molecular BAM file]()
>minimap2 -ax splice -t ... $BUILD.mmi molecule_consensus.fasta > molecule.sam
* [SNPs calling, fusion transcript detection]()

...

## Quick Start

