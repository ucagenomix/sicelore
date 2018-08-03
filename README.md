# Sicelore

- authors: Kevin Lebrigand & Rainer Waldmann, <a href="http://www.genomique.eu">UCAGenomiX</a>, Sophia-Antipolis, France
- current release: v0.1 (august 2018)

# summary
Acronyme for [Si]ngle [Ce]ll [Lo]ng [Re]ad is a suite of tools dedicated to the bioinformatics processing, analysis 
and exploration of highly multiplexed single cell droplets-based approach experiments sequenced with Oxford Nanopore 
Technology long reads.

# description
Typically starting with a possorted_genome.bam file from a 10xGenomics single cell cellranger pipeline, the workflow integrate several sequential steps for standard analysis and processing of sequencing runs of the same unfragmented library on Nanopore device. It requires java 8 installed, dropseq.jar

# requirements
- Java 8
- minimap2 heng li mapper (<a href="https://github.com/lh3/minimap2">https://github.com/lh3/minimap2</a>)
- Drop-seq tools v1.13 from MacCaroll's lab (<a href="http://mccarrolllab.com/download/1276/">http://mccarrolllab.com/download/1276/</a>)

# workflow
> 1. 10xGenomics possorted_genome_bam.bam preprocessing

> 2. minimap2 mapping of Nanopore reads
>minimap2 -ax splice -t ... $BUILD.mmi nanopore_reads.fastq > minimap.sam

> 3. Long reads SAM records Gene name tagging (GE) with dropseq.jar

> 4. Long reads SAM records read sequence tagging (US)
>java -jar -Xmx22g sicelor.jar TagReadWithSequence

> 5. Short reads 10xGenomics attributes cell barocde (BC) and UMI (U8) association to Long reads

> 6. Filtering of BC/U8 associated long reads records
>java -jar -Xmx22g sicelor.jar FilterGetBcUmiReads

> 7. Long reads Isoforms Expression assessment
>java -jar -Xmx22g sicelor.jar IsoformExpressionMatrix

> 8. Molecule Consensus sequence computation for long read error correction
>java -jar -Xmx22g sicelor.jar MoleculeConsensus

> 9. Mapping of molecule sequence to obtain molecular BAM file
>minimap2 -ax splice -t ... $BUILD.mmi molecule_consensus.fasta > molecule.sam

> 10. SNPs calling, fusion transcript detection

...
