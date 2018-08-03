# Sicelore

- authors: Kevin Lebrigand & Rainer Waldmann, <a href="http://www.genomique.eu">UCAGenomiX</a>, Sophia-Antipolis, France
- current release: v0.1 (august 2018)

# summary
Acronyme for [Si]ngle [Ce]ll [Lo]ng [Re]ad is a suite of tools dedicated to the bioinformatics processing, analysis 
and exploration of highly multiplexed single cell droplets-based approach experiments sequenced with Oxford Nanopore 
Technology long reads.

# features

USAGE: SiCeLoReMain <program name> [-h]

Available Programs:
--------------------------------------------------------------------------------------
SiCeLoRe Pipeline:                               SiCeLoRe Pipeline - [Si]ngle [Ce]ll [Lo]ng [Re]ads tools
    IsoformExpressionMatrix                      Produce Isoforms Expression Matrix
    MoleculeConsensus                            Procude molecule consensus sequence
    ReTagReadWithBarcodes                        Tag molecule bam file with IG/BC/U8 tags contained in molecule read name
    TagReadWithSequence                          Tag read with fastq sequence

--------------------------------------------------------------------------------------
SiCeLoRe Utils:                                  SiCeLoRe Utils                               
    CollapseToGeneExpressionMatrix               Collapse Isoforms Expression Matrix to Genes Expression Matrix
    FilterGetBcUmiReads                          Filter for reads associated with an Illumina cell barcode and UMI
    FilterGetMoleculeReads                       Filter for reads that have at least a SAMrecord containing GI, BC and a U8 tags
    SplitBamPerCell                              Split a bam cells-by-cells provided in .csv file
    SplitBamPerCluster                           Split a bam according to the cell types provided in .csv file

