#!/bin/bash

# need Java.1.8 in JAVA_HOME
java=`which $JAVA_HOME`
poa=`which poa`
racon=`which racon`
minimap2=`which minimap2`
samtools=`which samtools`
output_dir="${PWD}/output_dir"
tmp_dir="${output_dir}/tmp/"

if [ -z "$java" ] || [ -z "$poa" ] || [ -z "$samtools" ] || [ -z "$minimap2" ] || [ -z "$racon" ]
then
    echo -e "\nMissing path to required softwares:"
    echo -e "\tjava=$java"
    echo -e "\tpoa=$poa"
    echo -e "\tsamtools=$samtools"
    echo -e "\tminimap2=$minimap2"
    echo -e "\tracon=$racon"
    echo -e "\nPlease update your \$PATH and rerun.\n\n"
    exit
fi

# create output directory
mkdir $output_dir
mkdir $tmp_dir

# parse illumina bam file
$java -jar Jar/IlluminaParser-1.0.jar -i Data/190c.clta.illumina.bam -o $output_dir/190c.clta.illumina.bam.obj -t Barcodes/cellBC.190.tsv -b CB -g GN -u UB

# scan nanopore reads
$java -jar Jar/NanoporeReadScanner-0.5.jar -i Data/190c.clta.nanopore.reads.fq -o $output_dir

# map reads to genome
$minimap2 -ax splice -uf --MD --sam-hit-only -t 4 --junc-bed Gencode/gencode.v18.mm10.junctions.bed Data/chr4.fa.gz $output_dir/passed/190c.clta.nanopore.readsFWD.fq > $output_dir/minimap.sam
$samtools view -Sb $output_dir/minimap.sam -o $output_dir/minimap.unsorted.bam
$samtools sort $output_dir/minimap.unsorted.bam -o $output_dir/minimap.bam
$samtools index $output_dir/minimap.bam

# tag reads with gene name
$java -jar -Xmx4g Jar/Sicelore-1.0.jar AddGeneNameTag I=$output_dir/minimap.bam O=$output_dir/GE.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
$samtools index $output_dir/GE.bam

# tag reads with fastq sequence
$java -jar -Xmx4g Jar/Sicelore-1.0.jar AddBamReadSequenceTag I=$output_dir/GE.bam O=$output_dir/GEUS.bam FASTQ=$output_dir/passed/190c.clta.nanopore.readsFWD.fq VALIDATION_STRINGENCY=SILENT
$samtools index $output_dir/GEUS.bam

# tag reads with cellBC/UMI barcodes
$java -jar -Xmx4g Jar/NanoporeBC_UMI_finder-1.0.jar -i $output_dir/GEUS.bam -o $output_dir/GEUS10xAttributes.bam -k $output_dir/190c.clta.illumina.bam.obj --ncpu 8 --maxUMIfalseMatchPercent 1 --maxBCfalseMatchPercent 5 --logFile $output_dir/out.log
$samtools index $output_dir/GEUS10xAttributes.bam
$samtools index $output_dir/GEUS10xAttributes_umifound_.bam

# generate isoform matrix
#$java -jar -Xmx4g Jar/Sicelore-1.0.jar IsoformMatrix DELTA=2 METHOD=STRICT GENETAG=GE I=$output_dir/GEUS10xAttributes_umifound_.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt CSV=Barcodes/cellBC.190.tsv OUTDIR=$output_dir PREFIX=sicread VALIDATION_STRINGENCY=SILENT

# compute consensus sequence
$java -jar -Xmx4g Jar/Sicelore-1.0.jar ComputeConsensus T=10 I=$output_dir/GEUS10xAttributes_umifound_.bam O=$output_dir/consensus.fa TMPDIR=$tmp_dir

# map molecules to genome
$minimap2 -ax splice -uf --MD --sam-hit-only -t 4 --junc-bed Gencode/gencode.v18.mm10.junctions.bed Data/chr4.fa.gz $output_dir/consensus.fa > $output_dir/molecule.sam
$samtools view -Sb $output_dir/molecule.sam -o $output_dir/molecule.unsorted.bam
$samtools sort $output_dir/molecule.unsorted.bam -o $output_dir/molecule.bam
$samtools index $output_dir/molecule.bam

# add cellBC/UMI tags
$java -jar -Xmx4g Jar/Sicelore-1.0.jar AddBamMoleculeTags I=$output_dir/molecule.bam O=$output_dir/molecule.tags.bam
$samtools index $output_dir/molecule.tags.bam
	
# add gene name tag
$java -jar -Xmx4g Jar/Sicelore-1.0.jar AddGeneNameTag I=$output_dir/molecule.tags.bam O=$output_dir/molecule.tags.GE.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
$samtools index $output_dir/molecule.tags.GE.bam
	
# generate molecule isoform matrix
$java -jar -Xmx4g Jar/Sicelore-1.0.jar IsoformMatrix DELTA=2 METHOD=STRICT ISOBAM=true GENETAG=GE I=$output_dir/molecule.tags.GE.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt CSV=Barcodes/cellBC.190.tsv OUTDIR=$output_dir PREFIX=sicmol VALIDATION_STRINGENCY=SILENT
$samtools index $output_dir/sicmol_isobam.bam

# cleaning
cd $output_dir
rm -fr failed 190c.clta.illumina.bam.obj consensus.fa GEUS10xAttributes_umifound_.bam GEUS10xAttributes_umifound_.bam.bai molecule.tags.GE.bam molecule.tags.GE.bam.bai GE.bam GE.bam.bai GEUS.bam GEUS.bam.bai GEUS10xAttributes.bam GEUS10xAttributes.bam.bai minimap.bam minimap.bam.bai minimap.sam minimap.unsorted.bam molecule.bam molecule.bam.bai molecule.sam molecule.tags.bam molecule.tags.bam.bai molecule.unsorted.bam out.log passed tmp

