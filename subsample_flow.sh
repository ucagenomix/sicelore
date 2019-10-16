
# clone repository should have been done already
# git clone https://github.com/ucagenomix/sicelore.git
# cd sicelore
# ./subsample_flow.sh

# create output directory
mkdir output_dir

# parse illumina bam file
java -jar Jar/IlluminaParser-1.0.jar -i Data/190c.clta.illumina.bam -o output_dir/190c.clta.illumina.bam.obj -t Barcodes/cellBC.190.tsv -b CB -g GN -u UB

# scan nanopore reads
java -jar Jar/NanoporeReadScanner-0.5.jar -i Data/190c.clta.nanopore.reads.fq -o output_dir/190c.clta.nanopore.reads.scanned.fq

# map reads to genome
minimap2 -ax splice -uf --MD --sam-hit-only -t 4 --junc-bed Gencode/gencode.v18.mm10.junctions.bed Data/chr4.fa.gz output_dir/190c.clta.nanopore.reads.scanned.fq > output_dir/minimap.sam
samtools view -Sb output_dir/minimap.sam -o output_dir/minimap.unsorted.bam
samtools sort output_dir/minimap.unsorted.bam -o output_dir/minimap.bam
samtools index output_dir/minimap.bam

# tag reads with gene name
java -jar -Xmx12g Jar/Sicelore-1.0.jar AddGeneNameTag I=output_dir/minimap.bam O=output_dir/GE.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT
samtools index output_dir/GE.bam

# tag reads with fastq sequence
java -jar -Xmx12g Jar/Sicelore-1.0.jar AddBamReadSequenceTag I=output_dir/GE.bam O=output_dir/GEUS.bam FASTQ=Data/190c.clta.nanopore.reads.fq VALIDATION_STRINGENCY=SILENT
samtools index output_dir/GEUS.bam

# tag reads with cellBC/UMI barcodes
java -jar -Xmx12g Jar/NanoporeBC_UMI_finder-1.0.jar -i output_dir/GEUS.bam -o output_dir/GEUS10xAttributes.bam -k output_dir/190c.clta.illumina.bam.obj --maxUMIfalseMatchPercent 1 --maxBCfalseMatchPercent 5 --logFile output_dir/out.log

# generate isoform matrix
java -jar -Xmx12g Jar/Sicelore-1.0.jar IsoformMatrix DELTA=2 METHOD=STRICT GENETAG=GE I=output_dir/GEUS10xAttributes_umifound_.bam REFFLAT=Gencode/gencode.v18.mm10.refFlat.txt CSV=Barcodes/cellBC.190.tsv OUTDIR=output_dir PREFIX=sicelore VALIDATION_STRINGENCY=SILENT

# compute consensus sequence
java -jar -Xmx12g Jar/Sicelore-1.0.jar ComputeConsensus T=10 I=output_dir/GEUS10xAttributes.umifound_.bam O=output_dir/consensus.fa

# ...