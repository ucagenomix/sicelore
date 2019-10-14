#!/usr/bin/perl
###################################################
# Author:  Kevin Lebrigand, Rainer Waldmann
# github:  https://github.com/ucagenomix/sicelore/
# Date:    Mon Aug 05 11:07:08 2019
# Purpose: This script shows the steps required
#	   for Single Cell Long Reads processing
# 	   for running on CentOS PBS job scheduler.
#          Tested on CentOS release 6.5 (Final)
#          might not work on all distributions.
###################################################
use strict;
use warnings;
use POSIX;
 
my ($i,$j,$job_dir,$cmd);

my $core_unit 			= 20;
my $queue			= "node20c";
my $perlpath 			= "/usr/bin/perl";
my $qusbpath			= "/opt/torque/bin/qsub";
my $samtoolspath 		= "/share/apps/local/bin/samtools";
my $javapath 			= "/share/apps/local/java/jdk1.8.0_161/bin/java";
my $fastppath			= "/share/apps/local/bin/fastp";
my $minimap2path		= "/share/apps/local/minimap2/minimap2";
my $dropseqpath			= "/share/apps/local/Drop-seq_tools-1.12/jar/";
my $sicelorepath		= "/share/apps/local/sicelore/jars/Sicelore-1.0.jar";
my $jarpath			= "/share/apps/local/sicelore/jars/";
my $picardmerger		= "/share/apps/local/picard-tools-1.119/MergeSamFiles.jar";
my $tmpdir			= "/tmp/";

my $parameters = {};
$parameters->{chunks}		= 10;
$parameters->{build}		= "mm10";
$parameters->{PREFIX} 		= "sicelore";
$parameters->{fastq}		= "<path to scanned fastq file>";# NanoporeReadScanner-0.5.jar scanned fastq file (cf. https://github.com/ucagenomix/sicelore/blob/master/README.md)
$parameters->{minimap2mmi}	= "<path to build.mmi>";# minimap2 index file for species of interest
$parameters->{csv} 		= "<path to barcodes .tsv>";# Cell Ranger cell barcode .tsv file (remove "-1")
$parameters->{obj} 		= "<path to parsedForNanopore_v0.2.obj from IlluminaParser.jar>";# IlluminaParser-1.0.jar .obj file (cf. https://github.com/ucagenomix/sicelore/blob/master/README.md)
$parameters->{GTF} 		= "<path to cellranger genes.gtf>";# Cell Ranger genes.gtf reference file
$parameters->{REFFLAT} 		= "<path to gencode.v18.refflat.txt>";# gtfToGenePred refflat file (cf. https://github.com/ucagenomix/sicelore/blob/master/README.md)
$parameters->{juncbed} 		= "<path to gencode.vM18.junctions.bed>";# paftools.js junctions bed file (cf. https://github.com/ucagenomix/sicelore/blob/master/README.md)
$parameters->{results_dir}	= "<path to results directory>";

error_check($parameters);

my @jobs;

fastqsplitter();
# detect chunks jobs
detectJobs();

if(@jobs > 0){
	minimap2();
	tagBamWithGeneExon();
        tagBamWithReadSequence();
	tagBamWithBarcodes();
	minimap2merge();
	isoformMatrix();
	
	# redefine @jobs to work 
	# on independant chromosome
	# complete for full listing
	@jobs = ('1','2','3');
	chromospliter();
	consensusCall();
	consensusMapping();
	isoformMoleculeMatrix();
}
else{ print "\n\nNo .fastq files detected !!!\n\n"; exit; }

sub error_check {
    my($parameters) = @_;
    
    my $error_to_print; 
 	my $parameters_to_print;
 	
 	unless(defined $parameters->{fastq} && -e $parameters->{fastq}) { $error_to_print .=  "\tNo valid fastq file specified.\n";}
 	unless(defined $parameters->{minimap2mmi} && -e $parameters->{minimap2mmi}) { $error_to_print .=  "\tNo valid minimap2 index file specified.\n";}
 	unless(defined $parameters->{csv} && -e $parameters->{csv}) { $error_to_print .=  "\tNo valid .csv barcodes file specified.\n";}
 	unless(defined $parameters->{obj} && -e $parameters->{obj}) { $error_to_print .=  "\tNo valid Illumina parsed file specified (.obj).\n";}
 	unless(defined $parameters->{GTF} && -e $parameters->{GTF}) { $error_to_print .=  "\tNo valid GTF file specified.\n";}
 	unless(defined $parameters->{REFFLAT} && -e $parameters->{REFFLAT}) { $error_to_print .=  "\tNo valid REFFLAT file specified.\n";}
 	unless(defined $parameters->{juncbed} && -e $parameters->{juncbed}) { $error_to_print .=  "\tNo valid junctions bed file specified.\n";}

    unless(defined $parameters->{results_dir}) { $error_to_print .= "\tNo valid results_dir specified.\n";}
    if(-e $parameters->{results_dir} && -d $parameters->{results_dir}){ 
		$parameters_to_print .=  "\tresults_dir\t= $parameters->{results_dir} \n"; 
    } 
    else{
		mkdir $parameters->{results_dir};
		mkdir $parameters->{results_dir}."/scripts/";
		mkdir $parameters->{results_dir}."/chunks/";
		mkdir $parameters->{results_dir}."/chromosomes/";
		$parameters_to_print .=   "\t$parameters->{results_dir} did not exist so it was created \n"; 
    }

    if (defined $error_to_print){ print STDERR "ERROR(s):\n", $error_to_print, "\n"; exit; }
}
sub fastqsplitter{
	my @cmds=();
	my $index = 1;
	my $ext = "split";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	push(@cmds, "mkdir ".$parameters->{results_dir}."/chunks/");
	push(@cmds, "cd ".$parameters->{results_dir}."/chunks/");
	push(@cmds, $fastppath." -i ".$parameters->{fastq}." -Q -A --thread 1 --split_prefix_digits=4 --out1=sub.fastq --split=".$parameters->{chunks});
 	createCommonJob($ext, $core_unit, @cmds);
	print JOBS $index."\t".$parameters->{results_dir}."/scripts/".$ext.".sh\n";
	
	close(JOBS);
	doJob($ext);
}
sub detectJobs {
	opendir(IMD, $parameters->{results_dir}."/chunks/") || die("Cannot open directory");
	my @files = readdir(IMD);
	for($i=0;$i<@files;$i++){
		if($files[$i] =~ /.fastq$/){
			$files[$i] =~ s/.fastq$//;
			mkdir $parameters->{results_dir}."/".$files[$i];
			mkdir $parameters->{results_dir}."/".$files[$i]."/scripts";
			push(@jobs, $files[$i]);
		}
	}	
	closedir(IMD);
}
sub createJob {
	my($job_dir,$ext,$ppn,@cmds) = @_;
	my $job_name = $job_dir."_".$ext;
	
	open SH, "> ".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh" || die $_;
	print SH "#!/bin/bash\n";
	print SH "#PBS -N ".$job_name."\n";
	print SH "#PBS -r y\n";
	print SH "#PBS -o '".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".txt'\n";
	print SH "#PBS -j oe\n";
	print SH "#PBS -q ".$queue."\n";
	print SH "#PBS -l nodes=1:ppn=".$ppn."\n";
	print SH "date '+TS[START]: %Y-%m-%d %k:%M:%S.%N'\n";
	print SH "echo StartTime is `date`\n";
	print SH "#echo Working directory is \$PBS_O_WORKDIR\n";
	print SH "#cd \$PBS_O_WORKDIR\n";
	print SH "echo Directory is `pwd`\n";  
	print SH "echo Running on host \$PBS_O_HOST / `hostname`\n";
	print SH "echo Job \$PBS_O_JOBID - \$PBS_JOBNAME in Queue \$PBS_QUEUE\n";
	print SH "umask 0002\n";
	print SH "echo 'Preparing ".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".txt'\n";
	print SH "touch '".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".txt'\n";
	print SH "chmod 0664 '".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".txt'\n";
		
	for(my $i=0;$i<@cmds;$i++){
		print SH $cmds[$i]."\n";
		print SH "echo Time is `date`\n";
	}

	print SH "RETVAL=\$?\n";
	print SH "set -e\n";
	print SH "date '+TS[JOB_END]: %Y-%m-%d %k:%M:%S.%N'\n";
	print SH "set +x\n";
	print SH "echo EndTime is `date`\n";
	print SH "date '+TS[END]: %Y-%m-%d %k:%M:%S.%N'\n";
	print SH "exit \$RETVAL\n";
	close(SH);
	
	system("/bin/chmod 775 ".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh");
}
sub createCommonJob {
	my($ext,$ppn,@cmds) = @_;
	
	open SH, "> ".$parameters->{results_dir}."/scripts/".$ext.".sh" || die $_;
	print SH "#!/bin/bash\n";
	print SH "#PBS -N ".$ext."\n";
	print SH "#PBS -r y\n";
	print SH "#PBS -o '".$parameters->{results_dir}."/scripts/".$ext.".txt'\n";
	print SH "#PBS -j oe\n";
	print SH "#PBS -q ".$queue."\n";
	print SH "#PBS -l nodes=1:ppn=".$ppn."\n";
	print SH "date '+TS[START]: %Y-%m-%d %k:%M:%S.%N'\n";
	print SH "echo StartTime is `date`\n";
	print SH "#echo Working directory is \$PBS_O_WORKDIR\n";
	print SH "#cd \$PBS_O_WORKDIR\n";
	print SH "echo Directory is `pwd`\n";  
	print SH "echo Running on host \$PBS_O_HOST / `hostname`\n";
	print SH "echo Job \$PBS_O_JOBID - \$PBS_JOBNAME in Queue \$PBS_QUEUE\n";
	print SH "umask 0002\n";
	print SH "echo 'Preparing ".$parameters->{results_dir}."/scripts/".$ext.".txt'\n";
	print SH "touch '".$parameters->{results_dir}."/scripts/".$ext.".txt'\n";
	print SH "chmod 0664 '".$parameters->{results_dir}."/scripts/".$ext.".txt'\n";
		
	for(my $i=0;$i<@cmds;$i++){
		#print $cmds[$i]."\n";
		print SH $cmds[$i]."\n";
		print SH "echo Time is `date`\n";
	}

	print SH "RETVAL=\$?\n";
	print SH "set -e\n";
	print SH "date '+TS[JOB_END]: %Y-%m-%d %k:%M:%S.%N'\n";
	print SH "set +x\n";
	print SH "echo EndTime is `date`\n";
	print SH "date '+TS[END]: %Y-%m-%d %k:%M:%S.%N'\n";
	print SH "exit \$RETVAL\n";
	close(SH);
	
	system("/bin/chmod 775 ".$parameters->{results_dir}."/scripts/".$ext.".sh");
}
sub doJob {
	my($ext) = @_;
	$cmd = $qusbpath." -j ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt";
	system($cmd);
}
sub minimap2 {
	my @cmds=();
	my $index = 0;
	my $ext = "minimap2"; 
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	for(my $i=0; $i<@jobs; $i++){
		my $job_dir = $jobs[$i];
		my $job_name = $jobs[$i]."_".$ext;
		my $pathdir = $parameters->{results_dir}."/".$job_dir."/";
		
		@cmds=();
		
		# mm10 pseudosgenes 
 		my $Nflag = "";
 		if($parameters->{build} eq "mm10"){ $Nflag = " -N 100"; }
 		
 		push(@cmds, $minimap2path." -ax splice -uf --MD ".$Nflag." --sam-hit-only -t ".$core_unit." --junc-bed ".$parameters->{juncbed}." ".$parameters->{minimap2mmi}." ".$parameters->{results_dir}."/chunks/".$jobs[$i].".fastq > ".$pathdir.$jobs[$i].".minimap.sam");
		push(@cmds, $samtoolspath." view -Sb ".$pathdir.$jobs[$i].".minimap.sam -o ".$pathdir.$jobs[$i].".minimap.unsorted.bam");
		push(@cmds, $samtoolspath." sort ".$pathdir.$jobs[$i].".minimap.unsorted.bam -o ".$pathdir.$jobs[$i].".minimap.bam");
		push(@cmds, $samtoolspath." index ".$pathdir.$jobs[$i].".minimap.bam");
 		createJob($job_dir, $ext, $core_unit, @cmds);
		$index = $i+1;
		print JOBS $index."\t".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh\n";
	}
	close(JOBS);
	doJob($ext);
}
sub tagBamWithGeneExon {
	my @cmds=();
	my $index = 0;
	my $ext = "GE";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	for(my $i=0; $i<@jobs; $i++){
		my $job_dir = $jobs[$i];
		my $job_name = $jobs[$i]."_".$ext;
		my $pathdir = $parameters->{results_dir}."/".$job_dir."/";
		
		@cmds=();
		push(@cmds, $javapath." -jar -Xmx22g ".$parameters->{sicelore}." AddGeneNameTag I=".$pathdir.$jobs[$i].".minimap.bam O=".$pathdir.$jobs[$i].".GE.bam REFFLAT=".$parameters->{REFFLAT}." GENETAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT");
		push(@cmds, "samtools index ".$pathdir.$jobs[$i].".GE.bam");
 		createJob($job_dir, $ext, $core_unit, @cmds);
		$index = $i+1;
		print JOBS $index."\t".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh\n";
	}
	close(JOBS);
	doJob($ext);
}
sub tagBamWithReadSequence {
	my @cmds=();
	my $index = 0;
	my $ext = "US";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	for(my $i=0; $i<@jobs; $i++){
		my $job_dir = $jobs[$i];
		my $job_name = $jobs[$i]."_".$ext;
		my $pathdir = $parameters->{results_dir}."/".$job_dir."/";
		
		@cmds=();
		push(@cmds, $javapath." -jar -Xmx22g ".$parameters->{sicelore}." AddBamReadSequenceTag I=".$pathdir.$jobs[$i].".GE.bam O=".$pathdir.$jobs[$i].".GEUS.bam FASTQ=".$parameters->{results_dir}."/chunks/".$jobs[$i].".fastq VALIDATION_STRINGENCY=SILENT");
 		push(@cmds, $samtoolspath." index ".$pathdir.$jobs[$i].".GEUS.bam");
 		createJob($job_dir, $ext, $core_unit, @cmds);
		$index = $i+1;
		print JOBS $index."\t".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh\n";
	}
	close(JOBS);
	doJob($ext);
}
sub tagBamWithBarcodes {
	my @cmds=();
	my $index = 0;
	my $ext = "barcodes";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	for(my $i=0; $i<@jobs; $i++){
		my $job_dir = $jobs[$i];
		my $job_name = $jobs[$i]."_".$ext;
		my $pathdir = $parameters->{results_dir}."/".$job_dir."/";
		
		@cmds=();
 		push(@cmds, "cd ".$jarpath);
		push(@cmds, $javapath." -jar -Xmx22000m NanoporeBC_UMI_finder-1.0.jar -i ".$pathdir.$jobs[$i].".GEUS.bam -o ".$pathdir.$jobs[$i].".GEUS10xAttributes.bam -k ".$parameters->{obj}." --maxUMIfalseMatchPercent 1 --maxBCfalseMatchPercent 5 --logFile ".$pathdir.$jobs[$i].".log");
  		createJob($job_dir, $ext, $core_unit, @cmds);
		$index = $i+1;
		print JOBS $index."\t".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh\n";
	}
	close(JOBS);
	doJob($ext);
}
sub minimap2merge {
	my @cmds=();
	my $index = 0;
	my $ext = "minimap2merge";
	my $bam_list = "";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	for(my $i=0; $i<@jobs; $i++){
		my $pathdir = $parameters->{results_dir}."/".$jobs[$i]."/";
		$bam_list .= "INPUT=".$pathdir.$jobs[$i].".GEUS10xAttributes_umifound_.bam ";
	}
	push(@cmds, $javapath." -jar -Xmx44g ".$picardmerger." ".$bam_list." ASSUME_SORTED=true USE_THREADING=true TMP_DIR=".$tmpdir." MAX_RECORDS_IN_RAM=100000000 OUTPUT=".$parameters->{results_dir}."/GEUS10xAttributes.umifound.bam VALIDATION_STRINGENCY=SILENT");
	push(@cmds, $samtoolspath." index ".$parameters->{results_dir}."/GEUS10xAttributes.umifound.bam");
 	createCommonJob($ext, $core_unit, @cmds);
	$index = $i+1;
	print JOBS $index."\t".$parameters->{results_dir}."/scripts/".$ext.".sh\n";
	
	close(JOBS);
	doJob($ext);
}
sub isoformMatrix {
	my @cmds=();
	my $index = 1;
	my $ext = "isoformMatrix";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	push(@cmds, $javapath." -jar -Xmx44g ".$parameters->{sicelore}." IsoformMatrix DELTA=2 METHOD=STRICT AMBIGUOUS_ASSIGN=false I=".$parameters->{results_dir}."/GEUS10xAttributes.umifound.bam REFFLAT=".$parameters->{REFFLAT}." CSV=".$parameters->{csv}." OUTDIR=".$parameters->{results_dir}." PREFIX=".$parameters->{PREFIX}." VALIDATION_STRINGENCY=SILENT");
 	createCommonJob($ext, $core_unit, @cmds);
	print JOBS $index."\t".$parameters->{results_dir}."/scripts/".$ext.".sh\n";
	
	close(JOBS);
	doJob($ext);
}
sub chromospliter {
	my @cmds=();
	my $index = 0;
	my $ext = "chrspliter";
	
	mkdir $parameters->{results_dir}."/chromosomes/";
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	for(my $i=0; $i<@jobs; $i++){
		my $job_dir = $jobs[$i];
		my $job_name = $jobs[$i]."_".$ext;
		my $pathdir = $parameters->{results_dir}."/".$job_dir."/";
		mkdir $pathdir;
		mkdir $pathdir."/scripts/";
		
		@cmds=();
 		push(@cmds, $samtoolspath." view -Sb ".$parameters->{results_dir}."/GEUS10xAttributes.umifound.bam ".$jobs[$i]." -o ".$pathdir."GEUS10xAttributes.umifound.".$jobs[$i].".bam");
		push(@cmds, $samtoolspath." index ".$pathdir."GEUS10xAttributes.umifound.".$jobs[$i].".bam");
 		createJob($job_dir, $ext, $core_unit, @cmds);
		$index = $i+1;
		print JOBS $index."\t".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh\n";
	}
	close(JOBS);
	doJob($ext);
}
sub consensusCall {
	my @cmds=();
	my $index = 0; 
	my $ext = "conscall";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	for(my $i=0; $i<@jobs; $i++){
		my $job_dir = $jobs[$i];
		my $job_name = $jobs[$i]."_".$ext;
		my $pathdir = $parameters->{results_dir}."/".$job_dir."/";
		
		@cmds=();
		push(@cmds, $javapath." -jar -Xmx22g ".$parameters->{sicelore}." ComputeConsensus T=".$core_unit." I=".$pathdir."GEUS10xAttributes.umifound.".$jobs[$i].".bam O=".$pathdir."molecules.".$jobs[$i].".fa");
 		createJob($job_dir, $ext, $core_unit, @cmds);
		$index = $i+1;
		print JOBS $index."\t".$parameters->{results_dir}."/".$job_dir."/scripts/".$job_name.".sh\n";
	}
	close(JOBS);
	doJob($ext);
}
sub consensusMapping {
	my @cmds=();
	my $index = 0;
	my $ext = "consmap"; 
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	push(@cmds, "/bin/cat ".$parameters->{results_dir}."/*/molecules.*.fa > ".$parameters->{results_dir}."/del.fa");
	push(@cmds, $javapath." -jar -Xmx44g ".$parameters->{sicelore}." DeduplicateMolecule I=".$parameters->{results_dir}."/del.fa O=".$parameters->{results_dir}."/deduplicate.fa");
	push(@cmds, $minimap2path." -ax splice -uf --sam-hit-only -t ".$core_unit." --junc-bed ".$parameters->{juncbed}." ".$parameters->{minimap2mmi}." ".$parameters->{results_dir}."/deduplicate.fa > ".$parameters->{results_dir}."/molecules.sam");
	push(@cmds, $samtoolspath." view -Sb ".$parameters->{results_dir}."/molecules.sam -o ".$parameters->{results_dir}."/molecules.unsorted.bam");
	push(@cmds, $samtoolspath." sort ".$parameters->{results_dir}."/molecules.unsorted.bam -o ".$parameters->{results_dir}."/molecules.bam");
	push(@cmds, $samtoolspath." index ".$parameters->{results_dir}."/molecules.bam");	
 	createCommonJob($ext, $core_unit, @cmds);
	$index = $i+1;
	print JOBS $index."\t".$parameters->{results_dir}."/scripts/".$ext.".sh\n";

	close(JOBS);
	doJob($ext);
}
sub isoformMoleculeMatrix {
	my @cmds=();
	my $index = 1;
	my $ext = "isoformMoleculeMatrix";
	
	open JOBS, "> ".$parameters->{results_dir}."/JOB_LIST_".$ext.".txt" || die $_;
	push(@cmds, $javapath." -jar -Xmx22g ".$parameters->{sicelore}." AddBamMoleculeTags I=".$parameters->{results_dir}."/molecules.bam O=".$parameters->{results_dir}."/molecules.tags.bam");
	push(@cmds, "samtools index ".$parameters->{results_dir}."/molecules.tags.bam");
	push(@cmds, "cd ".$dropseqpath);
	push(@cmds, $javapath." -jar -Xmx22g dropseq.jar TagReadWithGeneExon I=".$parameters->{results_dir}."molecules.tags.bam O=".$parameters->{results_dir}."molecules.tags.GE.bam ANNOTATIONS_FILE=".$parameters->{GTF}." TAG=GE ALLOW_MULTI_GENE_READS=true USE_STRAND_INFO=true VALIDATION_STRINGENCY=SILENT");
	push(@cmds, $javapath." -jar -Xmx44g ".$parameters->{sicelore}." IsoformMatrix DELTA=2 METHOD=STRICT AMBIGUOUS_ASSIGN=false I=".$parameters->{results_dir}."/molecules.tags.GE.bam REFFLAT=".$parameters->{REFFLAT}." GENETAG=GE CSV=".$parameters->{csv}." OUTDIR=".$parameters->{results_dir}." PREFIX=".$parameters->{PREFIX}.".mol.strict VALIDATION_STRINGENCY=SILENT");
	push(@cmds, $javapath." -jar -Xmx44g ".$parameters->{sicelore}." IsoformMatrix DELTA=2 METHOD=SOFT AMBIGUOUS_ASSIGN=false I=".$parameters->{results_dir}."/molecules.tags.GE.bam REFFLAT=".$parameters->{REFFLAT}." GENETAG=GE CSV=".$parameters->{csv}." OUTDIR=".$parameters->{results_dir}." PREFIX=".$parameters->{PREFIX}.".mol.soft VALIDATION_STRINGENCY=SILENT");
 	createCommonJob($ext, $core_unit, @cmds);
	print JOBS $index."\t".$parameters->{results_dir}."/scripts/".$ext.".sh\n";
	
	close(JOBS);
	doJob($ext);
}



