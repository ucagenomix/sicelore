package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import java.io.*;
import htsjdk.samtools.util.*;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Compute consensus sequence per molecule.", oneLineSummary = "Compute consensus sequence per molecule.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class ComputeConsensus extends CommandLineProgram {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;

    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output .fasta file of consensus sequences")
    public File OUTPUT;
    @Argument(shortName = "T", doc = "The number of threads (default 20)")
    public int nThreads = 20;
    @Argument(shortName = "TMPDIR", doc = "Full path to TMPDIR")
    public String TMPDIR = "/share/data/scratch/sicelore/";
    @Argument(shortName = "CELLTAG", doc = "Cell barcode tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "TSOENDTAG", doc = "TSO end tag (default=TE)", optional=true)
    public String TSOENDTAG = "TE";
    @Argument(shortName = "UMIENDTAG", doc = "UMI end tag (default=UE)", optional=true)
    public String UMIENDTAG = "UE";
    @Argument(shortName = "POLYAENDTAG", doc = "PolyA end tag (default=PE)", optional=true)
    public String POLYAENDTAG = "PE";
    @Argument(shortName = "USTAG", doc = "Read sequence tag (default=US)", optional=true)
    public String USTAG = "US";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (default=150)", optional=true)
    public int MAXCLIP = 150;
    @Argument(shortName = "MAPQV0", doc = "Wether or not to keep mapqv=0 SAM records (default=false)", optional=true)
    public boolean MAPQV0 = false; // maybe better if true, need to be checkedfor futur projects, we might loose around 10% of UMI assign reads, not necessary molecule but depth
    //@Argument(shortName = "FASTQ", doc = "The .FASTQ file")
    //public File FASTQ;
    @Argument(shortName = "MAXREADS", doc = "Maximum number of reads per UMI to use for consensus sequence calling (default=20)", optional=true)
    public int MAXREADS = 20;
    @Argument(shortName = "MINPS", doc = "Default Phred Score for molecule having only 1 or 2 reads (default=3)", optional=true)
    public int MINPS = 3;
    @Argument(shortName = "MAXPS", doc = "Default Phred score for 100% base aggrement (default=20)", optional=true)
    public int MAXPS = 20;
    @Argument(shortName = "DEBUG", doc = "Debug mode, print consensus command and do not delete temp files (default=false)", optional=true)
    public boolean DEBUG = false;

    public ComputeConsensus() {
        log = Log.getInstance(ComputeConsensus.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        //IOUtil.assertFileIsReadable(FASTQ);
                
        String SPOAPATH = this.findExecutableOnPath("spoa");
        //String RACONPATH = this.findExecutableOnPath("racon");
        //String MINIMAP2PATH = this.findExecutableOnPath("minimap2");
        
        if(SPOAPATH == null)
            log.info(new Object[]{"Unable to find spoa, please add it to your PATH"});
        //else if(RACONPATH == null)
        //    log.info(new Object[]{"Unable to find racon, please add it to your PATH"});
        //else if(MINIMAP2PATH == null)
        //    log.info(new Object[]{"Unable to find minimap2, please add it to your PATH"});
        else{
            Consensus c = new Consensus();
            c.setStaticParams(MAXREADS,TMPDIR,SPOAPATH,DEBUG, MINPS);
            
            ConsensusMsa cmsa = new ConsensusMsa();
            cmsa.setStaticParams(MAXPS);
            
            LongreadRecord lrr = new LongreadRecord();
            lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,TSOENDTAG,UMIENDTAG,POLYAENDTAG,USTAG,MAXCLIP,RNTAG);
            
            //log.info(new Object[]{"loadFastq\tSTART..."});
            //FastqLoader fastq = new FastqLoader(FASTQ);
            //if(fastq.getMap().size()>0)
            //    lrr.setFastqLoader(fastq);
            
            boolean load_sequence = true;
            boolean is_gene_mandatory = false;
            boolean is_umi_mandatory = true;
            LongreadParser bam = new LongreadParser(INPUT, MAPQV0, load_sequence, is_gene_mandatory, is_umi_mandatory);
            
            MoleculeDataset dataset = new MoleculeDataset(bam);
            dataset.callConsensus(OUTPUT, nThreads);
        }
        
        return 0;
    }
    
    public static String findExecutableOnPath(String name)
    {
        for (String dirname : System.getenv("PATH").split(File.pathSeparator)) {
            File file = new File(dirname, name);
            if (file.isFile() && file.canExecute()) {
                return file.getAbsolutePath();
            }
        }
        return null;
    }
    
    public static void main(String[] args) {
        System.exit(new ComputeConsensus().instanceMain(args));
    }
}
