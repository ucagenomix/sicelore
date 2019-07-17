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
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "T", doc = "The number of threads (default 20)")
    public int nThreads = 20;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=2)")
    public int DELTA = 2;
    @Argument(shortName = "RACON", doc = "Racon path")
    public String RACONPATH = "/share/apps/local/racon/bin/";
    @Argument(shortName = "MINIMAP2PATH", doc = "Minimap2 path")
    public String MINIMAP2PATH = "/share/apps/local/minimap2/";
    @Argument(shortName = "TMPDIR", doc = "TMPDIR")
    public String TMPDIR = "/share/data/scratch/sicelore/";
    @Argument(shortName = "METHOD", doc = "Isoform assignment method (STRICT or SOFT)")
    public String METHOD = "SOFT";
    @Argument(shortName = "CELLTAG", doc = "Cell barcode tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "TSOENDTAG", doc = "TSO end tag (default=TE)", optional=true)
    public String TSOENDTAG = "TE";
    @Argument(shortName = "UMIENDTAG", doc = "Cell barcode tag (default=UE)", optional=true)
    public String UMIENDTAG = "UE";
    @Argument(shortName = "USTAG", doc = "Read sequence tag (default=US)", optional=true)
    public String USTAG = "US";
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (default=150)", optional=true)
    public int MAXCLIP = 150;
    @Argument(shortName = "AMBIGUOUS_ASSIGN", doc = "Wether or not to assign the UMI to an isoform if ambiguous (default=false)", optional=true)
    public boolean AMBIGUOUS_ASSIGN = false;

    public ComputeConsensus() {
        log = Log.getInstance(ComputeConsensus.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFFLAT);

	Molecule m = new Molecule();
	m.setStaticParams(TMPDIR,RACONPATH,MINIMAP2PATH);
	LongreadRecord lrr = new LongreadRecord();
	lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,TSOENDTAG,UMIENDTAG,USTAG,MAXCLIP);
        
        if(!"STRICT".equals(METHOD) && !"SOFT".equals(METHOD)){
            log.info(new Object[]{"\tIsoform method: [" + METHOD + "] not allowed, please choose STRICT or SOFT only"});
            return 0;
        }
        //System.out.println(System.getenv("PATH"));
        
        UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        LongreadParser bam = new LongreadParser(INPUT, true, false);
        MoleculeDataset dataset = new MoleculeDataset(bam);
        dataset.setIsoforms(model, DELTA, METHOD, AMBIGUOUS_ASSIGN);
        dataset.callConsensus(OUTPUT, nThreads);

        return 0;
    }

    public static void main(String[] args) {
        System.exit(new ComputeConsensus().instanceMain(args));
    }
}
