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

@CommandLineProgramProperties(summary = "Compute molecule consensus sequence", oneLineSummary = "Compute molecule consensus sequence", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
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
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=10)")
    public int DELTA = 5;
    @Argument(shortName = "SOFT", doc = "Transcripts exons can be smaller than LongReadRecord exons (detection of specific alternative exons like flip/flop gria2 of Pkm1/Pkm2)")
    public boolean SOFT = false;
    @Argument(shortName = "RACON", doc = "Racon path")
    public String RACONPATH = "/share/apps/local/racon/bin/";
    @Argument(shortName = "MINIMAP2PATH", doc = "Minimap2 path")
    public String MINIMAP2PATH = "/share/apps/local/minimap2/";
    @Argument(shortName = "TMPDIR", doc = "TMPDIR")
    public String TMPDIR = "/share/data/scratch/sicelore/";
    @Argument(shortName = "CELLTAG", doc = "Cell barcode tag (optional, default=BC)")
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (optional, default=U8)")
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (optional, default=IG)")
    public String GENETAG = "IG";
    @Argument(shortName = "TSOENDTAG", doc = "TSO end tag (optional, default=TE)")
    public String TSOENDTAG = "TE";
    @Argument(shortName = "UMIENDTAG", doc = "Cell barcode tag (optional, default=UE)")
    public String UMIENDTAG = "UE";
    @Argument(shortName = "USTAG", doc = "Read sequence tag (optional, default=US)")
    public String USTAG = "US";
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (optional, default=150)")
    public int MAXCLIP = 150;

    public ComputeConsensus() {
        log = Log.getInstance(ComputeConsensus.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
    }

    protected int doWork()
    {
        int nb = 0;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFFLAT);

	Molecule m = new Molecule();
	m.setStaticParams(TMPDIR,RACONPATH,MINIMAP2PATH);
	LongreadRecord lrr = new LongreadRecord();
	lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,TSOENDTAG,UMIENDTAG,USTAG,MAXCLIP);
        
        //System.out.println(System.getenv("PATH"));
        
        UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        LongreadParser bam = new LongreadParser(INPUT, true, false);
        MoleculeDataset dataset = new MoleculeDataset(bam);
        dataset.setIsoforms(model, DELTA, SOFT);
        dataset.callConsensus(OUTPUT, nThreads);

        return 0;
    }

    public static void main(String[] args) {
        System.exit(new ComputeConsensus().instanceMain(args));
    }
}
