package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import gnu.trove.THashMap;
import java.io.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.List;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "HaplotypeCaller from molecules Bam file", oneLineSummary = "HaplotypeCaller from molecules Bam file", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class HaplotypeCaller extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=2)")
    public int DELTA = 2;
    @Argument(shortName = "MINEVIDENCE", doc = "Minimum evidence for Novel isoforms to be kept (default=2 molecules)")
    public int MINEVIDENCE = 10;
    @Argument(shortName = "RNMIN", doc = "Minimum number of reads to consider the UMI for Novel isoforms identification (default=1 read)")
    public int RNMIN = 5;
    @Argument(shortName = "OUTDIR", doc = "The output directory")
    public File OUTDIR;
    @Argument(shortName = "PREFIX", doc = "Prefix for output file names (default=HaplotypeCaller)")
    public String PREFIX = "HaplotypeCaller";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "ISOFORMTAG", doc = "Isoform tag (default=IT)", optional=true)
    public String ISOFORMTAG = "IT";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "T", doc = "The number of threads (default 20)")
    public int nThreads = 20;
    @Argument(shortName = "MAXUMIS", doc = "Maximum number of UMIs per isoform to use for consensus sequence calling (default=20)", optional=true)
    public int MAXUMIS = 20;
    @Argument(shortName = "MINPS", doc = "Default Phred Score for molecule having only 1 or 2 reads (default=3)", optional=true)
    public int MINPS = 3;
    @Argument(shortName = "MAXPS", doc = "Default Phred score for 100% base aggrement (default=20)", optional=true)
    public int MAXPS = 20;
    @Argument(shortName = "DEBUG", doc = "Debug mode, print consensus command and do not delete temp files (default=false)", optional=true)
    public boolean DEBUG = false;
        
    public CellList cellList;
    private ProgressLogger pl;
    private final Log log;
    
    UCSCRefFlatParser refmodel;
    UCSCRefFlatParser mymodel;

    public HaplotypeCaller()
    {
        log = Log.getInstance(HaplotypeCaller.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        
        process();
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
    
    protected void process()
    {
        this.cellList = new CellList(CSV);
        log.info(new Object[]{"\tCells detected\t\t[" + this.cellList.size() + "]"});
       
        // initialize models
        refmodel = new UCSCRefFlatParser(REFFLAT);
        mymodel = new UCSCRefFlatParser(DELTA, MINEVIDENCE, RNMIN, refmodel);
        mymodel.loader(INPUT,this.cellList,CELLTAG,UMITAG,GENETAG,ISOFORMTAG,RNTAG);
        
        THashMap<String, List<TranscriptRecord>> mapGenesTranscripts = mymodel.getMapGenesTranscripts();
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            for(int i=0; i<tAll.size(); i++){
                TranscriptRecord t = tAll.get(i);
                
                //System.out.println(geneId+"|"+t.getTranscriptId()+ " --> " + t.getEvidenceList().size());
                
                if(!"undef".equals(t.getTranscriptId()) && t.getEvidenceList().size() >= MINEVIDENCE){
                    this.exportEvidenceFastaFile(geneId+"."+t.getTranscriptId(), t.getEvidenceList());
                }
            }
        }
    }
     
    public void exportEvidenceFastaFile(String name, List<LongreadRecord> evidenceList)
    {
        File OUT = new File(OUTDIR.getAbsolutePath() + "/" + name + ".rn" + RNMIN + ".e" + MINEVIDENCE + ".fa");
        BufferedOutputStream os = null;
     
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUT));
            
            for(int i=0; i<evidenceList.size(); i++) {
                LongreadRecord lrr = evidenceList.get(i);
                if(lrr.getRn() >= RNMIN)
                    os.write(lrr.printFas().getBytes());
            }
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
    
    
    
    public static void main(String[] args) {
        System.exit(new HaplotypeCaller().instanceMain(args));
    }
}
