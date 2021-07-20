package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.HashMap;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Export edit distance for each cellBC / UMI found", oneLineSummary = "Export edit distance for each cellBC / UMI found", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class EditDistance extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "OUTDIR", doc = "The output directory")
    public File OUTDIR;
    @Argument(shortName = "PREFIX", doc = "Prefix for output file names (default=ed)")
    public String PREFIX = "ed";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "CELLEDTAG", doc = "Cell edit distance tag (default=B1)", optional=true)
    public String CELLEDTAG = "B1";
    @Argument(shortName = "UMIEDTAG", doc = "UMI edit distance tag (default=U1)", optional=true)
    public String UMIEDTAG = "U1";
    
    public CellList cellList;
    private ProgressLogger pl;
    private final Log log;

    public EditDistance() {
        log = Log.getInstance(EditDistance.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        process();

        return 0;
    }

    protected void process()
    {
        HashMap<String, Integer> edBarcode = new HashMap<String, Integer>();
        HashMap<String, Integer> edUmi = new HashMap<String, Integer>();
        
        int nb_records = 0;
        BufferedOutputStream os = null;
        File OUT = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + ".txt");
        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"\tCells detected\t\t[" + this.cellList.size() + "]"});
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        
        try{
            for(SAMRecord r : samReader){
                pl.record(r);
                nb_records++;
                String molkey = (String)r.getAttribute(CELLTAG)+":"+(String)r.getAttribute(UMITAG);
                int b1 = ((Integer) r.getAttribute(CELLEDTAG) != null) ? (Integer) r.getAttribute(CELLEDTAG) : 10;
                int u1 = ((Integer) r.getAttribute(UMIEDTAG) != null) ? (Integer) r.getAttribute(UMIEDTAG) : 10; 
                
                if(b1<10 && u1<10){
                    edBarcode.put(molkey,b1);
                    edUmi.put(molkey,u1);
                }
            }
            samReader.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        log.info(new Object[]{"\tSAM records\t\t[" + nb_records + "]"});
        log.info(new Object[]{"\tMolecules detected\t[" + edBarcode.size() + "]"});
        
        try{
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUT));
            os.write(new String("Molecule\tB1\tU1\n").getBytes());
            for(String key : edBarcode.keySet()){
                os.write(new String(key + "\t" + (int)edBarcode.get(key) + "\t" + (int)edUmi.get(key) + "\n").getBytes());
            }
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close();} catch (Exception e3) { System.err.println("can not close stream");  } }
        
    }

    public static void main(String[] args) {
        System.exit(new EditDistance().instanceMain(args));
    }
}
