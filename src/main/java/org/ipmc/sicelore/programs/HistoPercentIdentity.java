package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import java.math.BigDecimal;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.CellList;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Histogram of percentage of identity for molecules or reads.", oneLineSummary = "Histogram of percentage of identity for molecules or reads.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class HistoPercentIdentity extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "CELLTAG", doc = "The cell barcode tag (default=BC)")
    public String CELLTAG = "BC";
    @Argument(shortName = "DVTAG", doc = "The divergence tag (default=de (minimap2.10=dv)")
    public String DVTAG = "de";
    @Argument(shortName = "O", doc = "The output file")
    public File OUTPUT;
    
    public CellList cellList;
    public HashMap<String, Double> identity;
    
    public HistoPercentIdentity() {
        log = Log.getInstance(HistoPercentIdentity.class);
        pl = new ProgressLogger(log);
        this.identity = new HashMap<String, Double>();
    }

    protected int doWork()
    {
        DataOutputStream os = null;
        
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"Cells detected\t\t\t[" + this.cellList.size() + "]"});
        
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        int records = 0;
        try {
            log.info(new Object[]{"Parsing bam file\t\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                records++;
                String readName = r.getReadName();
                String cb = ((String)r.getAttribute(CELLTAG) != null)?(String)r.getAttribute(CELLTAG):"nocelltag";
                double id = (Float) r.getAttribute(DVTAG);
                cb=cb.replace("-1","");
                
                if(this.cellList.contains(cb)){
                    if(this.identity.get(readName) == null)
                        this.identity.put(readName, id);
                    else if(this.identity.get(readName) < id)
                        this.identity.put(readName, id);
                }
            }
            inputSam.close();
        } catch (Exception e) { e.printStackTrace(); } finally { try {inputSam.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        
        HashMap<Double, Integer> histogram = new HashMap<Double, Integer>();
        Set cles = this.identity.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()){
            String readName = (String) it.next();
            double x = 100 * (1.0 - (double)this.identity.get(readName));
            
            BigDecimal bd = new BigDecimal(x);
            bd = bd.setScale(1,BigDecimal.ROUND_DOWN);
            x = bd.doubleValue();
            
            //System.out.println(readName + ":"+  x);
            
            if(histogram.get(x) != null)
                histogram.put(x, histogram.get(x)+1);
            else
                histogram.put(x, 1);
        }
        
        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("pctIdentity\tnumber\n");
            Set cles2 = histogram.keySet();
            Iterator it2 = cles2.iterator();
            while (it2.hasNext()){
                Double interval = (Double) it2.next();
                if(interval > 80.0)
                    os.writeBytes(interval+"\t"+histogram.get(interval) + "\n");
            }

            os.close();
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new HistoPercentIdentity().instanceMain(paramArrayOfString));
    }
}
