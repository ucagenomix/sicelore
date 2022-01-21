package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.HashMap;
import java.util.HashSet;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Clipping size Histogram", oneLineSummary = "Clipping size Histogram", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class HistoClipping extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output histogram", optional=true)
    public File OUTPUT;
    @Argument(shortName = "BIN", doc = "Bin window (default=50)")
    public int BIN=50;
    @Argument(shortName = "MIN", doc = "Minimum length (default=0)")
    public int MIN=0;
    @Argument(shortName = "MAX", doc = "Minimum length (default=5000)")
    public int MAX=5000;
    
    private ProgressLogger pl;
    private final Log log;

    public HashMap<Integer, Integer> histo;
    
    public HistoClipping() {
        log = Log.getInstance(HistoClipping.class);
        pl = new ProgressLogger(log);
        
        this.histo = new HashMap<Integer, Integer>();
    }

    protected int doWork()
    {
        int records=0;
        DataOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        try {
            SamReader sam = SamReaderFactory.makeDefault().open(INPUT);
            
            for(SAMRecord r : sam){
                pl.record(r);
                records++;
                String readName = r.getReadName();
                boolean isClipped = false;
                
                String cigar = r.getCigarString();
                String[] cigartype = cigar.split("[0-9]+");
                String[] cigarsize = cigar.split("[A-Z]");
                
                int atstart=0;
                int atend = 0;
                
                if (("H".equals(cigartype[1]) || "S".equals(cigartype[1])))
                    atstart = new Integer(cigarsize[1]).intValue();
                if (("H".equals(cigartype[cigartype.length - 1]) || "S".equals(cigartype[cigartype.length - 1])))
                    atend = new Integer(cigarsize[cigarsize.length - 1]).intValue();
                
                int max_clipping = (int)(Math.max(atstart, atend))/BIN;
                if(this.histo.containsKey(max_clipping))
                    this.histo.put(max_clipping,this.histo.get(max_clipping) + 1);
                else
                    this.histo.put(max_clipping,1);
            }
            sam.close();
            
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }

        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("length\tnumber\n");
            
            int min = MIN/BIN;
            int max = MAX/BIN;
            
            for(int i=min; i<=max; i++){
                int nn = (this.histo.get(i) != null)?this.histo.get(i):0;
                os.writeBytes((i*BIN)+"\t"+ nn + "\n");
            }

            os.close();
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
                
        return 0;
    }

    public static void main(String[] args) {
        System.exit(new HistoClipping().instanceMain(args));
    }
}
