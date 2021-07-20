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
import java.util.HashSet;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Export clipped reads", oneLineSummary = "Export clipped reads.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ExportClippedReads extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output FastQ file", optional=true)
    public File OUTPUT;
    @Argument(shortName = "MINCLIP", doc = "Hard or Soft Clipping size to call as clipped read (default=150)", optional=true)
    public int MINCLIP = 150;
    
    private ProgressLogger pl;
    private final Log log;

    public ExportClippedReads() {
        log = Log.getInstance(ExportClippedReads.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        int records=0;
        int clippedRecords=0;
        BufferedOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        HashSet<String> clippedReads = new HashSet<String>();
        
        try {
            SamReader sam = SamReaderFactory.makeDefault().open(INPUT);
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            
            for(SAMRecord r : sam){
                pl.record(r);
                records++;
                String readName = r.getReadName();
                boolean isClipped = false;
                
                String cigar = r.getCigarString();
                String[] cigartype = cigar.split("[0-9]+");
                String[] cigarsize = cigar.split("[A-Z]");
                
                if (("H".equals(cigartype[1]) || "S".equals(cigartype[1])) &&  new Integer(cigarsize[0]).intValue() > MINCLIP)
                    isClipped = true;
                if (("H".equals(cigartype[cigartype.length - 1]) || "S".equals(cigartype[cigartype.length - 1])) &&  new Integer(cigarsize[cigarsize.length - 1]).intValue() > MINCLIP)
                    isClipped = true;

                if(isClipped && !clippedReads.contains(readName)){
                    clippedRecords++;
                    //log.info(new Object[]{clippedRecords+"/"+records+"\t"+cigar});
                    
                    String US = (String)r.getAttribute("US");
                    String UQ = (String)r.getAttribute("UQ");
                    os.write(new String("@"+readName+"\n"+US+"\n+\n"+UQ+"\n").getBytes());
                    clippedReads.add(readName);
                }
            }
            sam.close();
            os.close();
            
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }

        return 0;
    }

    public static void main(String[] args) {
        System.exit(new ExportClippedReads().instanceMain(args));
    }
}
