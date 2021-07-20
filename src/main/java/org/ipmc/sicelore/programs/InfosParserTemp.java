package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 *
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "InfosParserTemp", oneLineSummary = "InfosParserTemp", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class InfosParserTemp extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output file")
    public File OUTPUT;
    
    public InfosParserTemp()
    {
        log = Log.getInstance(ExportMetrics.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        BufferedOutputStream os = null;
        IOUtil.assertFileIsWritable(OUTPUT);
        
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            
            for (SAMRecord r : inputSam) {
                pl.record(r);
                
                String readName = r.getReadName();
                String CB = (String)r.getAttribute("CB");
                String UB = (String)r.getAttribute("UB");
                String RG = (String)r.getAttribute("RG");
                String GN = (String)r.getAttribute("GN");
                int txStart = r.getAlignmentStart();
                int txEnd = r.getAlignmentEnd();
                
                log.info(new Object[]{readName + "\t" + txStart+ "\t" +txEnd+ "\t" +RG + "\t" +  GN + "\t" + CB + "\t" + UB });
                os.write(new String(readName + "\t"  + txStart+ "\t" +txEnd+ "\t" + RG + "\t" +  GN + "\t" + CB + "\t" + UB+"\n").getBytes());
                
            }
            inputSam.close();
            os.close();
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { inputSam.close(); os.close();} catch (Exception e) { System.err.println("can not close stream"); } }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new InfosParserTemp().instanceMain(paramArrayOfString));
    }
}