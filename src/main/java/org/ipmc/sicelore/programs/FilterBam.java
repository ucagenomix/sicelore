package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Filter Bam", oneLineSummary = "Filter Bam", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class FilterBam extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "TAG", doc = "The SAM Tag")
    public String TAG;
    @Argument(shortName = "VALUE", doc = "The value cutoff (TAG > value)")
    public int VALUE; 
    @Argument(shortName = "UNDEF", doc = "Wether or not the isoform is define (IT=\"undef\") (Default=false)")
    public boolean UNDEF=false;

    public FilterBam() {
        log = Log.getInstance(FilterBam.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String IT = (String)r.getAttribute("IT");
                int rn = ((Integer) r.getAttribute("R1") != null) ? (Integer) r.getAttribute("R1") : (Integer) r.getAttribute("RN");
                
                //if(rn > VALUE)
                if(rn > VALUE){
                    if(UNDEF && IT.equals("undef"))
                        localSAMFileWriter.addAlignment(r);
                    else if(! UNDEF && ! IT.equals("undef"))
                        localSAMFileWriter.addAlignment(r);
                }
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception localException) { localException.printStackTrace(); }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new FilterBam().instanceMain(paramArrayOfString));
    }
}
