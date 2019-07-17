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

@CommandLineProgramProperties(summary = "Molecule Bam filtering on UMI depth and isoform level annotation.", oneLineSummary = "Molecule Bam filtering on UMI depth and isoform level annotation.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class FilterMoleculeBam extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "RNTAG", doc = "The read number tag (default=RN)")
    public String RNTAG = "RN";
    @Argument(shortName = "VALUE", doc = "The value cutoff (TAG > value)")
    public int VALUE; 
    @Argument(shortName = "ISOTAG", doc = "The isoform tag (default=IT)")
    public String ISOTAG = "IT";
    @Argument(shortName = "UNDEF", doc = "Wether or not to keep molecule not define at the isoform level (ISOTAG=\"undef\") (Default=false)")
    public boolean UNDEF = false;

    public FilterMoleculeBam() {
        log = Log.getInstance(FilterMoleculeBam.class);
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
                String IT = (String)r.getAttribute(ISOTAG);
                int rn = (Integer) r.getAttribute(RNTAG);
                
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
        System.exit(new FilterMoleculeBam().instanceMain(paramArrayOfString));
    }
}
