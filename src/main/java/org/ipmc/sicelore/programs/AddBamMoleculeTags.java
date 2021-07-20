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

@CommandLineProgramProperties(summary = "SAMrecord IG/BC/U8 tagging from molecule bam read names.", oneLineSummary = "SAMrecord IG/BC/U8 tagging from molecule bam read names.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class AddBamMoleculeTags extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file with tags")
    public File OUTPUT;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "ISOTAG", doc = "Isoform name tag (default=IT)", optional=true)
    public String ISOTAG = "IT";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";

    public AddBamMoleculeTags() {
        log = Log.getInstance(AddBamMoleculeTags.class);
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
                String str = r.getReadName();
                String[] info = str.split("-");
                if(info.length == 1)
                    info = str.split("\\|");
                
                if(info.length == 3){ // BC|U8|NBREADS --> v1.02
                    r.setAttribute(CELLTAG, info[0]);
                    r.setAttribute(UMITAG, info[1]);
                    r.setAttribute(RNTAG, new Integer(info[2]).intValue());
                }
                localSAMFileWriter.addAlignment(r);
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception localException) {
            localException.printStackTrace();
        }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new AddBamMoleculeTags().instanceMain(paramArrayOfString));
    }
}
