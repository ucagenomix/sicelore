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

@CommandLineProgramProperties(summary = "Tag molecule bam file with IG/BC/U8 tags contained in molecule read name", oneLineSummary = "Tag molecule bam file with IG/BC/U8 tags contained in molecule read name", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class AddBamMoleculeTags extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file with tags")
    public File OUTPUT;

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
                String[] info = str.split("\\|");
                
                if(info.length == 5){ // GENEID|TRANSCRIPTID|BC|U8|NBREADS --> isoform set before consensus
                    r.setAttribute("IG", info[0]);
                    //r.setAttribute("IT", info[1]);
                    r.setAttribute("BC", info[2]);
                    r.setAttribute("U8", info[3]);
                    // molecule longreads
                    r.setAttribute("RN", new Integer(info[4]).intValue());
                }
                else if(info.length == 4){ // GENEID|BC|U8|NBREADS --> isoform set after consensus
                    r.setAttribute("IG", info[0]);
                    r.setAttribute("BC", info[1]);
                    r.setAttribute("U8", info[2]);
                    // molecule longreads
                    r.setAttribute("RN", new Integer(info[3]).intValue());
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
