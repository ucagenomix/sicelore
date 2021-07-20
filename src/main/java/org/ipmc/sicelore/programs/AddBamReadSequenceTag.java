package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SAMrecord read sequence and QV tagging.", oneLineSummary = "SAMrecord read sequence and QV tagging.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class AddBamReadSequenceTag extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output BAM file")
    public File OUTPUT;
    @Argument(shortName = "FASTQ", doc = "The .FASTQ file")
    public File FASTQ;
    @Argument(shortName = "SEQTAG", doc = "The sequence tag <default=US>")
    public String SEQTAG = "US";
    @Argument(shortName = "QVTAG", doc = "The QV tag <default=UQ>")
    public String QVTAG = "UQ";

    public AddBamReadSequenceTag() {
        log = Log.getInstance(AddBamReadSequenceTag.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsWritable(OUTPUT);

        log.info(new Object[]{"loadFastq\tSTART..."});
        FastqLoader localFastqLoader = new FastqLoader(FASTQ);
        
        log.info(new Object[]{"loadFastq\t" + localFastqLoader.getMap().size() + " reads loaded"});

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            log.info(new Object[]{"addTag\t\tSTART..."});
            for (SAMRecord localSAMRecord : localSamReader) {
                pl.record(localSAMRecord);
                String name = localSAMRecord.getReadName();
                
                //log.info(new Object[]{name});
                
                String seq = new String((byte[]) localFastqLoader.getMap().get(name));
                localSAMRecord.setAttribute(SEQTAG, seq);
                String qv = new String((byte[]) localFastqLoader.getMapQV().get(name));
                localSAMRecord.setAttribute(QVTAG, qv);
                
                localSAMFileWriter.addAlignment(localSAMRecord);
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception localException1) {
            localException1.printStackTrace();
            try {
                localSamReader.close();
                localSAMFileWriter.close();
            } catch (Exception localException2) {
                System.err.println("can not close stream");
            }
        } finally {
            try {
                localSamReader.close();
                localSAMFileWriter.close();
            } catch (Exception localException3) {
                System.err.println("can not close stream");
            }
        }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new AddBamReadSequenceTag().instanceMain(paramArrayOfString));
    }
}
