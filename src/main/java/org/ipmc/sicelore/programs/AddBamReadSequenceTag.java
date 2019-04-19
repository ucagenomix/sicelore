package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.*;
import gnu.trove.THashMap;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Tag read with FASTQ sequence (US) and QV value (UQ)", oneLineSummary = "Tag read with FASTQ sequence (US) and QV value (UQ)", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
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
    @Argument(shortName = "TAG", doc = "The tag <default=US>")
    public String TAG;
    @Argument(shortName = "TAGQV", doc = "The tag <default=UQ>")
    public String TAGQV;

    public AddBamReadSequenceTag() {
        log = Log.getInstance(AddBamReadSequenceTag.class);
        pl = new ProgressLogger(log);
        TAG = "US";
        TAGQV = "UQ";
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(FASTQ);
        IOUtil.assertFileIsWritable(OUTPUT);

        log.info(new Object[]{"loadFastq\tSTART..."});
        FastqLoader localFastqLoader = new FastqLoader(FASTQ);
        THashMap localTHashMap = localFastqLoader.getMap();
        THashMap localTHashMapQV = localFastqLoader.getMapQV();
        log.info(new Object[]{"loadFastq\t" + localTHashMap.size() + " reads loaded"});

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            log.info(new Object[]{"addTag\t\tSTART..."});
            for (SAMRecord localSAMRecord : localSamReader) {
                pl.record(localSAMRecord);
                String name = localSAMRecord.getReadName();
                String seq = new String((byte[]) localTHashMap.get(name));
                localSAMRecord.setAttribute(TAG, seq);
                String qv = new String((byte[]) localTHashMapQV.get(name));
                localSAMRecord.setAttribute(TAGQV, qv);
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
