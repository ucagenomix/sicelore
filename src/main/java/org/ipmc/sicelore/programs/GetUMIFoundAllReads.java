package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 *
 */
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Filter associated (BC and U8 tags) molecules all SAMrecords", oneLineSummary = "Filter associated (BC and U8 tags) molecules all SAMrecords", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class GetUMIFoundAllReads extends CommandLineProgram {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;
    private HashSet<String> reads;

    @Argument(shortName = "I", doc = "The input unfiltered BAM file")
    public File INPUT;
    @Argument(shortName = "U", doc = "The input UMI found BAM file")
    public File UMIFOUND;
    @Argument(shortName = "O", doc = "The output BAM file")
    public File OUTPUT;

    public GetUMIFoundAllReads() {
        log = Log.getInstance(GetUMIFoundAllReads.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);

        reads = new HashSet<String>();
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(UMIFOUND);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        htsjdk.samtools.SamReader umifoundSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(UMIFOUND);
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);

        htsjdk.samtools.SAMFileHeader header = inputSam.getFileHeader();
        SAMFileWriter outputSam = new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        try {
            log.info(new Object[]{"Parsing UMIfound file\tstart..."});
            for (SAMRecord r : umifoundSam) {
                pl.record(r);
                reads.add(r.getReadName());
            }
            umifoundSam.close();
            log.info(new Object[]{"Parsing UMIfound file\tend..."});

            pl = new htsjdk.samtools.util.ProgressLogger(log);
            log.info(new Object[]{"Generating Molecule Reads file\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                if (reads.contains(r.getReadName())) {
                    outputSam.addAlignment(r);
                }
            }
            inputSam.close();
            outputSam.close();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                umifoundSam.close();
                inputSam.close();
                outputSam.close();
            } catch (Exception e) {
                System.err.println("can not close stream");
            }
        }

        return 0;
    }

    public static void main(String[] args) {
        System.exit(new GetUMIFoundAllReads().instanceMain(args));
    }

}
