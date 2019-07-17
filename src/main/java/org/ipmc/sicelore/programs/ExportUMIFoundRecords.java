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

@CommandLineProgramProperties(summary = "Export cellBC/UMI found SAMrecords.", oneLineSummary = "Export cellBC/UMI found SAMrecords.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ExportUMIFoundRecords extends CommandLineProgram {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;

    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output BAM file")
    public File OUTPUT;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";

    public ExportUMIFoundRecords() {
        log = Log.getInstance(ExportUMIFoundRecords.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
    }

    protected int doWork() {
        int ww = 0;

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader header = inputSam.getFileHeader();
        SAMFileWriter outputSam = new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        try {
            for (SAMRecord r : inputSam) {
                pl.record(r);
                String read_id = r.getReadName();
                String barcode = (String) r.getAttribute(CELLTAG);
                String umi = (String) r.getAttribute(UMITAG);

                if (barcode != null && umi != null)
                    outputSam.addAlignment(r);
            }
            inputSam.close();
            outputSam.close();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                inputSam.close();
                outputSam.close();
            } catch (Exception e) {
                System.err.println("can not close stream");
            }
        }

        return 0;
    }

    public static void main(String[] args) {
        System.exit(new ExportUMIFoundRecords().instanceMain(args));
    }

}
