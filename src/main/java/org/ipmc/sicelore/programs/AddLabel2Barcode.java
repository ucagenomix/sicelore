package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import java.io.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Add sample label to cell barcode SAM Tag", oneLineSummary = "Add sample label to cell barcode SAM Tag", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class AddLabel2Barcode extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file with label")
    public File OUTPUT;
    @Argument(shortName = "LABEL", doc = "The label to add to SAM tag")
    public File LABEL;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    
    private ProgressLogger pl;
    private final Log log;

    public AddLabel2Barcode() {
        log = Log.getInstance(AddLabel2Barcode.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log, 1000000, "\tProcessed\t", "Records");
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        process();
        return 0;
    }

    protected void process()
    {
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = inputSam.getFileHeader();
        samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, OUTPUT);
        
        try{
            for(SAMRecord r : inputSam){
                pl.record(r);
                String cellBarcode = (String)r.getAttribute(CELLTAG);
                r.setAttribute(CELLTAG, cellBarcode + "-" + LABEL);
                samFileWriter.addAlignment(r);
            }
            inputSam.close();
            samFileWriter.close();
        } catch (Exception e) { e.printStackTrace(); }
    }
  
    public static void main(String[] args) {
        System.exit(new AddLabel2Barcode().instanceMain(args));
    }
}
