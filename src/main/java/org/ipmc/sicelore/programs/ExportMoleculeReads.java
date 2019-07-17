package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 *
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.HashSet;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Export reads for molecules in .csv file.", oneLineSummary = "Export reads for molecules in .csv file.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ExportMoleculeReads extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output BAM file")
    public File OUTPUT;
    @Argument(shortName = "CSV", doc = "The \"cellBC,UMI\" barcode list (.csv file)")
    public File CSV;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";

    public ExportMoleculeReads()
    {
        log = Log.getInstance(ExportMoleculeReads.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        HashSet<String> mols = new HashSet<String>();
        
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(CSV));
            String line = fichier.readLine();
            while(line != null) {
                mols.add(line);
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); }

        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader header = inputSam.getFileHeader();
        SAMFileWriter outputSam = new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);

        try {
            log.info(new Object[]{"Parsing bam file\t\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                
                String cell_barcode = (String)r.getAttribute(CELLTAG);
                String molecule_barcode = (String)r.getAttribute(UMITAG);

                if (mols.contains(cell_barcode+","+molecule_barcode) )
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

    public static void main(String[] paramArrayOfString) {
        System.exit(new ExportMoleculeReads().instanceMain(paramArrayOfString));
    }
}