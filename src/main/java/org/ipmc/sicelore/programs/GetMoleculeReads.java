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

@CommandLineProgramProperties(summary = "Export reads for molecules in file", oneLineSummary = "Export reads for molecules in file", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class GetMoleculeReads extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output BAM file")
    public File OUTPUT;
    @Argument(shortName = "CSV", doc = "The molecules cell and molecular barcodes .csv file (BC,U8 one per line)")
    public File CSV;

    public GetMoleculeReads()
    {
        log = Log.getInstance(GetMoleculeReads.class);
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
                String[] tmp = line.split(",");
                mols.add(tmp[0]+"|"+tmp[1]);
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
                
                String cell_barcode = (String)r.getAttribute("BC");
                String molecule_barcode = (String)r.getAttribute("U8");

                if (mols.contains(cell_barcode+"|"+molecule_barcode) ) {
                    outputSam.addAlignment(r);
                }
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
        System.exit(new GetMoleculeReads().instanceMain(paramArrayOfString));
    }
}