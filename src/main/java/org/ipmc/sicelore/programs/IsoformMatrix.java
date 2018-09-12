package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import java.io.*;
import java.util.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Compute expression matrices (isoforms and genes levels) and export metrics files", oneLineSummary = "Compute expression matrices (isoforms and genes levels) and export metrics files", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class IsoformMatrix extends CommandLineProgram
{
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "OUTDIR", doc = "The output directory")
    public File OUTDIR;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=10)")
    public int DELTA = 10;
    @Argument(shortName = "SOFT", doc = "Transcripts exons can be smaller than LongReadRecord exons (detection of specific alternative exons like flip/flop gria2 of Pkm1/Pkm2)")
    public boolean SOFT = false;
    @Argument(shortName = "PREFIX", doc = "Prefix for output file names (default=sicelore)")
    public String PREFIX = "sicelore";

    public HashSet<String> DTEcells;
    private final Log log;

    public IsoformMatrix() {
        log = Log.getInstance(IsoformMatrix.class);
        this.DTEcells = new HashSet<String>();
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        process();

        return 0;
    }

    protected void process()
    {
        File ISOMATRIX   = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isoforms_matrix.txt");
        File ISOMETRICS  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isoforms_metrics.txt");
        File GENEMATRIX  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_genes_matrix.txt");
        File GENEMETRICS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_genes_metrics.txt");
        File CELLMETRICS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_cells_metrics.txt");
        //File MOLMETRICS  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_ecarts_metrics.txt");

        loadDTEcells();
        log.info(new Object[]{"Cells loaded\t\t\t[" + DTEcells.size() + "]"});

        // 4mn and 9.6Gb for 1.450.000 SAMrecords [747.000 molecules]
        UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        LongreadParser bam = new LongreadParser(INPUT);
        MoleculeDataset dataset = new MoleculeDataset(bam, model, DELTA, SOFT);
        
        Matrix matrix = dataset.produceMatrix(model, DTEcells);
        matrix.writeIsoformMatrix(ISOMATRIX, ISOMETRICS);
        matrix.writeGeneMatrix(GENEMATRIX);
        matrix.writeCellMetrics(CELLMETRICS);
        matrix.writeGeneMetrics(GENEMETRICS);
        
        //dataset.writeMoleculeMetrics(MOLMETRICS);

        log.info(new Object[]{"\tMatrix cells\t\t\t[" + matrix.getCellMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix genes\t\t\t[" + matrix.getGeneMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix isoforms\t\t\t[" + matrix.getMatrice().size() + "]"});
        log.info(new Object[]{"\tMatrix total counts\t\t[" + matrix.getTotal_count() + "]"});
        log.info(new Object[]{"\tMatrix isoform def\t\t[" + matrix.getTotal_isoform_def() + "]"});
        log.info(new Object[]{"\tMatrix isoform undef\t\t[" + matrix.getTotal_isoform_undef() + "]"});
        
        //dataset.displayMetrics(METRICS);
    }

    public void loadDTEcells()
    {
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(CSV));
            String line = fichier.readLine();
            while (line != null) {
                DTEcells.add(line);
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        System.exit(new IsoformMatrix().instanceMain(args));
    }
}
