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
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import java.util.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Isoform and gene level expression matrices production.", oneLineSummary = "Isoform and gene level expression matrices production.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class IsoformMatrix extends CommandLineProgram
{
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=2)")
    public int DELTA = 2;
    @Argument(shortName = "OUTDIR", doc = "The output directory")
    public File OUTDIR;
    @Argument(shortName = "PREFIX", doc = "Prefix for output file names (default=sicelore)")
    public String PREFIX = "sicelore";
    @Argument(shortName = "ISOBAM", doc = "Wether or not to produce a bam file having geneId (IG) and TranscriptId (IT) SAM flags (default=false)", optional=true)
    public boolean ISOBAM = false;
    @Argument(shortName = "METHOD", doc = "Isoform assignment method (STRICT or SOFT)")
    public String METHOD = "SOFT";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "TSOENDTAG", doc = "TSO end tag (default=TE)", optional=true)
    public String TSOENDTAG = "TE";
    @Argument(shortName = "UMIENDTAG", doc = "Cell barcode tag (default=UE)", optional=true)
    public String UMIENDTAG = "UE";
    @Argument(shortName = "USTAG", doc = "Read sequence tag (default=US)", optional=true)
    public String USTAG = "US";
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (default=150)", optional=true)
    public int MAXCLIP = 150;
    @Argument(shortName = "AMBIGUOUS_ASSIGN", doc = "Wether or not to assign the UMI to an isoform if ambiguous (default=false)", optional=true)
    public boolean AMBIGUOUS_ASSIGN = false;
    
    public HashSet<String> DTEcells;
    private ProgressLogger pl;
    private final Log log;

    public IsoformMatrix() {
        log = Log.getInstance(IsoformMatrix.class);
        pl = new ProgressLogger(log);
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
        File outISOBAM = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isobam.bam");
        File MOLMETRICS  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_molecules_metrics.txt");

	LongreadRecord lrr = new LongreadRecord();
	lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,TSOENDTAG,UMIENDTAG,USTAG,MAXCLIP);

        loadDTEcells();
        log.info(new Object[]{"\tCells loaded\t\t[" + DTEcells.size() + "]"});
        
        if(!"STRICT".equals(METHOD) && !"SOFT".equals(METHOD)){
            log.info(new Object[]{"\tIsoform method: [" + METHOD + "] not allowed, please choose STRICT or SOFT only"});
            return;
        }
        
        UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        LongreadParser bam = new LongreadParser(INPUT, false, true);
        MoleculeDataset dataset = new MoleculeDataset(bam);
        dataset.setIsoforms(model, DELTA, METHOD, AMBIGUOUS_ASSIGN);
        
        Matrix matrix = dataset.produceMatrix(model, DTEcells);
        log.info(new Object[]{"\twriteIsoformMatrix\t[start]"});
        matrix.writeIsoformMatrix(ISOMATRIX, ISOMETRICS);
        log.info(new Object[]{"\twriteGeneMatrix\t\t[start]"});
        matrix.writeGeneMatrix(GENEMATRIX);
         log.info(new Object[]{"\twriteCellMetrics\t[start]"});
        matrix.writeCellMetrics(CELLMETRICS);
         log.info(new Object[]{"\twriteGeneMetrics\t[start]"});
        matrix.writeGeneMetrics(GENEMETRICS);

        log.info(new Object[]{"\tMatrix cells\t\t[" + matrix.getCellMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix genes\t\t[" + matrix.getGeneMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix isoforms\t\t[" + matrix.getMatrice().size() + "]"});
        log.info(new Object[]{"\tMatrix total counts\t[" + matrix.getTotal_count() + "]"});
        log.info(new Object[]{"\tMatrix isoform def\t[" + matrix.getTotal_isoform_def() + "]"});
        log.info(new Object[]{"\tMatrix isoform undef\t[" + matrix.getTotal_isoform_undef() + "]"});
        
        //dataset.displayMetrics(METRICS);
        
        if(ISOBAM){
            log.info(new Object[]{"\tProducing ISOBAM\t[true]"});
            SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
            htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
            samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, outISOBAM);
            try{
                for(SAMRecord r : samReader){
                    pl.record(r);
                    String isokey = (String)r.getAttribute("BC")+":"+(String)r.getAttribute("U8");
                    
                    Molecule molecule = dataset.getMolecule(isokey);
                    if(molecule!=null){
                        r.setAttribute("IG", (molecule.getGeneId()!=null)?molecule.getGeneId():"undef");
                        r.setAttribute("IT", (molecule.getTranscriptId()!=null)?molecule.getTranscriptId():"undef");
                    }
                    else{
                        r.setAttribute("IG", "undef");
                        r.setAttribute("IT", "undef");
                    }
                    samFileWriter.addAlignment(r);
                }
                samReader.close();
                samFileWriter.close();
            } catch (Exception e) { e.printStackTrace(); }
        }
    }

    public void loadDTEcells()
    {
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(CSV));
            String line = fichier.readLine();
            while (line != null) {
                line = line.replace("-1", "");
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
