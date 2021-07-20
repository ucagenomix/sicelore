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
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.concurrent.TimeUnit;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Isoform and gene level expression matrices production.", oneLineSummary = "Isoform and gene level expression matrices production.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class RAMtest extends CommandLineProgram
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
    @Argument(shortName = "METHOD", doc = "Isoform assignment method (default=STRICT)")
    public String METHOD = "STRICT";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=GE)", optional=true)
    public String GENETAG = "GE";
    @Argument(shortName = "TSOENDTAG", doc = "TSO end tag (default=TE)", optional=true)
    public String TSOENDTAG = "TE";
    @Argument(shortName = "UMIENDTAG", doc = "UMI end tag (default=UE)", optional=true)
    public String UMIENDTAG = "UE";
    @Argument(shortName = "POLYAENDTAG", doc = "PolyA end tag (default=PE)", optional=true)
    public String POLYAENDTAG = "PE";
    @Argument(shortName = "USTAG", doc = "Read sequence tag (default=US)", optional=true)
    public String USTAG = "US";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (default=150)", optional=true)
    public int MAXCLIP = 150;
    @Argument(shortName = "AMBIGUOUS_ASSIGN", doc = "Wether or not to assign the UMI to an isoform if ambiguous (default=false)", optional=true)
    public boolean AMBIGUOUS_ASSIGN = false;
    @Argument(shortName = "MAPQV0", doc = "Wether or not to keep mapqv=0 SAM records (default=false)", optional=true)
    public boolean MAPQV0 = false;
    
    public CellList cellList;
    private ProgressLogger pl;
    private final Log log;

    public RAMtest() {
        log = Log.getInstance(RAMtest.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);

        LongreadRecord lrr = new LongreadRecord();
	lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,TSOENDTAG,UMIENDTAG,POLYAENDTAG,USTAG,MAXCLIP, RNTAG);
        
        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"\t\tCells detected\t\t[" + this.cellList.size() + "]"});
        
        if(!"STRICT".equals(METHOD)){
            log.info(new Object[]{"\tIsoform method: [" + METHOD + "] not allowed, only STRICT method allowed (SCORE disabled)"});
            return 0;
        }

        process();

        return 0;
    }

    protected void process()
    {
        LongreadParser bam = new LongreadParser(INPUT, MAPQV0, false, true, true);
        MoleculeDataset dataset = new MoleculeDataset(bam);
        dataset.initModel(REFFLAT);
        
        log.info(new Object[]{"\tstart sleeping after dataset"});
        try{
            TimeUnit.SECONDS.sleep(30);
        }catch(Exception e){}
        
        dataset.setIsoforms(DELTA, METHOD, AMBIGUOUS_ASSIGN);

        log.info(new Object[]{"\tstart sleeping after setIsoforms"});
        try{
            TimeUnit.SECONDS.sleep(30);
        }catch(Exception e){}

        File ISOMATRIX   = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isomatrix.txt");
        File ISOMETRICS  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isometrics.txt");
        File JUNCMATRIX  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_juncmatrix.txt");
        File JUNCMETRICS  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_juncmetrics.txt");
        File GENEMATRIX  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_genematrix.txt");
        File GENEMETRICS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_genemetrics.txt");
        File CELLMETRICS = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_cellmetrics.txt");
        File MOLINFOS  = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_molinfos.txt");
        File outISOBAM = new File(OUTDIR.getAbsolutePath() + "/" + PREFIX + "_isobam.bam");

        Matrix matrix = dataset.produceMatrix(this.cellList);
        
        log.info(new Object[]{"\tstart sleeping after produceMatrix"});
        try{
            TimeUnit.SECONDS.sleep(30);
        }catch(Exception e){}
        
        log.info(new Object[]{"\twriteJunctionMatrix\t[start]"});
        matrix.writeJunctionMatrix(JUNCMATRIX, JUNCMETRICS);
        log.info(new Object[]{"\twriteIsoformMatrix\t[start]"});
        matrix.writeIsoformMatrix(ISOMATRIX, ISOMETRICS, MOLINFOS, dataset.getModel());
        log.info(new Object[]{"\twriteGeneMatrix\t\t[start]"});
        matrix.writeGeneMatrix(GENEMATRIX, GENEMETRICS);
        log.info(new Object[]{"\twriteCellMetrics\t[start]"});
        matrix.writeCellMetrics(CELLMETRICS);
        
        log.info(new Object[]{"\tMatrix cells size\t[" + matrix.getCellMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix genes size\t[" + matrix.getGeneMetrics().size() + "]"});
        log.info(new Object[]{"\tMatrix junctions size\t[" + matrix.getMatriceJunction().size() + "]"});
        log.info(new Object[]{"\tMatrix isoforms size\t[" + matrix.getMatrice().size() + "]"});
        log.info(new Object[]{"\tMatrix isoforms counts\t[" + matrix.getTotal_count() + "]"});
        log.info(new Object[]{"\tMatrix isoforms define\t[" + matrix.getTotal_isoform_def() + "]"});
        log.info(new Object[]{"\tMatrix isoforms undefine[" + matrix.getTotal_isoform_undef() + "]"});
        
        log.info(new Object[]{"\tstart sleeping 2"});
        try{
            TimeUnit.SECONDS.sleep(30);
        }catch(Exception e){}

        
        if(ISOBAM){
            log.info(new Object[]{"\tProducing ISOBAM\t[true]"});
            SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
            htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
            samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
            SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, outISOBAM);
            try{
                for(SAMRecord r : samReader){
                    pl.record(r);
                    String isokey = (String)r.getAttribute(CELLTAG)+":"+(String)r.getAttribute(UMITAG);
                    
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

    public static void main(String[] args) {
        System.exit(new RAMtest().instanceMain(args));
    }
}
