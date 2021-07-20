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
import java.util.HashSet;
import java.util.List;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Add isoform annotation to reads.", oneLineSummary = "Add isoform annotation to reads.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class AddIsoBam extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "METHOD", doc = "Isoform assignment method (STRICT or SCORE)")
    public String METHOD = "STRICT";
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=2)")
    public int DELTA = 2;
    @Argument(shortName = "AMBIGUOUS_ASSIGN", doc = "Wether or not to assign the UMI to an isoform if ambiguous (default=false)", optional=true)
    public boolean AMBIGUOUS_ASSIGN = false;
    @Argument(shortName = "MAXCLIP", doc = "Maximum cliping size at both read ends to call as chimeric read (default=150)", optional=true)
    public int MAXCLIP = 150;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    
    public CellList cellList;
    private ProgressLogger pl;
    private final Log log;

    public AddIsoBam() {
        log = Log.getInstance(IsoformMatrix.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(INPUT);
        process();

        return 0;
    }

    protected void process()
    {
        LongreadRecord lrr = new LongreadRecord();
        MoleculeDataset md = new MoleculeDataset();
        
	lrr.setStaticParams(CELLTAG,UMITAG,GENETAG,"TE","UE","PE","US",MAXCLIP,"RN");
        UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        htsjdk.samtools.SamReader samReader = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
        samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, OUTPUT);
            
        try {
            for(SAMRecord r : samReader)
            {
                pl.record(r);
                lrr = LongreadRecord.fromSAMRecord(r, false);
                Longread lr = new Longread(lrr.getName());
                lr.addRecord(lrr);
                Molecule molecule = new Molecule(lr.getBarcode(), lr.getUmi(), 1);
                molecule.addLongread(lr);
                
                List<TranscriptRecord> transcripts = model.select(molecule.getGeneIdsArray());
                //if("SCORE".equals(METHOD))
                //    md.setIsoformScore(molecule, transcripts, DELTA, AMBIGUOUS_ASSIGN);
                
                if("STRICT".equals(METHOD)) 
                    md.setIsoformStrictNew(molecule, DELTA);
                
                r.setAttribute("IT", molecule.getTranscriptId());
                samFileWriter.addAlignment(r);
            }
            samReader.close();
            samFileWriter.close();
        } catch (Exception e) { e.printStackTrace(); }
    }

    public static void main(String[] args) {
        System.exit(new AddIsoBam().instanceMain(args));
    }
}
