package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.LongreadRecord;
import org.ipmc.sicelore.utils.Longread;
import org.ipmc.sicelore.utils.Matrix;
import org.ipmc.sicelore.utils.Molecule;
import org.ipmc.sicelore.utils.TranscriptRecord;
import org.ipmc.sicelore.utils.UCSCRefFlatParser;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Output statistics for cellBC/UMI bam file", oneLineSummary = "Output statistics for cellBC/UMI bam file", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ReadBamStats extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "CELLTAG", doc = "The cell barcode tag (illumina=CB, long=BC)")
    public String CELLTAG = "CB";
    @Argument(shortName = "UMITAG", doc = "The UMI tag (illumina=UB, long=U8)")
    public String UMITAG = "UB";
    @Argument(shortName = "GENETAG", doc = "The gene tag (illumina=GN, long=IG or GE)")
    public String GENETAG = "GN";

    public HashSet<String> DTEcells;
    
    public HashSet<String> reads;
    public HashSet<String> unmapped;
    public HashSet<String> incells;
    public HashSet<String> incellsGN;
    public HashSet<String> outcells;
    
    public ReadBamStats() {
        log = Log.getInstance(ReadBamStats.class);
        pl = new ProgressLogger(log);
        this.DTEcells = new HashSet<String>();
        
        this.reads = new HashSet<String>();
        this.unmapped = new HashSet<String>();
        this.incells = new HashSet<String>();
        this.incellsGN = new HashSet<String>();
        this.outcells = new HashSet<String>();
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);

        loadDTEcells();        
        
        log.info(new Object[]{"Cells detected\t\t\t[" + this.DTEcells.size() + "]"});
        
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        int records = 0;
        try {
            log.info(new Object[]{"Parsing bam file\t\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                records++;
                String readName = r.getReadName();
                String referenceName = r.getReferenceName();
                String cb = ((String)r.getAttribute(CELLTAG) != null)?(String)r.getAttribute(CELLTAG):"nocelltag";
                String ub = (String)r.getAttribute(UMITAG);
                String gn = (String)r.getAttribute(GENETAG);
                cb=cb.replace("-1","");
                
                this.reads.add(readName);
                
                if("*".equals(referenceName))
                    this.unmapped.add(readName);
                
                if(this.DTEcells.contains(cb)){
                    this.incells.add(readName);
                    if(gn != null){
                        this.incellsGN.add(readName);
                    }
                }
                else
                    this.outcells.add(readName);
                
                if(records%1000000 == 0)
                    log.info(new Object[]{"[records="+records+", reads=" + this.reads.size() + ", unmapped=" + this.unmapped.size() + ", inCells=" + this.incells.size() + ", inCellsGene=" + this.incellsGN.size() + ", outCells=" + this.outcells.size() + "]"});
                
            }
            inputSam.close();
        } catch (Exception e) {
            e.printStackTrace();
        } finally { try {inputSam.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        log.info(new Object[]{"Reads\t\t\t[" + this.reads.size() + "]"});
        log.info(new Object[]{"Reads unmapped\t\t[" + this.unmapped.size() + "]"});
        log.info(new Object[]{"Reads in cells\t\t[" + this.incells.size() + "]"});
        log.info(new Object[]{"Reads in cells in gene\t[" + this.incellsGN.size() + "]"});
        log.info(new Object[]{"Reads out of cells\t[" + this.outcells.size() + "]"});
        
        return 0;
    }
    
    public void loadDTEcells()
    {
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(CSV));
            String line = fichier.readLine();
            while (line != null) {
                line=line.replace("-1","");
                DTEcells.add(line);
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new ReadBamStats().instanceMain(paramArrayOfString));
    }
}
