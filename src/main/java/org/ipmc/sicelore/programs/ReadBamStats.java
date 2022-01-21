package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import gnu.trove.THashMap;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.File;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.CellList;
import org.ipmc.sicelore.utils.ReadInfo;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Compute reads statistics for Illumina or Nanopore bam file.", oneLineSummary = "Compute reads statistics for Illumina or Nanopore bam file.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
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

    public CellList cellList;
    public THashMap<String, ReadInfo> hashReads;
    
    public ReadBamStats() {
        log = Log.getInstance(ReadBamStats.class);
        pl = new ProgressLogger(log);
        
        this.hashReads = new THashMap<String, ReadInfo>();
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);

        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"Cells detected\t\t\t[" + this.cellList.size() + "]"});
        
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        int records = 0;
        try {
            log.info(new Object[]{"Parsing bam file\t\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                records++;
                String readName = r.getReadName();String referenceName = r.getReferenceName();
                String cb = ((String)r.getAttribute(CELLTAG) != null)?(String)r.getAttribute(CELLTAG):"nocelltag";
                String ub = (String)r.getAttribute(UMITAG);
                String gn = (String)r.getAttribute(GENETAG);
                cb=cb.replace("-1","");
                
                if(! hashReads.containsKey(readName)){
                    hashReads.put(readName, new ReadInfo());
                }
                
                if(! "*".equals(referenceName)){
                    ((ReadInfo)hashReads.get(readName)).setIs_mapped(true);
                
                    if(this.cellList.contains(cb)){
                        ((ReadInfo)hashReads.get(readName)).setIs_incell(true);
                        
                        if(ub != null){
                            ((ReadInfo)hashReads.get(readName)).setIs_umi(true);
                            
                            if(gn != null)
                                ((ReadInfo)hashReads.get(readName)).setIs_ingene(true);
                        }
                    }
                    else
                        ((ReadInfo)hashReads.get(readName)).setIs_incell(false);
                }
            }
            inputSam.close();
        } catch (Exception e) {
            e.printStackTrace();
        } finally { try {inputSam.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        Integer mapped = 0;
        Integer incell = 0;
        Integer outcell = 0;
        Integer ingene = 0;
        Integer inumi = 0;
        for(String key : hashReads.keySet()){
            ReadInfo r = (ReadInfo)hashReads.get(key);
            if(r.getIs_mapped())
                 mapped++;
            if(r.getIs_incell())
                 incell++;
            else
                outcell++;
            
            if(r.getIs_ingene())
                 ingene++;
            if(r.getIs_umi())
                 inumi++;            
        }
        log.info(new Object[]{"Total reads\t" + this.hashReads.size()});
        log.info(new Object[]{"Reads mapped\t" + mapped});
        log.info(new Object[]{"AND in cells\t" + incell});
        log.info(new Object[]{"AND with UMI\t" + inumi});
        log.info(new Object[]{"AND in gene\t" + ingene});
        
        return 0;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new ReadBamStats().instanceMain(paramArrayOfString));
    }
}
