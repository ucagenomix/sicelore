package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 *
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.CellList;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Export cells and molecules metrics from bam file", oneLineSummary = "Export cells and molecules metrics from bam file", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ExportMetrics extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "OM", doc = "The output molecule metrics file")
    public File OM;
    @Argument(shortName = "OC", doc = "The output cell metrics file")
    public File OC;
    @Argument(shortName = "CELLTAG", doc = "The cell barcode tag (illumina=CB, long=BC)")
    public String CELLTAG = "CB";
    @Argument(shortName = "UMITAG", doc = "The UMI tag (illumina=UB, long=U8)")
    public String UMITAG = "UB";
    @Argument(shortName = "GENETAG", doc = "The gene tag (illumina=GN, long=IG or GE)")
    public String GENETAG = "GN";

    public CellList cellList;
    HashMap<String, HashMap<String, HashSet<String>>> mamap;
    HashMap<String, HashSet<String>> mamapcell;
    HashMap<String, String> mygene;

    public ExportMetrics()
    {
        log = Log.getInstance(ExportMetrics.class);
        pl = new ProgressLogger(log);
        mamap = new HashMap<String, HashMap<String, HashSet<String>>>();
        mamapcell = new HashMap<String, HashSet<String>>();
        mygene = new HashMap<String, String>();
    }

    protected int doWork()
    {
        BufferedOutputStream os = null;
        BufferedOutputStream os2 = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OM);
        IOUtil.assertFileIsWritable(OC);
        IOUtil.assertFileIsReadable(CSV);

        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"Cells detected\t\t[" + this.cellList.size() + "]"});
        
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        try {
            log.info(new Object[]{"Parsing bam file\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                
                String readName = r.getReadName();
                String BC = (String)r.getAttribute(CELLTAG);
                if(BC != null)
                    BC = BC.replace("-1","");
                String U8 = (String)r.getAttribute(UMITAG);
                String IG = (String)r.getAttribute(GENETAG);
                if(IG == null)
                    IG="nogene";
                
                if(this.cellList.contains(BC) && U8 != null){
                    if(! mamap.containsKey(BC)){
                        mamap.put(BC, new HashMap<String, HashSet<String>>());
                        mamapcell.put(BC, new HashSet<String>());
                    }

                    if(! (mamap.get(BC)).containsKey(U8))
                          (mamap.get(BC)).put(U8, new HashSet<String>());

                    ((mamap.get(BC)).get(U8)).add(readName);
                    mygene.put(BC+"-"+U8,IG);
                }
            }
            inputSam.close();
            
            log.info(new Object[]{"Parsing bam file\t...end"});
            log.info(new Object[]{"Number of cells\t\t"+mamap.size()});
            
            os = new BufferedOutputStream(new java.io.FileOutputStream(OM));
            os.write(new String("cell\tumi\tgene\tnb_read\n").getBytes());
            
            for(String key : mamap.keySet()){
                HashMap<String, HashSet<String>> cellmap = (HashMap<String, HashSet<String>>)this.mamap.get(key);
                
                for(String key2 : cellmap.keySet()){
                    HashSet<String> reads = (HashSet<String>)cellmap.get(key2);
                    os.write(new String(key+"\t"+key2+"\t"+(String)mygene.get(key+"-"+key2)+"\t"+reads.size()+"\n").getBytes());
                    
                    // add all reads in mamapcell
                    mamapcell.get(key).addAll(reads);
                }
            }
            os.close();

            int total_umis=0;
            int total_reads=0;
            os2 = new BufferedOutputStream(new java.io.FileOutputStream(OC));
            os2.write(new String("cell\tnb_read\tnb_umi\n").getBytes());
            for(String key : mamapcell.keySet()){
                HashSet<String> reads = (HashSet<String>)mamapcell.get(key);
                os2.write(new String(key+"\t"+reads.size()+"\t"+mamap.get(key).size()+"\n").getBytes());
                
                total_umis+=mamap.get(key).size();
                total_reads+=reads.size();
            }
            os2.close();
            
            log.info(new Object[]{"Number of UMIs\t\t"+total_umis});
            log.info(new Object[]{"Number of reads\t\t"+total_reads});

        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { inputSam.close();  os.close(); os2.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new ExportMetrics().instanceMain(paramArrayOfString));
    }
}