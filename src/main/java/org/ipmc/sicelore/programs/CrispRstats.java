package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.HashMap;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "CrispRstats", oneLineSummary = "CrispRstats", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class CrispRstats extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output histogram", optional=true)
    public File OUTPUT;
    @Argument(shortName = "COORD", doc = "Genomic coordinate (default=21:17608000-17610000)")
    public String COORD="21:17608000-17610000";
    
    private ProgressLogger pl;
    private final Log log;

    public HashMap<Integer, Integer> histo;
    
    public CrispRstats() {
        log = Log.getInstance(CrispRstats.class);
        pl = new ProgressLogger(log);
        
        this.histo = new HashMap<Integer, Integer>();
    }

    protected int doWork()
    {
        int records=0;
        int MAX=0;
        
        DataOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        try {
            SamReader sam = SamReaderFactory.makeDefault().open(INPUT);
            
            String chromosome = COORD.split(":")[0];
            int start = new Integer((COORD.split(":")[1]).split("-")[0]).intValue();
            int end = new Integer((COORD.split(":")[1]).split("-")[1]).intValue();
            
            SAMRecordIterator iter = sam.query(chromosome, start, end, false);
            while (iter.hasNext()) {
                SAMRecord r = iter.next();
                pl.record(r);
                records++;
                String cigar = r.getCigarString();
                String[] cigartype = cigar.split("[0-9]+");
                String[] cigarsize = cigar.split("[A-Z]");
                
                int maxdel=0;
                for(int i=0; i < cigarsize.length; i++) {
                    if("D".equals(cigartype[i])){
                        int dd = new Integer(cigarsize[i-1]).intValue();
                        if(dd > maxdel)
                            maxdel = dd; 
                    }
                }
                
                if(maxdel > MAX)
                    MAX = maxdel;
                
                //System.out.println(cigar + "\n\t--> " + maxdel + "(" + MAX +")");
                
                if(this.histo.containsKey(maxdel))
                    this.histo.put(maxdel,this.histo.get(maxdel) + 1);
                else
                    this.histo.put(maxdel,1);
            }
            sam.close();
            
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }

        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("length\tnumber\n");
            
            for(int i=0; i<=MAX; i++){
                int nn = (this.histo.get(i) != null)?this.histo.get(i):0;
                os.writeBytes(i + "\t" + nn + "\n");
            }

            os.close();
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
                
        return 0;
    }

    public static void main(String[] args) {
        System.exit(new HistoClipping().instanceMain(args));
    }
}
