package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.*;
import gnu.trove.THashMap;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Histogram of mean read QV from Fastq file.", oneLineSummary = "Histogram of mean read QV from Fastq file.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class HistoFastqMeanQV extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "FASTQ", doc = "The .fastq input file")
    public File FASTQ;
    @Argument(shortName = "OUTPUT", doc = "The output histogram file")
    public File OUTPUT;
 
    public HistoFastqMeanQV() {
        log = Log.getInstance(HistoFastqMeanQV.class);
    }

    protected int doWork()
    {
        DataOutputStream os = null;
        IOUtil.assertFileIsReadable(FASTQ);
        
        int[] mymeans = new int[51]; 
        
        log.info(new Object[]{"loadFastq\tSTART..."});
        FastqLoader fastqLoader = new FastqLoader(FASTQ, true);
        THashMap mapQV = fastqLoader.getMapQV();
        log.info(new Object[]{"loadFastq\t" + mapQV.size() + " reads loaded"});

        int nb=0;
        Set cles = mapQV.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String name = (String) it.next();
            char[] charArray = new String((byte[]) mapQV.get(name)).toCharArray();
            int[] asciiArray = new int[charArray.length]; //array to store the ascii codes
            
            for(int count = 0; count < charArray.length; count++)
                asciiArray[count] = (int)charArray[count] - 33;

            java.util.Arrays.sort(asciiArray);
            int medianQV = asciiArray[(int)(charArray.length/2)];
            
            mymeans[medianQV]++;
            nb++;
            if(nb%1000000 == 0)
                log.info(new Object[]{nb + " processed...\t" + new java.util.Date()});
        }
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(OUTPUT));
            os.writeBytes("meanQV\tnbReads\n");
            for(int i=1; i<31; i++)
                os.writeBytes(i+ "\t" + mymeans[i] + "\n");
            os.close();
            
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
   
            
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new HistoFastqMeanQV().instanceMain(paramArrayOfString));
    }
}
