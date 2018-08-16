package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin
 */

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.Longread;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Histogram of reference divergence function of number of longreads per molecule", oneLineSummary = "Histogram of reference divergence function of number of longreads per molecule", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class HistoDvXReads extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input retag molecule SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output histogram.txt")
    public File OUTPUT;

    HashMap<Integer, List<Double>> map;

    public HistoDvXReads() {
        log = Log.getInstance(HistoReadsPerMolecule.class);
        pl = new ProgressLogger(log);
        map = new HashMap<Integer, List<Double>>();
    }

    protected int doWork() {
        List<String> lst;
        DataOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        try {
            List<Double> l =null;
            
            log.info(new Object[]{"Parsing bam file\t\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                
                String read_name = r.getReadName();
                String BC = (String)r.getAttribute("BC");
                String U8 = (String)r.getAttribute("U8");
                int NN = ((Integer) r.getAttribute("NN") != null) ? (Integer) r.getAttribute("NN") : 0;
                Float dv = (Float) r.getAttribute("dv");
                
                double pctId = 1.0 - dv;
                
                //System.out.println(read_name + "\t" + dv);
                if((l = (List<Double>)map.get(new Integer(NN))) == null)
                    map.put(new Integer(NN), new ArrayList<Double>());

                ((List<Double>)map.get(new Integer(NN))).add(pctId);
            }
            inputSam.close();
            log.info(new Object[]{"Parsing bam file\t\t...end"});
            
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("xreads\tnumber\tdv\n");
            for (int i=1; i<20; i++) {
                if((l = (List<Double>)map.get(new Integer(i))) == null){
                    os.writeBytes(i + "\t0\t0\t0\n");
                    log.info(new Object[]{ i + "\t0\t0\t0" });
                }
                else{
                    int number = ((List<Double>)map.get(new Integer(i))).size();
                    float sum = 0;
                    for(Double f : ((List<Double>)map.get(new Integer(i))))
                        sum += f.doubleValue();
                    
                    float avg = sum / number;
                    os.writeBytes(i + "\t" + number + "\t" + sum + "\t"  +avg + "\n");
                    log.info(new Object[]{ i + "\t" + number + "\t" + sum + "\t"  +avg });
                }
             }
            os.close();

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                inputSam.close();
                os.close();
            } catch (Exception e) {
                System.err.println("can not close stream");
            }
        }

        return 0;

    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new HistoDvXReads().instanceMain(paramArrayOfString));
    }
}

