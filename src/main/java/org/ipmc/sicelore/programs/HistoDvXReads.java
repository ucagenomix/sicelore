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
import org.ipmc.sicelore.utils.*;
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

    HashMap<String, MoleculeMetrics> map;

    public HistoDvXReads()
    {
        log = Log.getInstance(HistoDvXReads.class);
        pl = new ProgressLogger(log);
        map = new HashMap<String, MoleculeMetrics>(); 
    }

    protected int doWork()
    {
        List<String> lst;
        DataOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        try {
            log.info(new Object[]{"Parsing bam file\t\tstart..."});
            for (SAMRecord r : inputSam) {
                pl.record(r);
                
                String molecule_name = r.getReadName();
                String BC = (String)r.getAttribute("BC");
                String U8 = (String)r.getAttribute("U8");
                int R1 = ((Integer) r.getAttribute("R1") != null) ? (Integer) r.getAttribute("R1") : 0;
                int R2 = ((Integer) r.getAttribute("R2") != null) ? (Integer) r.getAttribute("R2") : 0;
                int R3 = ((Integer) r.getAttribute("R3") != null) ? (Integer) r.getAttribute("R3") : 0;
                
                Float dv = (Float) r.getAttribute("dv");
                double pctId = 1.0 - dv;
                
                map.put(molecule_name, new MoleculeMetrics(R1, R2, R3, pctId));
            }
            inputSam.close();
            log.info(new Object[]{"Parsing bam file\t\t...end"});
            
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("geneId\ttranscriptId\tBC\tU8\txReads\txCleanReads\txConsensusReads\tpctId\n");
            Set cles = map.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String key = (String) it.next();
                MoleculeMetrics m = (MoleculeMetrics)this.map.get(key);
                
                // Cdk4|ENSMUST00000006911.11|CGCGTTTTCCAAACAC|TCACCCAGCA|120|108|10
                String[] keys = key.split("\\|");
                os.writeBytes(keys[0]+"\t"+keys[1]+"\t"+keys[2]+"\t"+keys[3]+"\t"+m.getXReads()+"\t"+m.getXCleanReads()+"\t"+m.getXConsensusReads()+"\t"+m.getPctId()+"\n");
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

