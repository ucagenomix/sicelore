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
import org.ipmc.sicelore.utils.*;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Histogram of reference divergence function of number of longreads per molecule", oneLineSummary = "Histogram of reference divergence function of number of longreads per molecule", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class GetMoleculeMetrics extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input retag molecule SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output molecule metrics file")
    public File OUTPUT;

    HashMap<String, MoleculeMetrics> map;

    public GetMoleculeMetrics()
    {
        log = Log.getInstance(GetMoleculeMetrics.class);
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
                int RN = ((Integer) r.getAttribute("RN") != null) ? (Integer) r.getAttribute("RN") : 0;
                
                Float de = (Float) r.getAttribute("de");
                double pctId = 1.0 - de;
                
                map.put(molecule_name, new MoleculeMetrics(RN, pctId));
            }
            inputSam.close();
            log.info(new Object[]{"Parsing bam file\t\t...end"});
            
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("barcode\tumi\tgeneId\ttranscriptId\ttotal_reads\tpctIdentity\n");
            for(String key : map.keySet()){
                MoleculeMetrics m = (MoleculeMetrics)this.map.get(key);
                
                // Cdk4|ENSMUST00000006911.11|CGCGTTTTCCAAACAC|TCACCCAGCA|120
                String[] keys = key.split("\\|");
                os.writeBytes(keys[2]+"\t"+keys[3]+"\t"+keys[0]+"\t"+keys[1]+"\t"+m.getXReads()+"\t"+m.getPctId()+"\n");
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
        System.exit(new GetMoleculeMetrics().instanceMain(paramArrayOfString));
    }
}

