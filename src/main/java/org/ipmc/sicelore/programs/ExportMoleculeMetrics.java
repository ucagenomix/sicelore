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

@CommandLineProgramProperties(summary = "Export molecule metrics from re-tag molecule bam file", oneLineSummary = "Export molecule metrics from re-tag molecule bam file", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ExportMoleculeMetrics extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input re-tag molecule SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output molecule metrics file")
    public File OUTPUT;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=IG)", optional=true)
    public String GENETAG = "IG";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "DVTAG", doc = "The divergence tag (default=de (minimap2.10=dv)")
    public String DVTAG = "de";

    HashMap<String, MoleculeMetrics> map;

    public ExportMoleculeMetrics()
    {
        log = Log.getInstance(ExportMoleculeMetrics.class);
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
                
                String BC = (String)r.getAttribute(CELLTAG);
                String U8 = (String)r.getAttribute(UMITAG);
                String IG = (String)r.getAttribute(GENETAG);
                int RN = (Integer) r.getAttribute(RNTAG);
                Float de = (Float) r.getAttribute(DVTAG);
                double pctId = 1.0 - de;
                
                map.put(BC+"|"+U8+"|"+IG, new MoleculeMetrics(RN, pctId));
            }
            inputSam.close();
            log.info(new Object[]{"Parsing bam file\t\t...end"});
            
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("barcode\tumi\tgeneId\ttotal_reads\tpctIdentity\n");
            for(String key : map.keySet()){
                MoleculeMetrics m = (MoleculeMetrics)this.map.get(key);
                
                String[] keys = key.split("\\|");
                os.writeBytes(keys[0]+"\t"+keys[1]+"\t"+keys[2]+"\t"+m.getXReads()+"\t"+m.getPctId()+"\n");
            }
            os.close();

        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { inputSam.close();  os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new ExportMoleculeMetrics().instanceMain(paramArrayOfString));
    }
}