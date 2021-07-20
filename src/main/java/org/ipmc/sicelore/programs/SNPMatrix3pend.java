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
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.*;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SNP/editing events detection/quantification (cellBC/UMI count matrix).", oneLineSummary = "SNP/editing event detection/quantification (cellBC/UMI count matrix).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SNPMatrix3pend extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    
    @Argument(shortName = "SNP", doc = "The SNP/editing events .csv file \n#-----\nchromosome,position,strand,name\n3,80692286,-,Gria2_RG\n3,80706912,-,Gria2_QR\n3,80692286|80706912,-,Gria2_RGQR\n#-----")
    public File SNP;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    
    public SNPMatrix3pend() {
        log = Log.getInstance(SNPMatrix3pend.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(SNP);
        IOUtil.assertFileIsReadable(REFFLAT);
        UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        
        try {
            int nb=0;
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(SNP));
            String line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split(",");
                String chromosome = tok[0];
                boolean strand = ("-".equals(tok[2]))?true:false;
                String gene = tok[3];
                int pos = new Integer(tok[1]).intValue();
                
                int somme=0;
                int count=0;
                
                List<TranscriptRecord> tr = model.select(new String[]{gene});
                for (int i = 0; i < tr.size(); i++) {
                    if(tr.get(i).getChrom().equals(chromosome)){
                        int d = tr.get(i).getDistanceTo3p(pos);
                        //System.out.println("\t" + d);
                        if(d >0){
                            somme+=d;
                            count++;
                        }
                    }
                }
                double avg = 0.0;
                if(count > 0){
                    avg = new Double(somme).doubleValue()/new Double(count).doubleValue();
                }
                System.out.println(line+","+(int)avg);
                
                line = fichier.readLine();
            }
        } catch (Exception e) { e.printStackTrace(); }

        return 0;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new SNPMatrix3pend().instanceMain(paramArrayOfString));
    }
}
