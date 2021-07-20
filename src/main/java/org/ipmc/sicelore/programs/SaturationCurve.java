package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.*;
import java.util.regex.Pattern;
import htsjdk.samtools.util.*;
import java.io.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Compute saturation curve", oneLineSummary = "Compute saturation curve", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SaturationCurve extends CommandLineProgram {

    private final Log log;

    @Argument(shortName = "I", doc = "The input deduplicate molecules .fastq file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The saturation metrics")
    public File OUTPUT;
    @Argument(shortName = "READS", doc = "Total reads generated")
    public int READS = 286458673;
    @Argument(shortName = "STEP", doc = "Read step between point")
    public int STEP = 10000000;
    
    SaturationPoint[] map;
    
    public SaturationCurve() {
        log = Log.getInstance(SaturationCurve.class);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        this.initMap();
        this.printMap();
        int points = READS/STEP;
        
        int count = 0;
        String line = null;
        BufferedOutputStream os = null;
        
        log.info(new Object[]{"starting\t..."});
        try{
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            line = fichier.readLine();
            while(line != null) {
                if(Pattern.matches("^@.*", line)){
                    fichier.readLine();
                    fichier.readLine();
                    fichier.readLine();
                    
                    String[] ids = line.split("-");
                    double rn = new Double(ids[2]).doubleValue();
                    
                    for(int i=0; i<=points; i++){
                        SaturationPoint sp = map[i];
                        double seuil = sp.getProba();
                        boolean notadded = true;
                        
                        for(int j=0; j<=rn; j++){
                            if(Math.random()<seuil){
                                sp.addCountRead();
                                if(notadded){
                                    sp.addCount();
                                    notadded = false;
                                }
                            }
                        }
                    }
                }
                count++;
                if(count%1000000 == 0){
                    log.info(new Object[]{"processed " + count + " molecules..."});
                    //this.printMap();
                }
                
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); System.out.println(line); }
        log.info(new Object[]{"ending\t..."});
        
        this.printMap();
        
        return 1;
    }
    
    protected void initMap()
    {
        int points = READS/STEP;
        map = new SaturationPoint[points+1] ;
        log.info(new Object[]{"init map..."});
        
        for(int i=1; i<=points; i++){
            int reads = i*STEP;
            double proba = new Double(reads)/new Double(READS);
            map[i-1] = new SaturationPoint(i, reads, proba);
            //System.out.println(i + "," + reads + "," + proba);
        }
        map[points] = new SaturationPoint(points+1, READS, 1.0);
        //System.out.println(points + "," + READS + ",1.0");
    }
    protected void printMap()
    {
        int points = READS/STEP;
        log.info(new Object[]{"nb,reads,proba,dedup_reads,total_reads,saturation"});
        for(int i=0; i<=points; i++){
            log.info(new Object[]{map[i].toString()});
        }
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new SaturationCurve().instanceMain(paramArrayOfString));
    }
}
