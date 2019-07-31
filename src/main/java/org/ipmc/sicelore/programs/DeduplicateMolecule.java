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

@CommandLineProgramProperties(summary = "Remove duplicate molecule from Fasta file.", oneLineSummary = "Remove duplicate molecule from Fasta file.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class DeduplicateMolecule extends CommandLineProgram {

    private final Log log;

    @Argument(shortName = "I", doc = "The .fasta input file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output deduplicate fasta file")
    public File OUTPUT;
 
    public DeduplicateMolecule() {
        log = Log.getInstance(DeduplicateMolecule.class);
    }

    protected int doWork()
    {
        int count=0;
        BufferedOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        HashMap<String, Molecule> map = new HashMap<String, Molecule>();
        
        log.info(new Object[]{"loadFasta\tSTART..."});
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            String id = fichier.readLine();
            
            while(id != null && !"".equals(id)) {
                String seq = fichier.readLine();
                
                if(!"null".equals(seq)){
                    count++;

                    id = id.replace(">", "");
                    id = id.replace("\\|", "-");
                    String[] ids = id.split("-");
                    String key = ids[0]+ids[1];
                    int rn = new Integer(ids[2]).intValue();

                    if(! map.containsKey(ids[0]+ids[1]))
                        map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                    else{
                        if(map.get(key).getRn() < rn)
                            map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                        else if(map.get(key).getRn() == rn){
                            if(map.get(key).getConsensus().length() < seq.length())
                                map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                        }
                    }
                }
                id = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); }

        log.info(new Object[]{"loadFasta\tEND..."});
        log.info(new Object[]{"loadFasta\t" + count + " sequences loaded"});
        log.info(new Object[]{"loadFasta\t" + map.size() + " molecules"});
        
        log.info(new Object[]{"writeFasta\tSTART..."});
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            Set cles = map.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String key = (String) it.next();
                Molecule m = (Molecule) map.get(key);
                os.write(new String(">" + m.getBarcode() + "|" + m.getUmi() + "|" + m.getRn() + "\n" + m.getConsensus() + "\n").getBytes());
            }
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }

        log.info(new Object[]{"writeFasta\tEND..."});

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new DeduplicateMolecule().instanceMain(paramArrayOfString));
    }
}
