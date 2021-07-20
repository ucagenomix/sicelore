package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.util.regex.Pattern;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Annotate junctions", oneLineSummary = "Annotate junctions", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class JunctionAnnotate extends CommandLineProgram {

    private final Log log;

    @Argument(shortName = "I", doc = "The input genome build .fasta file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The junctions .csv file \n#-----\nname,chromosome,strand,start,end\nNovel.174584,4,-,136892563,1368927943\n#-----")
    public File CSV;
    @Argument(shortName = "O", doc = "The output file")
    public File OUTPUT;
    
    public JunctionAnnotate() {
        log = Log.getInstance(JunctionAnnotate.class);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        process();
        return 0;
    }
    
    protected void process()
    {
        int count=0;
        String line = null;
        DataOutputStream os = null;
        HashMap<String, String> build = new HashMap<String, String>();
        
        log.info(new Object[]{"loadFasta\tSTART..."});
        
        try {
            int index=0;
            Scanner sc = new Scanner(INPUT, "UTF-8");
            line = sc.nextLine();
            String[] tmp = line.split(" ");
            String ref = tmp[0].replaceAll(">","");
            String seq="";
            while (sc.hasNextLine()) {
                line = sc.nextLine();
                index++;
                //System.out.println(index + "\t" + line);
                if(Pattern.matches("^$", line)){}
                else if(Pattern.matches("^>.*", line)){
                    System.out.println(ref + " loaded size = " + seq.length());
                    build.put(ref,seq);
                    seq="";
                    tmp = line.split(" ");
                    ref = tmp[0].replaceAll(">","");
                }
                else
                    seq += line;
            }
            build.put(ref,seq);
        } catch (Exception e) { e.printStackTrace(); }
            
            /*
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            line = fichier.readLine();
            String[] tmp = line.split(" ");
            String ref = tmp[0].replaceAll(">","");
            line = fichier.readLine();
            String seq = "";
            while(line != null) {
                if(Pattern.matches("^>.*", line)){
                    
                    System.out.println(ref + " loaded");
                    
                    build.put(ref,seq);
                    seq="";
                    tmp = line.split(" ");
                    ref = tmp[0].replaceAll(">","");
                }
                else
                    seq += line;
                
                line = fichier.readLine();
            }
            build.put(ref,seq);
            
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); System.out.println(line); }
        */
        
        
        
        log.info(new Object[]{"loadFasta\tEND..."});
        log.info(new Object[]{"loadFasta\t" + build.size() + " references"});
        log.info(new Object[]{"loadFasta\tEND..."});
        
        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(CSV));
            line = fichier.readLine();
            
            // Novel.168397,7,+,262,223,30729854-30732414|30732761-30734704|30734768-30735816|30735962-30754809
            
            while(line != null && !"".equals(line)){
                String[] tok = line.split(",");
                String[] tokall = tok[5].split("\\|");
                boolean is_all_canonical = true;
                String info = "";
                
                for(int i=0; i<tokall.length; i++){
                    String[] pos = tokall[i].split("-");
                    String donor = "na";
                    String acceptor = "na";
                    
                    if("+".equals(tok[2])){
                        donor = ((String)build.get("chr"+tok[1])).substring(new Integer(pos[0]).intValue(),new Integer(pos[0]).intValue()+2);
                        acceptor = ((String)build.get("chr"+tok[1])).substring(new Integer(pos[1]).intValue()-3,new Integer(pos[1]).intValue()-1);
                    }
                    else if("-".equals(tok[2])){
                        donor = ((String)build.get("chr"+tok[1])).substring(new Integer(pos[1]).intValue()-3,new Integer(pos[1]).intValue()-1);
                        acceptor = ((String)build.get("chr"+tok[1])).substring(new Integer(pos[0]).intValue(),new Integer(pos[0]).intValue()+2);
                        donor = revComp(donor);
                        acceptor = revComp(acceptor);
                    }
                    donor = donor.toUpperCase();
                    acceptor = acceptor.toUpperCase();
                    
                    info += donor+"-"+acceptor+"|";
                    
                    if(!"GT-AG".equals(donor+"-"+acceptor))
                        is_all_canonical=false;
                }
                info = info.substring(0,info.length()-1);
                
                System.out.println(line + "\t" + is_all_canonical + "\t" + info);
                os.writeBytes(line + "\t" + is_all_canonical + "\t" + info + "\n");
                
                line = fichier.readLine();
            }
            
            fichier.close();
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }

        log.info(new Object[]{"writeFasta\tEND..."});
    }
    
    protected static String revComp(String s)
    {
        s = new StringBuilder(s).reverse().toString();
        s = s.replaceAll("C","Z");
        s = s.replaceAll("G","C");
        s = s.replaceAll("Z","G");
        s = s.replaceAll("A","Z");
        s = s.replaceAll("T","A");
        s = s.replaceAll("Z","T");
        
        return s;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new JunctionAnnotate().instanceMain(paramArrayOfString));
    }
}
