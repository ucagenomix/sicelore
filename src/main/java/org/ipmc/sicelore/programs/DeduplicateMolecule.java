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
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Remove duplicate molecule from Fasta/Fastq file. Default is select for RN optimization.", oneLineSummary = "Remove duplicate molecule from Fasta/Fastq file. Default is select for RN optimization.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class DeduplicateMolecule extends CommandLineProgram {

    private final Log log;

    @Argument(shortName = "I", doc = "The input .fasta or .fastq file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output deduplicate .fasta or .fastq file")
    public File OUTPUT;
    @Argument(shortName = "TSO", doc = "TSO sequence (default=AACGCAGAGTACATGG)")
    public String TSO = "AACGCAGAGTACATGG";
    @Argument(shortName = "MAXPOS", doc = "Number of nucleotides to search for TSO sequence if still in consensus (default=100 first nt.)")
    public int MAXPOS = 100;
    @Argument(shortName = "SELECT", doc = "Wether or not select molecule optimizing RN and length.")
    public boolean SELECT = true;
 
    public int tsocut=0;
    
    public DeduplicateMolecule() {
        log = Log.getInstance(DeduplicateMolecule.class);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        if(Pattern.matches(".*\\.fq", INPUT.getName().toLowerCase()) || Pattern.matches(".*\\.fastq", INPUT.getName().toLowerCase())){
            if(SELECT)
                processFastq();
            else
                processFastqNoSelection();
        }
        else
            processFasta();
        
        return 0;
    }
    
    protected String cleanTso(String ss)
    {
        int region = (ss.length() < MAXPOS)?ss.length():MAXPOS;
        // remove remaining TSO if found in the first 100nt.
        // TSO = AAGCAGTGGTATCAACGCAGAGTACATGG
        int xx = 0;
        if((xx = ss.substring(0,region).indexOf(TSO)) > 0){
            ss = ss.substring(xx+TSO.length());
            tsocut++;
        }
        else if((xx = ss.substring(0,region).indexOf(TSO.substring(0,10))) > 0){
            ss = ss.substring(xx+16);
            tsocut++;
        }
        else if((xx = ss.substring(0,region).indexOf(TSO.substring(6,16))) > 0){
            ss = ss.substring(xx+10);
            tsocut++;
        }
        return ss;
    }
    
    protected FastqRecord cleanTso(FastqRecord f)
    {
        int region = (f.getSeq().length() < MAXPOS)?f.getSeq().length():MAXPOS;
        // remove remaining TSO if found in the first 100nt.
        // TSO = AAGCAGTGGTATCAACGCAGAGTACATGG
        int xx = 0;
        if((xx = f.getSeq().substring(0,region).indexOf(TSO)) > 0){
            f.setSeq(f.getSeq().substring(xx+TSO.length()));
            f.setQual(f.getQual().substring(xx+TSO.length()));
            tsocut++;
        }
        else if((xx = f.getSeq().substring(0,region).indexOf(TSO.substring(0,10))) > 0){
            f.setSeq(f.getSeq().substring(xx+16));
            f.setQual(f.getQual().substring(xx+16));
            tsocut++;
        }
        else if((xx = f.getSeq().substring(0,region).indexOf(TSO.substring(6,16))) > 0){
            f.setSeq(f.getSeq().substring(xx+10));
            f.setQual(f.getQual().substring(xx+10));
            tsocut++;
        }
        return f;
    }
    
    protected void processFasta()
    {
        int count=0;
        String line = null;
        BufferedOutputStream os = null;
        HashMap<String, Molecule> map = new HashMap<String, Molecule>();
        
        log.info(new Object[]{"loadFasta\tSTART..."});
        try{
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            line = fichier.readLine();
            while(line != null) {
                if(Pattern.matches("^>.*", line)){
                    String seq = fichier.readLine();
                
                    if(!"null".equals(seq)){
                        count++;
                        seq = cleanTso(seq);
                        
                        line = line.replace(">", "");
                        line = line.replace("\\|", "-");
                        String[] ids = line.split("-");
                        String key = ids[0]+ids[1];
                        int rn = new Integer(ids[2]).intValue();

                        if(! map.containsKey(key))
                            map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                        else{
                            if(map.get(key).getRn() < rn)
                                map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                            else if(map.get(key).getRn() == rn){
                                if(map.get(key).getConsensusLength() < seq.length())
                                    map.put(key, new Molecule(ids[0], ids[1], seq, rn));
                            }
                        }
                    }
                }
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); System.out.println(line); }

        log.info(new Object[]{"loadFasta\tEND..."});
        log.info(new Object[]{"loadFasta\t" + count + " sequences loaded"});
        log.info(new Object[]{"loadFasta tso\t"+tsocut});
        log.info(new Object[]{"loadFasta\t" + map.size() + " molecules"});
        log.info(new Object[]{"loadFasta\tEND..."});
        
        log.info(new Object[]{"writeFasta\tSTART..."});
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            Set cles = map.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String key = (String) it.next();
                Molecule m = (Molecule) map.get(key);
                //os.write(new String(">" + m.getBarcode() + "-" + m.getUmi() + "-" + m.getRn() + "\n" + m.getConsensus() + "\n").getBytes());
                // change 31/08/2021
                os.write(new String(">" + m.getBarcode() + "-" + m.getUmi() + "-" + m.getRn() + "\n").getBytes());
                os.write(m.getConsensus());
                os.write(new String("\n").getBytes());
                
            }
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }

        log.info(new Object[]{"writeFasta\tEND..."});
    }
    
    protected void processFastq()
    {
        int count=0;
        int tsocut=0;
        String line = null;
        BufferedOutputStream os = null;
        HashMap<String, Molecule> map = new HashMap<String, Molecule>();
        
        log.info(new Object[]{"loadFastQ\tSTART..."});
        try{
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            line = fichier.readLine();
            while(line != null) {
                if(Pattern.matches("^@.*", line)){
                    String seq = fichier.readLine();
                    String tmp = fichier.readLine();
                    String qual = fichier.readLine();
                    
                    if(!"null".equals(seq)){
                        count++;
                        
                        if(count%500000 == 0){
                            log.info(new Object[]{count + " sequences processed [" + map.size() + "]"});
                        }
                        
                        FastqRecord f = new FastqRecord("x",seq,qual);
                        //f = cleanTso(f);
                        
                        line = line.replace("@", "");
                        line = line.replace("\\|", "-");
                        String[] ids = line.split("-");
                        String key = ids[0]+ids[1];
                        int rn = new Integer(ids[2]).intValue();

                        if(! map.containsKey(key))
                            map.put(key, new Molecule(ids[0], ids[1], f.getSeq(), f.getQual(), rn));
                        else{
                            if(map.get(key).getRn() < rn)
                                map.put(key, new Molecule(ids[0], ids[1], f.getSeq(), f.getQual(), rn));
                            else if(map.get(key).getRn() == rn){
                                if(map.get(key).getConsensusLength() < f.getSeq().length())
                                    map.put(key, new Molecule(ids[0], ids[1], f.getSeq(), f.getQual(), rn));
                            }
                        }
                    }
                }
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); System.out.println(line); }

        log.info(new Object[]{"loadFastQ\tEND..."});
        log.info(new Object[]{"loadFastQ\t" + count + " sequences loaded"});
        log.info(new Object[]{"loadFastQ tso\t"+tsocut});
        log.info(new Object[]{"loadFastQ\t" + map.size() + " molecules"});
        log.info(new Object[]{"loadFastQ\tEND..."});
        
        log.info(new Object[]{"writeFastQ\tSTART..."});
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            Set cles = map.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String key = (String) it.next();
                Molecule m = (Molecule) map.get(key);
                
                // change 31/08/2021
                os.write(new String("@" + m.getBarcode() + "-" + m.getUmi() + "-" + m.getRn() + "\n").getBytes());
                os.write(m.getConsensus());
                os.write(new String("\n+\n").getBytes());
                os.write(m.getConsensusQV());
                os.write(new String("\n").getBytes());
            }
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }

        log.info(new Object[]{"writeFastQ\tEND..."});
    }
    
    protected void processFastqNoSelection()
    {
        int count=0;
        int tsocut=0;
        String line = null;
        BufferedOutputStream os = null;
        HashMap<String, Molecule> map = new HashMap<String, Molecule>();
        
        log.info(new Object[]{"load/write FastQ\tSTART..."});
        try{
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            os = new BufferedOutputStream(new java.io.FileOutputStream(OUTPUT));
            
            line = fichier.readLine();
            while(line != null) {
                if(Pattern.matches("^@.*", line)){
                    String seq = fichier.readLine();
                    String tmp = fichier.readLine();
                    String qual = fichier.readLine();
                    
                    if(!"null".equals(seq)){
                        count++;
                        FastqRecord f = new FastqRecord("x",seq,qual);
                        f = cleanTso(f);
                        
                        line = line.replace("@", "");
                        line = line.replace("\\|", "-");
                        String[] ids = line.split("-");
                        String key = ids[0]+ids[1];
                        int rn = new Integer(ids[2]).intValue();

                        if(! map.containsKey(key)){
                            map.put(key, new Molecule(ids[0], ids[1], "", "", rn));
                            os.write(new String("@" + ids[0] + "-" + ids[1] + "-" + rn + "\n" + f.getSeq() + "\n+\n" + f.getQual() + "\n").getBytes());
                        }
                    }
                }
                line = fichier.readLine();
            }
            fichier.close();
            os.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }

        log.info(new Object[]{"load/write FastQ\tEND..."});
        log.info(new Object[]{"loadFastQ\t" + count + " sequences loaded"});
        log.info(new Object[]{"loadFastQ tso\t"+tsocut});
        log.info(new Object[]{"loadFastQ\t" + map.size() + " molecules"});
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new DeduplicateMolecule().instanceMain(paramArrayOfString));
    }
}
