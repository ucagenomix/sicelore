package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 *  
 */
import java.util.*;
import org.ipmc.common.utils.ExecuteCmd;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.util.concurrent.*;

public class Consensus implements Callable<String>
{
    private List<SubConsensus> sequences;
    private String name;
    private String cons;
    private String qv;
    
    protected static int MAX;
    protected static String TMPDIR;
    protected static String SPOAPATH;
    protected static boolean DEBUG;
    protected static int MINPS;
    //protected static String RACONPATH;
    //protected static String MINIMAP2PATH;
    
    private final static HashMap<Character, Integer> encode;
    private final static char[] decode;

    public Consensus() {}
    
    public Consensus(String name, List<LongreadRecord> evidenceList)
    {
        this.name = name;
        this.sequences = new ArrayList<SubConsensus>();
        
        // sort LongreadRecords on minimap2 "de" max to min value
        Collections.sort(evidenceList);
        int max = this.MAX;
        if(evidenceList.size() < this.MAX)
            max = evidenceList.size();
        
        for(int i=0; i<max; i++){
            LongreadRecord lrr = evidenceList.get(i);
            this.sequences.add(new SubConsensus(lrr.getName(),new String(lrr.getCdna())));
            //System.out.println(i+"\t"+name+"\t"+lrr.getName() + "\t" + lrr.getDe());
        }
    }
    
    public Consensus(String name, List<Longread> evidenceList, boolean from_molecule)
    {
        this.name = name;
        this.sequences = new ArrayList<SubConsensus>();
        
        // sort LongreadRecords on minimap2 "de" max to min value (hope so)
        Collections.sort(evidenceList);
        int max = this.MAX;
        if(evidenceList.size() < this.MAX)
            max = evidenceList.size();
        
        for(int i=0; i<max; i++){
            Longread lr = evidenceList.get(i);
            LongreadRecord lrr = lr.getBestRecord();
            //System.out.println(i+"\t"+name+"\t"+lrr.getName() + "\t" + lrr.getDe());
            this.sequences.add(new SubConsensus(lrr.getName(),new String(lrr.getCdna())));
        }
    }
    
    static {
        encode = new HashMap<Character, Integer>();
        encode.put('-', 0);
        encode.put('A', 1);
        encode.put('T', 2);
        encode.put('C', 3);
        encode.put('G', 4);
        encode.put('a', 1);
        encode.put('t', 2);
        encode.put('c', 3);
        encode.put('g', 4);

        decode = new char[5];
        decode[0] = '-';
        decode[1] = 'A';
        decode[2] = 'T';
        decode[3] = 'C';
        decode[4] = 'G';
    }

    public void setStaticParams(int MAX, String tmp, String spoa, boolean debug, int MINPS)
    {
        this.MAX = MAX;
        this.TMPDIR = tmp;
        this.SPOAPATH = spoa;
        this.DEBUG = debug;
        this.MINPS = MINPS;
        //this.RACONPATH = racon;
        //this.MINIMAP2PATH = minimap2;
    }

    public static int getIndexForChar(char base) { return (int) encode.get(base); }
    public static char getCharForIndex(int index) { return (char) decode[index]; }

    public List<SubConsensus> getSequences() { return sequences; }
    public String getName() { return this.name; }
    public String getCons() { return this.cons; }
    
    public String toString(){ return "name:"+this.name+",nbSeq="+ this.sequences.size(); }
    public String toFasta(){ return ">"+this.name+"\n"+ this.cons+"\n"; }
    public String toFastq(){ return "@"+this.name+"\n"+ this.cons+"\n+\n"+ this.qv+"\n"; }
    
    /*
    public String call() throws Exception
    {
        String line;
        DataOutputStream os=null;
        
        try {
            // if only one evidence -> get it
            if(this.sequences.size() == 1){
                this.cons = new String(this.sequences.get(0).getSequence());
            }
            // if 2 reads, take the longest
            else if(this.sequences.size() == 2){
                String s1 = this.sequences.get(0).getSequence();
                String s2 = this.sequences.get(1).getSequence();
                this.cons = (s1.length() > s2.length())? s1: s2;
            }
            // if more do consensus calling
            else{
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+this.name+"_reads.fa"));
                Iterator<SubConsensus> iterator2 = this.sequences.iterator();
                while(iterator2.hasNext()){
                    SubConsensus sc = (SubConsensus) iterator2.next();
                    os.writeBytes(">"+sc.getName()+"\n"+new String(sc.getSequence())+"\n");
                }
                os.close();

                // compute consensus with POA
                String[] commande = {"bash", "-c" , ""};
                String poadir = POAPATH.replaceAll("/poa$","");
                commande[2] = POAPATH + " -read_fasta "+TMPDIR+"/"+this.name+"_reads.fa -pir "+TMPDIR+"/"+this.name+".pir " + poadir + "/blosum80.mat -hb -best";
                //System.out.println(commande[2]);
                ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                // get CONSENS0 in fasta file .pir
                this.cons = "";
                BufferedReader file = new BufferedReader(new FileReader(TMPDIR+"/"+ this.name + ".pir"));
                line = file.readLine();
                while(! ">CONSENS0".equals(line)){ line = file.readLine(); }
                while(line != null && !"".equals(line)){ line = file.readLine(); this.cons += line; }                
                file.close();

                this.cons = this.cons.replaceAll("-", "");
                this.cons = this.cons.replaceAll("\\n", "");
                this.cons = this.cons.replaceAll("null", "");

                // do we still need to do that after poa consensus ?
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+this.name+"_consensus.fa"));
                os.writeBytes(">"+this.name+"\n"+this.cons+"\n");
                os.close();

                commande[2] = MINIMAP2PATH + " --secondary=no -ax map-ont "+TMPDIR+"/"+this.name+"_consensus.fa "+TMPDIR+"/"+this.name+"_reads.fa > "+TMPDIR+"/"+this.name+"_overlap.sam";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                commande[2] = RACONPATH + " "+TMPDIR+"/"+this.name+"_reads.fa "+TMPDIR+"/"+this.name+"_overlap.sam "+TMPDIR+"/"+this.name+"_consensus.fa > "+TMPDIR+"/"+this.name+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                BufferedReader fichier = new BufferedReader(new FileReader(TMPDIR+"/"+ this.name + "_corrected_consensus.fa"));
                line = fichier.readLine();
                this.cons = fichier.readLine();
                fichier.close();
                
                commande[2] = "rm "+TMPDIR+"/"+this.name+"_reads.fa "+TMPDIR+"/"+this.name+".pir "+TMPDIR+"/"+this.name+"_overlap.sam "+TMPDIR+"/"+this.name+"_consensus.fa "+TMPDIR+"/"+this.name+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                //executeCmd.run();
            }
        }catch(Exception e){ e.printStackTrace(); }
        
        return this.toFasta();
    }
    */
    
    // new method using spoa for QV calues export
    public String call() throws Exception
    {
        String line;
        DataOutputStream os=null;
        
        try {
            // if only one evidence -> get it
            if(this.sequences.size() == 1){
                this.cons = new String(this.sequences.get(0).getSequence());
                this.qv = getQvString(this.cons.length());
            }
            // if 2 reads, take the longest
            else if(this.sequences.size() == 2){
                String s1 = this.sequences.get(0).getSequence();
                String s2 = this.sequences.get(1).getSequence();
                this.cons = (s1.length() > s2.length())? s1: s2;
                this.qv = getQvString(this.cons.length());
            }
            // if more do consensus calling
            else{
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+this.name+"_reads.fa"));
                Iterator<SubConsensus> iterator2 = this.sequences.iterator();
                while(iterator2.hasNext()){
                    SubConsensus sc = (SubConsensus) iterator2.next();
                    os.writeBytes(">"+sc.getName()+"\n"+new String(sc.getSequence())+"\n");
                }
                os.close();

                // compute consensus with SPOA
                String[] commande = {"bash", "-c" , ""};
                commande[2] = SPOAPATH + " -r 2 "+TMPDIR+"/"+this.name+"_reads.fa > "+TMPDIR+"/"+this.name+".msa";
                if(this.DEBUG)
                    System.out.println(commande[2]);
                ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
                
                ConsensusMsa consensusMsa = new ConsensusMsa(TMPDIR+"/"+this.name+".msa");
                this.cons = consensusMsa.getCons();
                this.qv = consensusMsa.getQv();
                
                if(!this.DEBUG){
                    commande[2] = "rm "+TMPDIR+"/"+this.name+"_reads.fa "+TMPDIR+"/"+this.name+".msa";
                    executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                    executeCmd.run();
                }
            }
        }catch(Exception e){ e.printStackTrace(); }
        
        return this.toFastq();
    }
    
    public String getQvString(int size)
    {
        String ret = "";
        for(int i=0; i<size; i++){
            String s = Character.toString((char)(33+MINPS)); // this is the defautl QV value for 1 or 2 reads consensus, set to MINPS 
            ret += s;
        }
        return ret;
    }
}
