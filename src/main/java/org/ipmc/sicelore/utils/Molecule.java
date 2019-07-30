package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.util.*;
import org.ipmc.common.utils.ExecuteCmd;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.util.concurrent.*;

public class Molecule implements Callable<String>
{
    private List<Longread> longreads;
    private HashSet<String> geneIds;
    private String barcode;
    private String umi;
    private int rn;
    private String consensus = "";
    private String geneId = "undef";
    private String transcriptId = "undef";
    
    protected static String TMPDIR;
    protected static String POAPATH;
    protected static String RACONPATH;
    protected static String MINIMAP2PATH;
    
    //private int nbReads = 0;
    //private int nbCleanReads = 0;
    //private int nbConsensusReads = 0;
    //private String pctEcartMedian = "";
    
    private final static HashMap<Character, Integer> encode;
    private final static char[] decode;

    public Molecule() {}
    
    public Molecule(String barcode, String umi) {
        this.longreads = new ArrayList<Longread>();
        this.geneIds = new HashSet<String>();
        this.barcode = barcode;
        this.umi = umi;
    }
    
    public Molecule(String barcode, String umi, String consensus, int rn) {
        //this.longreads = new ArrayList<Longread>();
        //this.geneIds = new HashSet<String>();
        this.barcode = barcode;
        this.consensus = consensus;
        this.umi = umi;
        this.rn=rn;
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

    public void setStaticParams(String tmp, String poa, String racon, String minimap2){
        this.TMPDIR = tmp;
        this.POAPATH = poa;
        this.RACONPATH = racon;
        this.MINIMAP2PATH = minimap2;
        //System.out.println(TMPDIR);
    }

    public static int getIndexForChar(char base) {
        return (int) encode.get(base);
    }

    public static char getCharForIndex(int index) {
        return (char) decode[index];
    }

    public List<Longread> getLongreads() {
        return longreads;
    }

    public HashSet<String> getGeneIds() {
        return geneIds;
    }
    
    //public String getPctEcartMedian() {
    //    return pctEcartMedian;
    //}
        
    public String[] getGeneIdsArray() {
        String[] array = new String[this.geneIds.size()];
        this.geneIds.toArray(array);
        return array;
    }

    public String getBarcode() {
        return this.barcode;
    }

    public String getUmi() {
        return this.umi;
    }

    public int getRn() {
        return this.rn;
    }

    public String getConsensus() {
        return this.consensus;
    }

    public String getGeneId() {
        return this.geneId;
    }

    public void setRn(int rn) {
        this.rn = rn;
    }

    public String getTranscriptId() {
        return this.transcriptId;
    }
    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public void setTranscriptId(String transcriptId) {
        this.transcriptId=transcriptId;
    }
    
    public String getLabel() {
        return this.barcode + "-" + this.umi + "-" + this.longreads.size();
    }
    
    public String toString() {
        return getLabel();
    }
    /*
    public int getNbReads() {
        return nbReads;
    }
    public int getNbCleanReads() {
        return nbCleanReads;
    }
    public int getNbConsensusReads() {
        return nbConsensusReads;
    }
    */
    public void addLongread(Longread lr)
    {
        this.longreads.add(lr);
        
        Iterator<String> iterator = lr.getGeneIds().iterator();
        while (iterator.hasNext())
            this.geneIds.add((String)iterator.next());
    }

    public char getConsensusBase(int[] counts) {
        int retVal = 0;
        int max = 0;
        for (int i = 0; i < 5; i++) {
            if (counts[i] > max) {
                max = counts[i];
                retVal = i;
            }
        }
        return getCharForIndex(retVal);
    }

    public String call() throws Exception
    {
        String line;
        DataOutputStream os=null;
        //DataOutputStream os10=null;
        String prefix = this.barcode + "" + this.umi;
        
        // if only one read -> get best record as consensus
        if(this.longreads.size() == 1){
            LongreadRecord lrr = this.longreads.get(0).getBestRecord();
            this.consensus = new String(lrr.getCdna());
        }
        // if 2 reads, take the one with lowest de
        else if(this.longreads.size() == 2){
            LongreadRecord lrr1 = this.longreads.get(0).getBestRecord();
            LongreadRecord lrr2 = this.longreads.get(1).getBestRecord();
            this.consensus = (lrr1.getDe() < lrr2.getDe())?new String(lrr1.getCdna()):new String(lrr2.getCdna());
        }
        // if more do consensus calling
        else{            
            try {
                //int count=0;
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_reads.fa"));
                //os10 = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_reads10.fa"));
                Iterator<Longread> iterator2 = this.longreads.iterator();
                while(iterator2.hasNext()){
                    Longread lr = (Longread) iterator2.next();
                    LongreadRecord lrr = lr.getBestRecord();
                    os.writeBytes(">"+lr.getName()+"\n"+new String(lrr.getCdna())+"\n");
                    //if(count++<11)
                    //    os10.writeBytes(">"+lr.getName()+"\n"+new String(lrr.getCdna())+"\n");
                }
                os.close();
                //os10.close();

                // compute consensus with POA
                String[] commande = {"bash", "-c" , ""};
                commande[2] = POAPATH + "/poa -read_fasta "+TMPDIR+"/"+prefix+"_reads.fa -pir "+TMPDIR+"/"+prefix+".pir " + POAPATH + "/blosum80.mat -hb -best";
                //System.out.println(commande[2]);
                ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                // get CONSENS0 in fasta file .pir
                this.consensus = "";
                BufferedReader file = new BufferedReader(new FileReader(TMPDIR+"/"+prefix + ".pir"));
                line = file.readLine();
                while(! ">CONSENS0".equals(line)){ line = file.readLine(); }
                while(line != null && !"".equals(line)){ line = file.readLine(); this.consensus += line; }                
                file.close();
                this.consensus = this.consensus.replaceAll("-", "");
                this.consensus = this.consensus.replaceAll("\\n", "");
                this.consensus = this.consensus.replaceAll("null", "");
                
                // do we still need to do that after poa consensus ?
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_consensus.fa"));
                os.writeBytes(">"+this.barcode + this.umi+"\n"+this.consensus+"\n");
                os.close();

                commande[2] = MINIMAP2PATH + "/minimap2 --secondary=no -ax map-ont "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_reads.fa > "+TMPDIR+"/"+prefix+"_overlap.sam";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                commande[2] = RACONPATH + "/racon "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa > "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
                
                BufferedReader fichier = new BufferedReader(new FileReader(TMPDIR+"/"+prefix + "_corrected_consensus.fa"));
                line = fichier.readLine();
                this.consensus = fichier.readLine();
                fichier.close();
                
                commande[2] = "rm "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+".pir "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
           }
            catch(Exception e){ e.printStackTrace(); }
        }
        
        return ">" + this.getLabel() + "\n" + this.consensus + "\n";
    }
    
    /*
    public String call() throws Exception
    {
        double deMin = 1.0;
        String bestRead = "";
        List<DNASequence> lst = new ArrayList<DNASequence>();
        
        //System.out.println(this.getLabel());
        
        // should be an argument for consensus calling
        int nb_max_best_reads = 10;
        
        Collections.sort(this.longreads);
        Iterator<Longread> iterator = this.longreads.iterator();
        while (iterator.hasNext() && lst.size() < nb_max_best_reads) {
            Longread lr = (Longread) iterator.next();
            LongreadRecord lrr = lr.getBestRecord();
            String cdna = new String(lrr.getCdna());

            if (lrr.getDe() < deMin) {
                deMin = lrr.getDe();
                bestRead = cdna;
            }
            lst.add(new DNASequence(cdna));
        }

        try {
            if (lst.size() < 3) { this.consensus = bestRead; }
            else{
                //SimpleGapPenalty gapP = new SimpleGapPenalty((short) 5, (short) 2);
                //System.out.printf("Number of sequence to align\t" + lst.size() + "(" + getSizesOfSequence(lst) + ")\n");
                Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
                //System.out.printf("end of alignment\n");
                List<AlignedSequence<DNASequence, NucleotideCompound>> alignedSequence = profile.getAlignedSequences();
                //System.out.printf("Clustalw:%n%s%n", profile);
                Vector posConsVector = new Vector();
                int nSeqLength = profile.getLength();
                int nSeq = alignedSequence.size();

                for (int x = 0; x < nSeqLength; x++) { posConsVector.add(new int[5]); }

                for (int i = 0; i < nSeq; i++) {
                    String s = ((AlignedSequence) alignedSequence.get(i)).getSequenceAsString();
                    for (int x = 0; x < nSeqLength; x++) {
                        ((int[]) posConsVector.get(x))[getIndexForChar(s.charAt(x))]++;
                    }
                }
                //ConcurrencyTools.shutdown();
                for (int x = 0; x < nSeqLength; x++) {
                    this.consensus += getConsensusBase(((int[]) posConsVector.get(x)));
                }
                //System.out.println(this.consensus);
            }
        }
        catch (Exception e) { e.printStackTrace(); }
        
        if(this.consensus == null){ this.consensus = bestRead; }
        this.consensus = this.consensus.replaceAll("-", "");

        if (lst.size() > 2) {
            DataOutputStream os=null;
            try{
                String prefix = this.barcode + "" + this.umi;
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_consensus.fa"));
                os.writeBytes(">"+this.barcode + this.umi+"\n"+this.consensus+"\n");
                os.close();

                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_reads.fa"));
                Iterator<Longread> iterator2 = this.longreads.iterator();
                while(iterator2.hasNext()){
                    Longread lr = (Longread) iterator2.next();
                    LongreadRecord lrr = lr.getBestRecord();
                    os.writeBytes(">"+lr.getName()+"\n"+new String(lrr.getCdna())+"\n");
                }
                os.close();
                
                String[] commande = {"bash", "-c" , ""};
                commande[2] = MINIMAP2PATH + "/minimap2 --secondary=no -ax map-ont "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_reads.fa > "+TMPDIR+"/"+prefix+"_overlap.sam";
                ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();

                commande[2] = RACONPATH + "/racon "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa > "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
                        
                BufferedReader fichier = new BufferedReader(new FileReader(TMPDIR+"/"+prefix + "_corrected_consensus.fa"));
                String line = fichier.readLine();
                this.consensus = fichier.readLine();
                fichier.close();

                //commande[2] = "rm "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                //executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                //executeCmd.run();
           }
            catch(Exception e){ e.printStackTrace(); }
            finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        }
        
        if(this.consensus == null){ this.consensus = bestRead; }
        return ">" + this.getLabel() + "\n" + this.consensus + "\n";
    }

    public int getMedianExonBases()
    {
        int[] numArray = new int[this.longreads.size()];
        int i=0;
        Iterator<Longread> iterator = this.longreads.iterator();
        while (iterator.hasNext()) {
            Longread lr = (Longread) iterator.next();
            LongreadRecord lrr = lr.getAssociatedRecord();
            
            numArray[i++] = lrr.getExonBases();
        }
        
        Arrays.sort(numArray);
        int median;
        if (numArray.length % 2 == 0)
            median = ((int)numArray[numArray.length/2] + (int)numArray[numArray.length/2 - 1])/2;
        else
            median = (int) numArray[numArray.length/2];
        
        return median;
    }
    
    public double getMedianCna()
    {
        int[] numArray = new int[this.longreads.size()];
        int i=0;
        Iterator<Longread> iterator = this.longreads.iterator();
        while (iterator.hasNext()) {
            Longread lr = (Longread) iterator.next();
            LongreadRecord lrr = lr.getAssociatedRecord();
            
            numArray[i++] = lrr.getCdna().length();
        }
        
        Arrays.sort(numArray);
        double median;
        if (numArray.length % 2 == 0)
            median = ((int)numArray[numArray.length/2] + (int)numArray[numArray.length/2 - 1])/2;
        else
            median = (int) numArray[numArray.length/2];
        
        return median;
    }
    
    public String getExonBasesString()
    {
        String exonBasesList = "";
        Iterator<Longread> iterator = this.longreads.iterator();
        while (iterator.hasNext()) {
            Longread lr = (Longread) iterator.next();
            LongreadRecord lrr = lr.getAssociatedRecord();
            
            exonBasesList += new Integer(lrr.getExonBases()).toString() + ",";
        }
                
        return exonBasesList;
    }
    
    public String getCdnaString()
    {
        String str = "";
        Iterator<Longread> iterator = this.longreads.iterator();
        while (iterator.hasNext()) {
            Longread lr = (Longread) iterator.next();
            LongreadRecord lrr = lr.getAssociatedRecord();
            
            str += new Integer(lrr.getCdna().length()).toString() + ",";
        }
                
        return str;
    }
    */
}
