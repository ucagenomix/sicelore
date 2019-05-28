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
import org.biojava.nbio.core.sequence.DNASequence;
import java.util.concurrent.*;
import org.apache.commons.lang3.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.alignment.template.AlignedSequence;

public class Molecule implements Callable<String>
{
    private List<Longread> longreads;
    private HashSet<String> geneIds;
    private String barcode;
    private String umi;
    private String consensus = "";
    private String geneId = "undef";
    private String transcriptId = "undef";
    
    protected static String TMPDIR;
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

    public void setStaticParams(String tmp, String racon, String minimap2){
	this.TMPDIR = tmp;
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

    public String getConsensus() {
        return this.consensus;
    }

    public String getGeneId() {
        return this.geneId;
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
        return this.geneId + "|" + this.transcriptId + "|" + this.barcode + "|" + this.umi + "|" + this.longreads.size();
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

    public void setIsoform(List<TranscriptRecord> transcripts, int DELTA, boolean SOFT)
    {
        List list1=null;
        
        HashMap<String, Integer> candidates = new HashMap<String, Integer>();
        
        for(Longread lr : this.longreads){
            List<LongreadRecord> records = lr.getLongreadrecords();
            //LongreadRecord lrr = lr.getBestRecord();

            for(LongreadRecord lrr : records){
                List list = junctionsFromExons(lrr.getExons());

                for(TranscriptRecord transcriptrecord : transcripts){                        
                    list1 = junctionsFromExons(transcriptrecord.getExons());

                    if(map(list, list1, DELTA, SOFT)) {
                        String key = transcriptrecord.getTranscriptId() + "|" + transcriptrecord.getGeneId();

                        // case SOFT pocess put nb exons as the value to get the more complex transcript
                        if(SOFT)
                            candidates.put(key, transcriptrecord.getExons().size());
                        else{
                            if(candidates.containsKey(key))
                                candidates.put(key, candidates.get(key) + 1);
                             else
                                candidates.put(key, 1);
                        }
                    }
                }
            }
        }
        
        if(candidates.size() > 0){
            HashSet<String> bestCandidates = new HashSet<String>();
            int maxValueInMap=(Collections.max(candidates.values()));  // This will return max value in the Hashmap
            for(Map.Entry<String, Integer> entry : candidates.entrySet()) {  // Iterate through hashmap
                if (entry.getValue() == maxValueInMap) { bestCandidates.add(entry.getKey()); }
            }
            int index = new Random().nextInt(bestCandidates.size());
            Iterator<String> iter = bestCandidates.iterator();
            for (int i = 0; i < index; i++) { iter.next(); }
            String g = (String) iter.next();
            this.transcriptId = g.split("\\|")[0];
            this.geneId = g.split("\\|")[1];
        }
        else{
            int index = new Random().nextInt(geneIds.size());
            Iterator<String> iter = geneIds.iterator();
            for (int i = 0; i < index; i++) { iter.next(); }
            this.transcriptId = "undef";
            this.geneId = (String) iter.next();
        }
        
        // only 1 mono-exonic transcript in the model --> we set the isoform
        if(transcripts.size() == 1 && list1.size() == 0){
           this.transcriptId = transcripts.get(0).getTranscriptId();
           this.geneId = transcripts.get(0).getGeneId();
        }
    }
    
    public LinkedHashMap<String,Integer> create_rules(List<TranscriptRecord> transcriptrecord,int DELTA)
    {
        LinkedHashMap<String, Integer> map_result = new LinkedHashMap<>();
        List<List<int[]>> list_transcript = new ArrayList<>();
        for (int cpt = 0; cpt < transcriptrecord.size(); cpt++) {
            list_transcript.add(junctionsFromExons(transcriptrecord.get(cpt).getExons()));
        }
        for (int i = 0; i < list_transcript.size() - 1; i++) {
            for (int k = 0; k < list_transcript.get(i).size(); k++) {
                boolean result = false;
                int[] element_test = list_transcript.get(i).get(k);
                for (int ii = 0; ii < list_transcript.size(); ii++) {
                    for (int jj = 0; jj < list_transcript.get(ii).size(); jj++) {
                        if ((Math.abs(element_test[0] - list_transcript.get(ii).get(jj)[0]) <= DELTA) && (Math.abs(element_test[1] - list_transcript.get(ii).get(jj)[1]) <= DELTA) && i != ii) {
                            result = true;
                            list_transcript.get(ii).remove(jj);
                            jj--;
                        }
                    }
                }
                if (result) {
                    list_transcript.get(i).remove(k);
                }
            }
        }
        for (int iii = 0; iii < list_transcript.size(); iii++) {
            for (int j = 0; j < list_transcript.get(iii).size(); j++) {
                String key = list_transcript.get(iii).get(j)[0] + " " + list_transcript.get(iii).get(j)[1];
                map_result.put(key, iii);
            }
        }
        return map_result;
    }
    
    public TranscriptRecord asign_transcript(List<TranscriptRecord> transcriptrecord, LongreadRecord lrr, int DELTA, LinkedHashMap<String,Integer> map_result)
    {
        List<int[]> list_lrr = junctionsFromExons(lrr.getExons());
        for (int j = 0; j < list_lrr.size(); j++) {
            Set<String> keys = map_result.keySet();
            for (String key : keys) {
                String[] split_result = key.split(" ");
                if ((Math.abs(list_lrr.get(j)[0] - Integer.valueOf(split_result[0])) <= DELTA) && (Math.abs(list_lrr.get(j)[1] - Integer.valueOf(split_result[1])) <= DELTA)) {
                    return transcriptrecord.get(map_result.get(key));
                }
            }
        }
        return null;
    }

    public List<int[]> junctionsFromExons(List<int[]> exons) {
        ArrayList localArrayList = new ArrayList();

        for (int i = 1; i < exons.size(); i++) {
            int j = ((int[]) exons.get(i - 1))[1];
            int k = ((int[]) exons.get(i))[0];
            localArrayList.add(new int[]{j, k});
        }

        return localArrayList;
    }

    public void print(List<int[]> junctions) {
        for (int i = 0; i < junctions.size(); i++) {
            System.out.println(junctions.get(i)[0] + "-" + junctions.get(i)[1]);
        }
    }

    public boolean map(List<int[]> lrr_exons, List<int[]> tr_exons, int DELTA, boolean SOFT)
    {
        boolean bool = true;
	if(SOFT){
            if(tr_exons.size() <= lrr_exons.size()){
	   	for(int i=0; i<tr_exons.size(); i++){
                    if(!isIn((int[])tr_exons.get(i), lrr_exons, DELTA))
	        	bool = false;
	       }
            }
	    else{ bool = false; }
	}
        else{
            if(tr_exons.size() > 0 && tr_exons.size() == lrr_exons.size()) {
                for (int i = 0; i < tr_exons.size(); i++) {
                    if (!isIn((int[]) tr_exons.get(i), lrr_exons, DELTA)) {
                        bool = false;
                        //System.out.println(tr_exons.get(i)[0] +"-" + tr_exons.get(i)[1] + " not in read !");
                    }
                }
            } else {
                bool = false;
                //System.out.println("not same exons size");
            }
        }
        
        //System.out.println(tr_exons.size() + "\t" + lrr_exons.size() + "-->"+bool);
        return bool;
    }

    public boolean isAlreadyIn(int[] paramArrayOfInt, List<int[]> paramList) {
        boolean bool = false;
        for (int[] arrayOfInt : paramList) {
            if ((arrayOfInt[0] == paramArrayOfInt[0]) && (arrayOfInt[1] == paramArrayOfInt[1])) {
                bool = true;
            }
        }
        return bool;
    }

    public boolean isIn(int[] paramArrayOfInt, List<int[]> paramList, int paramInt) {
        boolean bool = false;
        for (int[] arrayOfInt : paramList) {
            if ((Math.abs(arrayOfInt[0] - paramArrayOfInt[0]) <= paramInt) && (Math.abs(arrayOfInt[1] - paramArrayOfInt[1]) <= paramInt)) {
                bool = true;
            }
        }
        return bool;
    }

    private static int[] toIntArray(String paramString) throws NumberFormatException {
        paramString = StringUtils.stripEnd(paramString, ",");
        String[] arrayOfString = paramString.split(",");
        int[] arrayOfInt = new int[arrayOfString.length];

        for (int i = 0; i < arrayOfString.length; i++) {
            arrayOfInt[i] = Integer.valueOf(arrayOfString[i]).intValue();
        }
        return arrayOfInt;
    }

    public String getSizesOfSequence(List<DNASequence> lst) {
        String t = "";
        for (int i = 0; i < lst.size(); i++) {
            t += new Integer(((DNASequence) lst.get(i)).getSequenceAsString().length()) + ",";
        }
        return t;
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
        double deMin = 1.0;
        String bestRead = "";
        List<DNASequence> lst = new ArrayList<DNASequence>();
        
        //System.out.println(this.getLabel());
        
        // should be an argument for consensus calling
        int nb_max_best_reads = 10;

        try {
            Collections.sort(this.longreads);
            
            Iterator<Longread> iterator = this.longreads.iterator();
            while (iterator.hasNext() && lst.size() < nb_max_best_reads) {
                Longread lr = (Longread) iterator.next();
                LongreadRecord lrr = lr.getBestRecord();
                
                String cdna = new String(lrr.getCdna());
                //System.out.println(cdna);
                
                if (lrr.getDe() < deMin) {
                    deMin = lrr.getDe();
                    bestRead = cdna;
                }
                lst.add(new DNASequence(cdna));
            }

            //if("ACACCAAGTCGCGTGT".equals(this.barcode) && "CTGGGATTAC".equals(this.umi))
            //    System.out.println(this.longreads.size()+"|"+lst.size());
            //this.nbConsensusReads = lst.size();
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
        catch (Exception e) { 
            e.printStackTrace(); 
            //System.out.println(this.getLabel()+"\n"+lst);
            this.consensus = bestRead;
        }
        
        this.consensus = this.consensus.replaceAll("-", "");
        
        /*
            Now need to polish the consensus with racon
        */
        if (lst.size() > 2) {
            DataOutputStream os=null;
            try{
                String prefix = this.barcode + "" + this.umi;
                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_consensus.fa"));
                os.writeBytes(">"+this.barcode + this.umi+"\n"+this.consensus+"\n");
                os.close();

                os = new DataOutputStream(new FileOutputStream(TMPDIR+"/"+prefix+"_reads.fa"));
                Iterator<Longread> iterator = this.longreads.iterator();
                while(iterator.hasNext()){
                    Longread lr = (Longread) iterator.next();
                    LongreadRecord lrr = lr.getBestRecord();
                    os.writeBytes(">"+lr.getName()+"\n"+new String(lrr.getCdna())+"\n");
                }
                os.close();
                
                String[] commande = {"bash", "-c" , ""};
                commande[2] = MINIMAP2PATH + "/minimap2 --secondary=no -ax map-ont "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_reads.fa > "+TMPDIR+"/"+prefix+"_overlap.sam";
                ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
                
                //System.out.println(System.getenv("PATH"));
                //executeCmd.getError();

                commande[2] = RACONPATH + "/racon "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa > "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
                
                        
                BufferedReader fichier = new BufferedReader(new FileReader(TMPDIR+"/"+prefix + "_corrected_consensus.fa"));
                String line = fichier.readLine();
               
                //System.out.println("avant\t"+this.getLabel() + "\n" + this.consensus);
                this.consensus = fichier.readLine();
                //System.out.println("apres\t"+this.getLabel() + "\n" + this.consensus);
                
                fichier.close();

                commande[2] = "rm "+TMPDIR+"/"+prefix+"_reads.fa "+TMPDIR+"/"+prefix+"_overlap.sam "+TMPDIR+"/"+prefix+"_consensus.fa "+TMPDIR+"/"+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], TMPDIR);
                executeCmd.run();
           }
            catch(Exception e){ e.printStackTrace(); }
            finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        }
        
        return ">" + this.getLabel() + "\n" + this.consensus + "\n";
    }
    /*
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
