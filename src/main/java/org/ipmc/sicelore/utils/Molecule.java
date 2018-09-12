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

public class Molecule implements Callable<String> {

    private List<Longread> longreads;
    private HashSet<String> geneIds;
    private String barcode;
    private String umi;
    private String consensus = "";
    private String geneId = "undef";
    private String transcriptId = "undef";
    private int nbReads = 0;
    private int nbCleanReads = 0;
    private int nbConsensusReads = 0;
    private String pctEcartMedian = "";
    
    private final static HashMap<Character, Integer> encode;
    private final static char[] decode;

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
    
    public String getPctEcartMedian() {
        return pctEcartMedian;
    }
    
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
    
    public String getLabel() {
        return this.geneId + "|" + this.transcriptId + "|" + this.barcode + "|" + this.umi + "|" + this.nbReads + "|" + this.nbCleanReads + "|" + this.nbConsensusReads;
    }
    
    public int getNbReads() {
        return nbReads;
    }
    public int getNbCleanReads() {
        return nbCleanReads;
    }
    public int getNbConsensusReads() {
        return nbConsensusReads;
    }

    public void addLongread(Longread lr) {
        this.longreads.add(lr);
        if (lr.getIs_associated()) {
            this.geneIds.add(lr.getGeneId());
            this.geneId = lr.getGeneId();
        }

        // random peaking case of multigenes molecules
        if (geneIds.size() > 1) {
            int index = new Random().nextInt(geneIds.size());
            Iterator<String> iter = geneIds.iterator();
            for (int i = 0; i < index; i++) {
                iter.next();
            }
            this.geneId = (String) iter.next();
            //System.out.println(geneIds.toString() + "\t" + this.geneId);
        }
    }

    public void setIsoform(List<TranscriptRecord> transcripts, int DELTA, boolean SOFT)
    {
        Collections.sort(this.longreads);
        
        for (Longread lr : this.longreads) {
            //List<LongreadRecord> longreadrecords = lr.getLongreadrecords();
            LongreadRecord lrr = lr.getAssociatedRecord();
            
            if (SOFT){
                LinkedHashMap<String,Integer> map_result = create_rules(transcripts,DELTA);
                //for (LongreadRecord lrr : longreadrecords) {
                    TranscriptRecord transcriptrecord = asign_transcript(transcripts,lrr,DELTA,map_result);
                    if (transcriptrecord != null) {
                        if ("undef".equals(this.transcriptId)) {
                            this.transcriptId = transcriptrecord.getTranscriptId();
                            this.geneId = transcriptrecord.getGeneId();
                        }
                    }
                //}
            }
            else { // ############### standard attribution here ###############
                //for (LongreadRecord lrr : longreadrecords) {
                    List list = junctionsFromExons(lrr.getExons());
                    for (TranscriptRecord transcriptrecord : transcripts)
                    {                        
                        List list1 = junctionsFromExons(transcriptrecord.getExons());
                        if (map(list, list1, DELTA)) {
                            if ("undef".equals(this.transcriptId)) {
                                this.transcriptId = transcriptrecord.getTranscriptId();
                                this.geneId = transcriptrecord.getGeneId();
                            }
                            //else if(!this.transcriptId.equals(transcriptrecord.getTranscriptId())){
                            //    System.out.println("Other reads saying other transcriptId\t"+this.transcriptId+" --> "+transcriptrecord.getTranscriptId());
                            //}
                        }
                    }
                //}
            }
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

    public boolean map(List<int[]> lrr_exons, List<int[]> tr_exons, int DELTA)
    {
        boolean bool = true;
        if (tr_exons.size() == lrr_exons.size()) {
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
        //System.out.println(tr_exons.size() + "\t" + lrr_exons.size() + "-->"+bool);
        return bool;
    }

    public void removeChimeriaReads()
    {
        double min_pct = 0.90;
        double max_pct = 1.10;
        
        this.nbReads = this.longreads.size();
        Collections.sort(this.longreads);
        // this order to get the better a the last one
        Collections.reverse(this.longreads);
        
        double median = this.getMedianCna();
        
        Iterator<Longread> it = this.longreads.iterator();
        while (it.hasNext()){
            Longread lr = (Longread) it.next();
            LongreadRecord lrr = lr.getAssociatedRecord();
            
            double pct = ((lrr.getCdna().length() - median) * 100) / median;
            pctEcartMedian += String.format("%.02f", pct) + "\n";
            //System.out.println(pctEcartMedian);
            
            // require to be longer than 90% length of median ExonBases of molecules records
            //if(lrr.getExonBases() > (max_pct * medianExonBases) || lrr.getExonBases() < (min_pct * medianExonBases)){
            if(lrr.getCdna().length() > (max_pct * median) || lrr.getCdna().length() < (min_pct * median)){
                if(this.longreads.size() > 1){
                    it.remove();
                }
            }
        }
        this.nbCleanReads = this.longreads.size();
        Collections.sort(this.longreads);
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
        double dvMin = 1.0;
        String bestRead = "";
        //int maxExonBases=0;
        double min_pct = 0.90;
        double max_pct = 1.10;
        List<DNASequence> lst = new ArrayList<DNASequence>();

        // should be an argument for consensus calling
        int nb_max_best_reads = 10;

        try {
            Collections.sort(this.longreads);
            int medianExonBases = this.getMedianExonBases();
            
            Iterator<Longread> iterator = this.longreads.iterator();
            while (iterator.hasNext() && lst.size() < nb_max_best_reads) {
                Longread lr = (Longread) iterator.next();
                LongreadRecord lrr = lr.getAssociatedRecord();
                
                //if(maxExonBases == 0)
                //    maxExonBases = lrr.getExonBases();
                
                // require to be longer than 90% length of median ExonBases of molecules records
                if(lrr.getExonBases() > (min_pct * medianExonBases) && lrr.getExonBases() < (max_pct * medianExonBases))
                {
                    if (lrr.getDv() < dvMin) {
                        dvMin = lrr.getDv();
                        bestRead = lrr.getCdna();
                    }
                    lst.add(new DNASequence(lrr.getCdna()));
                }
            }

            //if("ACACCAAGTCGCGTGT".equals(this.barcode) && "CTGGGATTAC".equals(this.umi))
            //    System.out.println(this.longreads.size()+"|"+lst.size());
            this.nbConsensusReads = lst.size();
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
            }
        } catch (Exception e) { e.printStackTrace(); }
        
        this.consensus = this.consensus.replaceAll("-", "");
        
        /*
            Now need to polish the consensus with racon
        */
        if (lst.size() > 2) {
            DataOutputStream os=null;
            try{
                String prefix = this.barcode + "" + this.umi;
                os = new DataOutputStream(new FileOutputStream("/share/data/scratch/sicelore/"+prefix+"_consensus.fa"));
                os.writeBytes(">"+this.barcode + this.umi+"\n"+this.consensus+"\n");
                os.close();

                os = new DataOutputStream(new FileOutputStream("/share/data/scratch/sicelore/"+prefix+"_reads.fa"));
                Iterator<Longread> iterator = this.longreads.iterator();
                while(iterator.hasNext()){
                    Longread lr = (Longread) iterator.next();
                    LongreadRecord lrr = lr.getAssociatedRecord();
                    os.writeBytes(">"+lrr.getName()+"\n"+lrr.getCdna()+"\n");
                }
                os.close();

                String[] commande = {"bash", "-c" , ""};
                commande[2] = "/share/apps/local/minimap2/minimap2 --secondary=no -ax map-ont "+prefix+"_consensus.fa "+prefix+"_reads.fa > "+prefix+"_overlap.sam";
                ExecuteCmd executeCmd = new ExecuteCmd(commande, new String[0], "/share/data/scratch/sicelore/");
                executeCmd.run();

                commande[2] = "/share/apps/local/racon/bin/racon "+prefix+"_reads.fa "+prefix+"_overlap.sam "+prefix+"_consensus.fa > "+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], "/share/data/scratch/sicelore/");
                executeCmd.run();

                BufferedReader fichier = new BufferedReader(new FileReader("/share/data/scratch/sicelore/"+prefix+"_corrected_consensus.fa"));
                String line = fichier.readLine();
                //System.out.println(this.getLabel() + "\n" + this.consensus);
                this.consensus = fichier.readLine();
                //System.out.println(this.getLabel() + "\n" + this.consensus);
                fichier.close();

                commande[2] = "rm "+prefix+"_reads.fa "+prefix+"_overlap.sam "+prefix+"_consensus.fa "+prefix+"_corrected_consensus.fa";
                executeCmd = new ExecuteCmd(commande, new String[0], "/share/data/scratch/sicelore/");
                executeCmd.run();
           }
            catch(Exception e){ e.printStackTrace(); }
            finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        }
        
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
    
}
