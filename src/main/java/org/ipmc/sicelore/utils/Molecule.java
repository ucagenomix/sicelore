package org.ipmc.sicelore.utils;

import java.util.*;
import org.ipmc.sicelore.utils.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.*;
import java.util.concurrent.*;
import org.apache.commons.lang3.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.alignment.SimpleAlignedSequence;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.alignment.SimpleGapPenalty;

public class Molecule implements Callable<String> {

    private List<Longread> longreads;
    private HashSet<String> geneIds;
    private String barcode;
    private String umi;
    private String consensus = "";
    private String geneId = "undef";
    private String transcriptId = "undef";

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
        return this.geneId + "|" + this.transcriptId + "|" + this.barcode + "|" + this.umi + "|" + this.longreads.size();
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

    public void setIsoform(List<TranscriptRecord> transcripts, int DELTA, boolean SOFT) {
        for (Longread lr : this.longreads) {
            List<LongreadRecord> longreadrecords = lr.getLongreadrecords();
            if (SOFT) {
                LinkedHashMap<String,Integer> map_result = create_rules(transcripts,DELTA);
                for (LongreadRecord lrr : longreadrecords) {
                    TranscriptRecord transcriptrecord = asign_transcript(transcripts,lrr,DELTA,map_result);
                    if (transcriptrecord != null) {
                        if ("undef".equals(this.transcriptId)) {
                            this.transcriptId = transcriptrecord.getTranscriptId();
                            this.geneId = transcriptrecord.getGeneId();
                        }
                    }
                }
            } else {
                for (LongreadRecord lrr : longreadrecords) {
                    List list = junctionsFromExons(lrr.getExons());
                    
                    //System.out.println(lrr.getName());
                    //this.print(lrr.getExons());
                    //this.print(list);

                    for (TranscriptRecord transcriptrecord : transcripts) {
                        
                        List list1 = junctionsFromExons(transcriptrecord.getExons());
                        //System.out.println(transcriptrecord.getTranscriptId());
                        int start_abs = Math.abs(transcriptrecord.getTxStart() - lrr.getTxStart());
                        int  end_abs = Math.abs(transcriptrecord.getTxEnd() - lrr.getTxEnd());
                        
                        if (map(list, list1, DELTA, SOFT) && start_abs < DELTA && end_abs < DELTA) {
                            if ("undef".equals(this.transcriptId)) {
                                this.transcriptId = transcriptrecord.getTranscriptId();
                                this.geneId = transcriptrecord.getGeneId();
                            }
                        }
                    }
                }
            }
        }
    }
    public LinkedHashMap<String,Integer> create_rules(List<TranscriptRecord> transcriptrecord,int DELTA) {
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
    public TranscriptRecord asign_transcript(List<TranscriptRecord> transcriptrecord, LongreadRecord lrr, int DELTA, LinkedHashMap<String,Integer> map_result) {
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

    public boolean map(List<int[]> lrr_exons, List<int[]> tr_exons, int DELTA, boolean SOFT) {
        boolean bool = true;
        
        if (SOFT) {
            if (tr_exons.size() <= lrr_exons.size()) {
                for (int i = 0; i < tr_exons.size(); i++) {
                    if (!isIn((int[]) tr_exons.get(i), lrr_exons, DELTA)) {
                        bool = false;
                    }
                }
            } else {
                bool = false;
            }
        } else {
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
        }

        //System.out.println(tr_exons.size() + "\t" + lrr_exons.size() + "-->"+bool);
        return bool;
    }

    public void cleanYourChimeriaReads() {

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

    public String call() throws Exception {
        int nb = 0;
        double dvMin = 1.0;
        String bestRead = "";

        // should be an argument for consensus calling
        int nb_max_best_reads = 10;

        try {
            int nn = 0;
            String str = "";
            Collections.sort(this.longreads);
            List<DNASequence> lst = new ArrayList<DNASequence>();

            Iterator<Longread> iterator = this.longreads.iterator();
            while (iterator.hasNext() && nn < nb_max_best_reads) {
                Longread lr = (Longread) iterator.next();
                LongreadRecord lrr = lr.getAssociatedRecord();

                // in case we don't have an associated read, should not be possible 
                // now that we remove unassociated molecule at new MoleculeDataset(bam) step
                if (lrr != null) {
                    if (lrr.getDv() < dvMin) {
                        dvMin = lrr.getDv();
                        bestRead = lrr.getCdna();
                    }
                    lst.add(new DNASequence(lrr.getCdna()));
                    nn++;
                    str += "," + lrr.getDv();
                }
            }

            if (lst.size() < 3) {
                this.consensus = bestRead;
            } else {

                //if(lst.size() > 3){
                //  System.out.println(lst.size() + "\t" + str);
                //}
                //SimpleGapPenalty gapP = new SimpleGapPenalty((short) 5, (short) 2);
                //System.out.printf("Number of sequence to align\t" + lst.size() + "(" + getSizesOfSequence(lst) + ")\n");
                Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
                //System.out.printf("end of alignment\n");
                List<AlignedSequence<DNASequence, NucleotideCompound>> alignedSequence = profile.getAlignedSequences();
                //System.out.printf("Clustalw:%n%s%n", profile);
                Vector posConsVector = new Vector();
                int nSeqLength = profile.getLength();
                int nSeq = alignedSequence.size();

                for (int x = 0; x < nSeqLength; x++) {
                    posConsVector.add(new int[5]);
                }

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
        } catch (Exception e) {
            e.printStackTrace();
        }

        this.consensus = this.consensus.replaceAll("-", "");
        //System.out.println(">"+this.getLabel()+"\n"+this.consensus);
        return ">" + this.getLabel() + "\n" + this.consensus + "\n";
    }
}
