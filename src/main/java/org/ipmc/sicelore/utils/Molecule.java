package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 *  
 */
import java.util.*;

public class Molecule
{
    private List<Longread> longreads;
    private HashSet<String> geneIds;
    private HashSet<Junction> junctionSet;
    private String barcode;
    private String umi;
    private int rn;
    private Float pctId;
    private byte[] consensus;
    private byte[] consensusQV;
    private int consensusLength;
    private String snpPhredScore = "";
    
    private String geneId = "undef";
    private String transcriptId = "undef";
    private int supporting_reads=0;

    public Molecule() {}
    
    public Molecule(String barcode, String umi, int rn) {
        this.longreads = new ArrayList<Longread>();
        this.geneIds = new HashSet<String>();
        this.junctionSet = new HashSet<Junction>();
        this.barcode = barcode;
        this.umi = umi;
        this.rn = rn;
    }
    
    public Molecule(String barcode, String umi, String consensus, int rn) {
        //this.longreads = new ArrayList<Longread>();
        //this.geneIds = new HashSet<String>();
        this.junctionSet = new HashSet<Junction>();
        this.barcode = barcode;
        this.consensus = consensus.getBytes();;
        this.umi = umi;
        this.rn=rn;
    }
    
    public Molecule(String barcode, String umi, String consensus, String consensusQV, int rn) {
        //this.longreads = new ArrayList<Longread>();
        //this.geneIds = new HashSet<String>();
        this.junctionSet = new HashSet<Junction>();
        this.barcode = barcode;
        this.consensusLength = consensus.length();
        this.consensus = consensus.getBytes();
        this.consensusQV = consensusQV.getBytes();
        this.umi = umi;
        this.rn=rn;
    }
    
    public List<Longread> getLongreads() { return longreads; }
    public HashSet<String> getGeneIds() { return geneIds; }
    public HashSet<Junction> getJunctionSet() { return this.junctionSet; }
    
    public String getBarcode() { return this.barcode; }
    public String getUmi() { return this.umi; }
    public Float getPctId() { return this.pctId; }
    
    public int getConsensusLength() { return this.consensusLength; }
    public void setConsensusLength(int consensusLength) { this.consensusLength = consensusLength; }

    public byte[] getConsensus() { return this.consensus; }
    public void setConsensus(byte[] consensus) { this.consensus = consensus; }
    
    public byte[] getConsensusQV() { return this.consensusQV; }
    public void setConsensusQV(byte[] consensusQV) { this.consensusQV = consensusQV; }
    
    public String getSnpPhredScore() { return this.snpPhredScore; }
    public void setSnpPhredScore(String snpPhredScore) { this.snpPhredScore = snpPhredScore; }
    
    public int getRn() { return this.rn; }
    public void setRn(int rn) { this.rn = rn; }

    public String getGeneId() { return this.geneId; }
    public void setGeneId(String geneId) { this.geneId = geneId; }
    
    public String getTranscriptId() { return this.transcriptId; }
    public void setTranscriptId(String transcriptId) { this.transcriptId=transcriptId; }
    
    public int getSupporting_reads() { return this.supporting_reads; }
    public void setSupporting_reads(int r) { this.supporting_reads=r; }
    
    public String[] getGeneIdsArray() {
        String[] array = new String[this.geneIds.size()];
        this.geneIds.toArray(array);
        return array;
    }

    public String getLabel() {
        return this.barcode + "-" + this.umi + "-" + this.longreads.size();
    }
    
    public int getNumberOfReads() {
        if(this.rn > 1)
            return this.rn;
        else
            return this.longreads.size();
    }
    
    public String toString()
    {
        String str = "cell="+this.barcode+",umi="+ this.umi + ", reads=" + this.longreads.size() + ", " + this.geneIds.toString();
        Iterator<Longread> it = this.longreads.iterator();
        while(it.hasNext()){
            Longread lr = (Longread) it.next();
            List<LongreadRecord> lrr = lr.getLongreadrecords();
            str += "\n\t" + lr.getName() + ", SAMrecords="+lrr.size();
        }
        
        return str;
    }
    
    public void addLongread(Longread lr)
    {
        this.longreads.add(lr);
        this.pctId = new Float(1) - lr.getLongreadrecords().get(0).getDe();
        
        Iterator<String> iterator = lr.getGeneIds().iterator();
        while (iterator.hasNext())
            this.geneIds.add((String)iterator.next());
    }
    
    public void addJunction(Junction j)
    {
        this.junctionSet.add(j);
    }
}
