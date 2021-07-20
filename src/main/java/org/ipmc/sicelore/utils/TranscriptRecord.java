package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */ 
import htsjdk.tribble.annotation.Strand;
import java.util.*;
import org.apache.commons.lang3.StringUtils;

public class TranscriptRecord implements Comparable<TranscriptRecord> {

    private String transcriptId;
    private String geneId;
    private String chrom;
    private Strand strand;
    private int txStart;
    private int txEnd;
    private int cdsStart;
    private int cdsEnd;
    private int exonCount;
    private int exonBases;
    private int cdsExonBases;
    private int[] exonStarts;
    private int[] exonEnds;
    private int[] exonFrames;
    
    private List<int[]> codingExons;
    private List<int[]> exons;
    private List<Junction> junctions;
    private double rpkm;
    
    // collapseModel pipeline requirements
    private List<LongreadRecord> evidenceList;
    private String categorie = "undef";
    private String subcategorie = "undef2";
    private List<int[]> novelJunctions;
    boolean is_novel = false;
    boolean is_known = true;
    private byte[] representative;
    public int nbUmis;
    public int nbCells;
    
    // annotateModel pipeline requirements
    public boolean is_valid = false;
    public boolean is_valid_cage = false;
    public boolean is_valid_polya = false;
    public boolean is_valid_junction = false;
    public int dist_cage = 0;
    public int dist_polya = 0;
    public int junctionReads = 0;
    
    // consensus pipeline requirements
    //public Consensus consensus;
    
    public TranscriptRecord() {
        evidenceList = new ArrayList<LongreadRecord>();
        novelJunctions = new ArrayList<int[]>();
    }
    
    public TranscriptRecord(String geneId, String transcriptId){
        evidenceList = new ArrayList<LongreadRecord>();
        novelJunctions = new ArrayList<int[]>();
        this.geneId = geneId;
        this.transcriptId = transcriptId;
    }

    public boolean equals (Object o){
        TranscriptRecord x = (TranscriptRecord) o;
        
        if(!this.geneId.equals(x.geneId)) return false;
        if(!this.transcriptId.equals(x.transcriptId)) return false;
        
        return true;
    }
    
    public int hashCode(){
      int hash = 5;
      hash = 59 * hash + this.geneId.hashCode();
      hash = 59 * hash + this.transcriptId.hashCode();
      return hash;
   }
    
    public int compareTo(TranscriptRecord tr) {
        // max exons to min
        return ((TranscriptRecord) tr).getExons().size() - this.exons.size();
        // min exons to max
        //return this.exonCount - ((TranscriptRecord)tr).getExonCount();
    }

    public static TranscriptRecord fromRefFlat(String[] fields) throws GTFParseException
    {
        if (fields.length < 11) {
            throw new RuntimeException(
                    "Invalid RefGene file. records should have at least 11 fields but found only: "
                    + fields.length);
        }

        TranscriptRecord record = new TranscriptRecord();
        record.junctions = new ArrayList();
        record.geneId = fields[0];
        record.transcriptId = fields[1];
        record.chrom = fields[2];
        record.strand = Strand.toStrand(fields[3]);

        try {
            record.txStart = Integer.valueOf(fields[4]);
            record.txEnd = Integer.valueOf(fields[5]);
            record.cdsStart = Integer.valueOf(fields[6]);
            record.cdsEnd = Integer.valueOf(fields[7]);
            record.exonCount = Integer.valueOf(fields[8]);
            record.exonStarts = TranscriptRecord.toIntArray(fields[9]);
            record.exonEnds = TranscriptRecord.toIntArray(fields[10]);
            //record.exonFrames = TranscriptRecord.toIntArray(fields[15]);
            record.exonBases = 0;
            record.exons = new ArrayList<int[]>();
            record.codingExons = new ArrayList<int[]>();

            for (int i = 0; i < record.exonStarts.length; i++) {
                int start = record.exonStarts[i];
                int end = record.exonEnds[i];

                record.exonBases += end - start;
                record.exons.add(new int[]{start+1, end});

                // Compute coding exons
                if (start > record.cdsEnd) {
                    continue;
                }
                if (end < record.cdsStart) {
                    continue;
                }

                if (start >= record.cdsStart && end <= record.cdsEnd) {
                    record.codingExons.add(new int[]{start, end});
                    record.cdsExonBases += end - start;
                } else if (start <= record.cdsStart && record.cdsStart <= end && end <= record.cdsEnd) {
                    record.codingExons.add(new int[]{record.cdsStart, end});
                    record.cdsExonBases += end - record.cdsStart;
                } else if (start >= record.cdsStart && record.cdsStart <= end && end >= record.cdsEnd) {
                    record.codingExons.add(new int[]{start, record.cdsEnd});
                    record.cdsExonBases += record.cdsEnd - start;
                } else if (start < record.cdsStart && end > record.cdsEnd) {
                    record.codingExons.add(new int[]{record.cdsStart, record.cdsEnd});
                    record.cdsExonBases += record.cdsEnd - record.cdsStart;
                }
            }
            
            for (int i = 1; i < record.exons.size(); i++) {
                int j = ((int[]) record.exons.get(i - 1))[1];
                int k = ((int[]) record.exons.get(i))[0];
                record.junctions.add(new Junction(j, k));
            }

        } catch (NumberFormatException e) {
            throw new GTFParseException(
                    "Invalid RefGene file. Can't parse integer value: ", e);
        }

        record.rpkm = 0.0;

        return record;
    }
    
    public List<LongreadRecord> getEvidenceList(){ return this.evidenceList; }
    
    public void add(LongreadRecord lrr) {
        evidenceList.add(lrr);
    }
    
    public String getGeneId() { return geneId; }
    public String getTranscriptId() { return transcriptId; }
    public String getChrom() { return chrom; }
    public Strand getStrand() { return strand; }
    public int getTxStart() { return txStart; }
    public int getTxEnd() { return txEnd; }
    public int getCdsStart() { return cdsStart; }
    public int getCdsEnd() { return cdsEnd; }
    public int getExonCount() { return exonCount; }
    public int getExonBases() { return exonBases; }
    public int getCdsExonBases() { return cdsExonBases; }
    public int[] getExonStarts() { return exonStarts; }
    public int[] getExonEnds() { return exonEnds; }
    public int[] getExonFrames() { return exonFrames; }

    public int getNbUmis() { return nbUmis; }
    public int getNbCells() { return nbCells; }
    
    public void setJunctionReads(int junctionReads) { this.junctionReads = junctionReads; }
    public int getJunctionReads() { return this.junctionReads; }
    
    public void setDist_cage(int dist_cage) { this.dist_cage = dist_cage; }
    public int getDist_cage() { return this.dist_cage; }
    
    public void setIs_valid_cage(boolean is_valid_cage) { this.is_valid_cage = is_valid_cage; }
    public boolean getIs_valid_cage() { return this.is_valid_cage; }
    
    public void setDist_polya(int dist_polya) { this.dist_polya = dist_polya; }
    public int getDist_polya() { return this.dist_polya; }
    
    public void setIs_valid_polya(boolean is_valid_polya) { this.is_valid_polya = is_valid_polya; }
    public boolean getIs_valid_polya() { return this.is_valid_polya; }
    
    public void setIs_valid_junction(boolean is_valid_junction) { this.is_valid_junction = is_valid_junction; }
    public boolean getIs_valid_junction() { return this.is_valid_junction; }
    
    public void setIs_valid(boolean is_valid) { this.is_valid = is_valid; }
    public boolean getIs_valid() { return this.is_valid; }
    
    public void setRPKM(double rpkm) { this.rpkm = rpkm; }
    public double getRPKM() { return this.rpkm; }
    
    public void setCategorie(String categorie) { this.categorie=categorie; }
    public String getCategorie() { return categorie; }
    
    public void setSubcategorie(String subcategorie) { this.subcategorie=subcategorie; }
    public String getSubcategorie() { return subcategorie; }
    
    public List<int[]> getNovelJunctions() { return this.novelJunctions; }
    public void setNovelJunctions(List<int[]> novelJunctions) { this.novelJunctions = novelJunctions; }
    
    public void setIs_novel(boolean is_novel) { this.is_novel=is_novel; }
    public boolean getIs_novel() { return this.is_novel; }
    
    public void setIs_known(boolean is_known) { this.is_known=is_known; }
    public boolean getIs_known() { return this.is_known; }

    public List<int[]> getCodingExons() { return this.codingExons; }
    public List<int[]> getExons() { return this.exons; }
    public void setExons(List<int[]> e) { this.exons = e; }
    public List<Junction> getJunctions() { return this.junctions; }
    
    public List<int[]> getExons(boolean cdsExonsOnly) {
        return cdsExonsOnly ? this.codingExons : this.exons;
    }

    //public Consensus getConsensus() { return this.consensus; }
 
    public void addNovelJunction(int[] j) {
        this.novelJunctions.add(j);
    }
    
    public String toString() {
        return "[ transcriptId=" + transcriptId + ", exonCount=" + exonCount + "]";
    }
    
    public String printLegendTxt()
    {
        // standard annotation
        String str = "geneId\ttranscriptId\tchrom\tstrand\ttxStart\ttxEnd\texons\tUMIs\tCells";
        // classification
        str += "\tcategorie\tsubcategorie\tnovelJunctions";
        //validation
        str += "\tnovelJunctions_reads\tis_valid_allNovelJunctions\tdist_cage\tis_valid_cage\tdist_polya\tis_valid_polya\tis_valid\n";
        
        return str;
    }
    
    public String printTxt()
    {    
        // standard annotation
        String str = geneId+"\t"+transcriptId+"\t"+chrom+"\t"+strand+"\t"+txStart+"\t"+txEnd+"\t"+exons.size()+"\t"+nbUmis+"\t"+nbCells;
        // classification
        str += "\t"+categorie+"\t"+subcategorie+"\t"+this.getNovelJunctionsString();
        //validation
        str += "\t"+junctionReads+"\t"+is_valid_junction+"\t"+dist_cage+"\t"+is_valid_cage+"\t"+dist_polya+"\t"+is_valid_polya+"\t"+is_valid+"\n";
        
        return str;
    }
    
    public String printRefflat()
    {    
        // RP23-333I7.1    ENSMUST00000194643.1    1       -       3905738 3986215 3986215 3986215 3       3905738,3985159,3986146,        3906134,3985351,3986215, 
        return geneId+"\t"+transcriptId+"\t"+chrom+"\t"+strand+"\t"+txStart+"\t"+txEnd+"\t"+cdsStart+"\t"+cdsEnd+"\t"+exons.size()+"\t"+this.getExonStartsString()+"\t"+this.getExonEndsString()+"\n";
    }
    
    public String getExonStartsString()
    {
        String str = "";
        for(int i=0;i<this.exons.size();i++)
            str+= (this.exons.get(i)[0]-1) + ",";
        return str;
    }
    
    public String getExonEndsString()
    {
        String str = "";
        for(int i=0;i<this.exons.size();i++)
            str+= this.exons.get(i)[1] + ",";
        return str;
    }
    
    public String getNovelJunctionsString(){
        String str = "";
        if(novelJunctions.size()>0){
            for (int[] a : novelJunctions)
                str += a[0] + "-" + a[1] + ",";
            str=str.substring(0,str.length()-1);
        }
        else
            str += "-";
        
        return str;
    }
    
    public String printGff()
    {
/*        
1       ONT     transcript      4774308 4785726 .       -       .       gene_id "PB.135"; transcript_id "ONT.135.19"; ensembl_id "NA"; 
1       ONT     exon            4774308 4774516 .       -       .       gene_id "PB.135"; transcript_id "ONT.135.19"; ensembl_id "NA"; category "full-splice_match"; subcategory "multi-exon"; color "#ff1aff";
1       ONT     exon            4777525 4777648 .       -       .       gene_id "PB.135"; transcript_id "ONT.135.19"; ensembl_id "NA"; category "full-splice_match"; subcategory "multi-exon"; color "#ff1aff";
1       ONT     exon            4782568 4782733 .       -       .       gene_id "PB.135"; transcript_id "ONT.135.19"; ensembl_id "NA"; category "full-splice_match"; subcategory "multi-exon"; color "#ff1aff";
1       ONT     exon            4783951 4784105 .       -       .       gene_id "PB.135"; transcript_id "ONT.135.19"; ensembl_id "NA"; category "full-splice_match"; subcategory "multi-exon"; color "#ff1aff";
1       ONT     exon            4785573 4785726 .       -       .       gene_id "PB.135"; transcript_id "ONT.135.19"; ensembl_id "NA"; category "full-splice_match"; subcategory "multi-exon"; color "#ff1aff";
*/
        String str = chrom + "\tsicelore\ttranscript\t" + txStart + "\t" + txEnd + "\t.\t" + strand + "\t.\tgene_id \"" + geneId + "\"; transcript_id \""+ transcriptId +"\"; category \"" + categorie + "\"; subcategory \""+ subcategorie +"\"; UMIs \"" + nbUmis + "\"; Cells \"" + nbCells + "\"; ";
        str += "novelJunctions \"" + this.getNovelJunctionsString() +"\"; ";
        str += "supportingReads \"" + this.junctionReads +"\"; ";
        str += "CAGEdist \"" + dist_cage +"\"; ";
        str += "POLYAdist \"" + dist_polya +"\"; ";
        str += "color \""+ this.getColor(subcategorie) + "\";\n";
        for(int i=0;i<exons.size(); i++)
            str += chrom + "\tsicelore\texon\t" + exons.get(i)[0] + "\t" + exons.get(i)[1] + "\t.\t" + strand + "\t.\tgene_id \"" + geneId + "\"; transcript_id \""+ transcriptId +"\";\n";
        
        return str;
    }
    
    protected static String getColor(String cat)
    {
        String col = "#000000";
        if (cat.equals("gencode")) col = "#014e8e";
        if (cat.equals("combination_of_known_junctions")) col = "#9dd122";
        if (cat.equals("combination_of_known_splicesites")) col = "#c594e1";
        if (cat.equals("at_least_one_novel_splicesite")) col = "#e65802";
        
        return col;
    }
    
    public String printFas() {
        return ">" + geneId + "|" + transcriptId + "\n" + new String(representative) +"\n";
    }
    
    public String printExonsList() {
        String str="";
        for (int i=0; i<exons.size(); i++) 
            str += new String(exons.get(i)[0]+"-"+exons.get(i)[1]+",");
        
        str = str.substring(0,str.length()-1);
        return str;
    }
    
    // init categorie and subcategorie
    // set representative sequence as the longest one
    // correct start first and last exons to min txStart and max txEnd
    // observed among all lrr showing the exonic structure of the isoform
    public void initialize()
    {
        this.representative = "A".getBytes();
        int minStart = Integer.MAX_VALUE;
        int maxEnd = Integer.MIN_VALUE;
        
        //System.out.println(transcriptId+"\t"+evidenceList.size());
        
        for (int i=0; i<evidenceList.size(); i++){
            LongreadRecord lrr = evidenceList.get(i);
            strand = lrr.getStrand();
            chrom = lrr.getChrom();
            if(minStart > lrr.getTxStart()) { minStart = lrr.getTxStart(); }
            if(maxEnd < lrr.getTxEnd()) { maxEnd = lrr.getTxEnd(); }
            
            // keep the longest lrr cDNA as representative
            // next step would be to compute consensus with POA just as molecule, need an object doing that
            if(representative.length < lrr.getCdna().length)
                representative = lrr.getCdna();
        }
        
        if(is_novel){
            categorie="undef";
            subcategorie="undef2";
            exons.get(0)[0]=minStart;
            exons.get(exons.size()-1)[1]=maxEnd;
            txStart = minStart;
            cdsStart = minStart;
            txEnd = maxEnd;
            cdsEnd = maxEnd;
        }
        if(is_known){
            categorie="full_splice_match";
            subcategorie="gencode";
        }
        
        HashSet<String> hash = new HashSet<String>();
        for (int i=0; i<evidenceList.size(); i++)
            hash.add(evidenceList.get(i).getBarcode());
        
        this.nbUmis = evidenceList.size();
        this.nbCells = hash.size();
    }
    
    private static int[] toIntArray(String str) throws NumberFormatException
    {
        str = StringUtils.stripEnd(str, ",");
        String[] vals = str.split(",");
        int[] numbers = new int[vals.length];

        for (int i = 0; i < vals.length; i++) {
            numbers[i] = Integer.valueOf(vals[i]);
        }
        return numbers;
    }
    
    public int getDistanceTo3p(int genomic_pos)
    {
        int dist = 0;
        
        if(strand.equals(Strand.NEGATIVE)){
            for (int i = 0; i<this.exons.size(); i++) {
                int[] ex = this.exons.get(i);
                
                
                if(ex[0] < genomic_pos){
                    if(ex[1] > genomic_pos)
                        dist += genomic_pos-ex[0];
                    else
                        dist += ex[1]-ex[0];
                }
            }
        }
        else{
            for (int i = 0; i<this.exons.size(); i++) {
                int[] ex = this.exons.get(i);
                
                if(ex[1] > genomic_pos){
                    if(ex[0] < genomic_pos)
                        dist += ex[1]-genomic_pos;
                    else
                        dist += ex[1]-ex[0];
                }
            }
        }
        //System.out.println(transcriptId + "-->" + dist);
    	return dist; 
    }
    
}
