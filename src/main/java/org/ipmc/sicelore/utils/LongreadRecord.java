package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.util.*;
import org.apache.commons.lang3.StringUtils;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

public class LongreadRecord implements Comparable<LongreadRecord>
{
    private String name;
    private String barcode;
    private String umi;
    //private String chrom;
    //private Strand strand;
    //private int txStart;
    //private int txEnd;
    //private int exonCount;
    //private int exonBases;
    //private int[] exonStarts;
    //private int[] exonEnds;
    private String geneId;
    //private String[] geneIds;
    private List<int[]> exons;
    //private double rpkm;
    private byte[] cdna;
    //private String orientation;
    //private int umiEnd;
    //private int tsoEnd;
    private Float dv;
    //private boolean isSoftOrHardClipped = false;
    //private boolean isSecondaryOrSupplementary = false;
    //private int sizeStartToClip;
    //private int sizeEndToClip;
    //private boolean is_associated = false;
    private boolean isChimeria = false;
    private boolean isReversed = false;
    
    private LongreadRecord() { }

    public int compareTo(LongreadRecord lr){
        Float obj1 = new Float(((LongreadRecord) lr).getDv());
        Float obj2 = new Float(this.dv);
        int retval = obj2.compareTo(obj1);
        return retval;
    }

    public static LongreadRecord fromSAMRecord(SAMRecord r, boolean load_sequence) throws LongreadParseException
    {
        LongreadRecord record = new LongreadRecord();
        
        record.geneId = (String) r.getAttribute("IG");
        record.barcode = (String) r.getAttribute("BC");
        record.umi = (String) r.getAttribute("U8");
       
        if (record.geneId == null || record.barcode == null || record.umi == null)
            return null;
        
        try {
            record.name = r.getReadName();
            // IG flag is Illumina gene Name (from Rainer)
            //if ((String) r.getAttribute("GE") != null) {
            //    record.geneIds = ((String) r.getAttribute("GE")).split(",");
            //}
            
            //record.isIlluminaGeneFound = ((String) r.getAttribute("B0") != null) ? true : false; 
               
            
            //record.chrom = r.getReferenceName();
            //record.txStart = r.getAlignmentStart();
            //record.txEnd = r.getAlignmentEnd();
            //record.isSecondaryOrSupplementary = r.isSecondaryOrSupplementary();
            record.dv = (Float) r.getAttribute("dv");
            
            //Strand strand = Strand.fromString((r.getReadNegativeStrandFlag()) ? "-" : "+");
            //String orientation = (String) r.getAttribute("AR");  // if exists US is "TSO------------------------AAAA-UMI-BC-ADAPTOR"
            int umiEnd = ((Integer) r.getAttribute("UE") != null) ? (Integer) r.getAttribute("UE") : 0; 		// +1 is start of polyA with --------------TTTT read orientation !!!
            int tsoEnd = ((Integer) r.getAttribute("TE") != null) ? (Integer) r.getAttribute("TE") : 0;
            
            //boolean isSoftOrHardClipped = false;
            int sizeStartToClip = 0;
            int sizeEndToClip = 0;
            String exon_starts = "";
            String exon_ends = "";
            List<AlignmentBlock> blocks = r.getAlignmentBlocks();
            String cigar = r.getCigarString();
            String[] cigartype = cigar.split("[0-9]+");
            String[] cigarsize = cigar.split("[A-Z]");
            cigar = cigar.replaceAll("[0-9]+[IDS]","");
            
            // detect large softclipping starting or ending reads
            if ("H".equals(cigartype[1]) || "S".equals(cigartype[1])){
                //if (new Integer(cigarsize[0]).intValue() > 250) {
                //    isSoftOrHardClipped = true;
                    sizeStartToClip = new Integer(cigarsize[0]).intValue();
                //}
            }
            if ("H".equals(cigartype[cigartype.length - 1]) || "S".equals(cigartype[cigartype.length - 1])) {
                //if (new Integer(cigarsize[cigarsize.length - 1]).intValue() > 250) {
                //    isSoftOrHardClipped = true;
                    sizeEndToClip = new Integer(cigarsize[cigarsize.length - 1]).intValue();
                //}
            }
                
            // is possible to set IG if we have a GE but it will induce bad correlation with illumina profiling
            // seems to be difference in term of GTF used for both analysis: cellranger and dropseq.jar add GE tag
            // 50% are Gm.... in this case, very few real interesting genes, few exeample are:
            // Amdhd1, Atxn7l1os1, Bcas3os1, Boll, Ccdc42os, Dnah17, Dnm3os, Greb1, Il4, March10, Mixl1, Prg4, Prr15l, Slc39a5, Tm4sf20
            
            if (record.geneId != null && record.barcode != null && record.umi != null){
                String str = null;
                String readSequence = (String)r.getAttribute("US");
                //record.is_associated = true;
                
                // Strand "+" process
                if(! r.getReadNegativeStrandFlag()){
                    if(sizeEndToClip < 150){
                        if(load_sequence){
                            // need to trim sizeStartToClip bases at start in case of chimeria
                            if(sizeStartToClip > 150){
                                //System.out.println(record.name+"\t"+record.geneId+"\t"+sizeStartToClip+"\t"+tsoEnd+"\n"+readSequence);
                                if(sizeStartToClip > umiEnd)
                                    record.isReversed = true;
                                else
                                    str = readSequence.substring(sizeStartToClip, umiEnd);
                                //System.out.println(str+"\n");
                            }
                            else{
                               if(tsoEnd > umiEnd)
                                   record.isReversed = true;
                               else
                                   str = readSequence.substring((tsoEnd != 0) ? tsoEnd : 0, umiEnd);
                            }
                        }
                    }
                    else{ record.isChimeria = true; }
                }
                // Strand "-" process
                else{
                    if(sizeStartToClip < 150){
                        if(load_sequence){
                            // need to trim sizeStartToClip bases at start in case of chimeria
                            if(sizeEndToClip > 150){
                                //System.out.println(record.name+"\t"+record.geneId+"\t"+sizeEndToClip+"\t"+tsoEnd+"\n"+readSequence);
                                if(sizeEndToClip > umiEnd)
                                    record.isReversed = true;
                                else
                                   str = readSequence.substring(sizeEndToClip, umiEnd);
                                //System.out.println(str+"\n");
                            }
                            else{
                               if(tsoEnd > umiEnd)
                                   record.isReversed = true;
                               else
                                   str = readSequence.substring((tsoEnd != 0) ? tsoEnd : 0, umiEnd);
                            }
                        }
                    }
                    else{ record.isChimeria = true; }
                }
                    /*
                    if(orientation == null) {
                        str = readSequence.substring(umiEnd, (tsoEnd != 0) ? tsoEnd : readSequence.length());
                        str = complementWC(str);
                    }
                    else {
                        str = readSequence.substring((tsoEnd != 0) ? tsoEnd : 0, umiEnd);
                    }
                    if(sizeStartToClip > 0){
			if(sizeStartToClip >= str.length()){ return null; }
			str = str.substring(sizeStartToClip, str.length());
                    }
                    else if(sizeEndToClip > 0){
			if(sizeEndToClip >= str.length()){ return null; }
			str = str.substring(sizeEndToClip, str.length());
                    }
                    */
                if(load_sequence && !record.isChimeria && !record.isReversed)
                    record.cdna = str.getBytes();

                cigar = cigar.replaceAll("[0-9]+[IDS]","");
                cigartype = cigar.split("[0-9]+");
                cigarsize = cigar.split("[A-Z]");

                // init exons
                int block_index = 0;
                int s = blocks.get(0).getReferenceStart();
                int e = blocks.get(0).getLength();
                for (int i = 0; i < cigarsize.length - 1; i++) {
                    AlignmentBlock currBlock = blocks.get(block_index);

                    if ("M".equals(cigartype[i + 1])) {
                        block_index++;
                    }
                    else if ("N".equals(cigartype[i + 1])) {
                        exon_starts += s + ",";
                        exon_ends += e + ",";
                        s = currBlock.getReferenceStart();
                    }
                    e = currBlock.getReferenceStart() + currBlock.getLength();
                }
                exon_starts += s + ",";
                exon_ends += e + ",";

                int[] exonStarts = LongreadRecord.toIntArray(exon_starts);
                int[] exonEnds = LongreadRecord.toIntArray(exon_ends);
                record.exons = new ArrayList<int[]>();
                for (int i = 0; i < exonStarts.length; i++)
                    record.exons.add(new int[]{exonStarts[i], exonEnds[i]});
            }
        } catch (Exception e) { throw new LongreadParseException("Invalid Bam file. " + record.name + ", Can't parse: ", e); }

        return record;
    }
    /*
    public static String complementWC(String dna) {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < dna.length(); i++) {
            char c = dna.charAt(i);
            if (dna.charAt(i) == 'T') {
                builder.append('A');
            }
            if (dna.charAt(i) == 'A') {
                builder.append('T');
            }
            if (dna.charAt(i) == 'C') {
                builder.append('G');
            }
            if (dna.charAt(i) == 'G') {
                builder.append('C');
            }
        }
        return builder.reverse().toString();
    }
    */
    private static int[] toIntArray(String str) throws NumberFormatException {
        str = StringUtils.stripEnd(str, ",");
        String[] vals = str.split(",");
        int[] numbers = new int[vals.length];

        for (int i = 0; i < vals.length; i++) {
            numbers[i] = Integer.valueOf(vals[i]);
        }
        return numbers;
    }

    public String toString()
    {
        //String str = "[" + chrom + ':' + txStart + '-' + txEnd + " --> " + geneId + ", "+exons.size()+" ex ]\n"+ new String(cdna) +"\n";
        String str = "lrr[" + geneId + ", " + exons.size() + " ex]\n";
       // str += "Exon Starts: " + ArrayUtils.toString(exonStarts) + "\n";
       // str += "Exon Ends: " + ArrayUtils.toString(exonEnds) + "\n";
        return str;
    }

    public String getGeneId() {
        return geneId;
    }

    //public String[] getGeneIds() {
    //    return geneIds;
    //}

    public byte[] getCdna() {
        return this.cdna;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public float getDv() {
        return this.dv;
    }
    /*
    public int getSizeStartToClip(){
        return sizeStartToClip;
    }
    
    public int getSizeEndToClip(){
        return sizeEndToClip;
    }
    
    public boolean getIsSoftOrHardClipped() {
        return isSoftOrHardClipped;
    }
    */
    public boolean getIsChimeria() {
        return isChimeria;
    }
    public boolean getIsReversed() {
        return isReversed;
    }
    
    //public boolean getIsSecondaryOrSupplementary() {
    //    return isSecondaryOrSupplementary;
    //}

    public List<int[]> getExons() {
        return this.exons;
    }

    //public boolean getIs_associated() {
    //    return this.is_associated;
    //}

    public String getName() { return name; }
    public void setName(String name) { this.name=name; }

    public String getBarcode() { return barcode; }
    public void setBarcode(String barcode) { this.barcode=barcode; }

    public String getUmi() {return umi; }
    public void setUmi(String umi) { this.umi=umi; }

    //public String getChrom() {
    //    return chrom;
    //}

    //public Strand getStrand() {
    //    return strand;
    //}

    //public int getTxStart() {
    //    return txStart;
    //}

    //public int getTxEnd() {
    //    return txEnd;
    //}
    /*
    public int getExonCount() {
        return exonCount;
    }

    public int getExonBases() {
        return exonBases;
    }

    public int[] getExonStarts() {
        return exonStarts;
    }

    public int[] getExonEnds() {
        return exonEnds;
    }
    */
}
