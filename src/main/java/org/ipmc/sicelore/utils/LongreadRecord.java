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
    private Float de;
    private int mapqv;
    //private boolean isSoftOrHardClipped = false;
    private boolean isSecondaryOrSupplementary = false;
    //private int sizeStartToClip;
    //private int sizeEndToClip;
    //private boolean is_associated = false;
    private boolean isChimeria = false;
    private boolean isReversed = false;
    
    protected static String CELLTAG = "BC"; // BC
    protected static String UMITAG = "U8"; // U8
    protected static String GENETAG = "IG"; //IG
    protected static String TSOENDTAG = "TE"; //TE
    protected static String UMIENDTAG = "UE"; //UE
    protected static String USTAG = "US"; // US
    protected static int MAXCLIP = 150; // 150

    public LongreadRecord() { }

    public void setStaticParams(String celltag, String umitag, String genetag, String tsoendtag, String umiendtag, String ustag, int maxclip){
	this.CELLTAG = celltag;
        this.UMITAG = umitag;
        this.GENETAG = genetag;
        this.TSOENDTAG = tsoendtag;
        this.UMIENDTAG = umiendtag;
        this.USTAG = ustag;
        this.MAXCLIP = maxclip;
    }

    public int compareTo(LongreadRecord lr){
        Float obj1 = new Float(((LongreadRecord) lr).getDe());
        Float obj2 = new Float(this.de);
        int retval = obj2.compareTo(obj1);
        return retval;
    }

    public static LongreadRecord fromSAMRecord(SAMRecord r, boolean load_sequence) throws LongreadParseException
    {
        LongreadRecord record = new LongreadRecord();
        
        record.geneId = (String) r.getAttribute(GENETAG);
        record.barcode = (String) r.getAttribute(CELLTAG);
        record.umi = (String) r.getAttribute(UMITAG);
        record.mapqv = r.getMappingQuality();
        
        if (record.barcode == null || record.umi == null || r.getReadUnmappedFlag())
             return null;
        
        try {
            record.name = r.getReadName();
            //record.chrom = r.getReferenceName();
            //record.txStart = r.getAlignmentStart();
            //record.txEnd = r.getAlignmentEnd();
            record.isSecondaryOrSupplementary = r.isSecondaryOrSupplementary();
            record.de = ((Float) r.getAttribute("de") != null) ? (Float) r.getAttribute("de") : (Float) r.getAttribute("df"); // minimap 2.17 (de) versus 2.10 (df)
            
            //Strand strand = Strand.fromString((r.getReadNegativeStrandFlag()) ? "-" : "+");
            //String orientation = (String) r.getAttribute("AR");  // if exists US is "TSO------------------------AAAA-UMI-BC-ADAPTOR"
            int umiEnd = ((Integer) r.getAttribute(UMIENDTAG) != null) ? (Integer) r.getAttribute(UMIENDTAG) : 0; 		// +1 is start of polyA with --------------TTTT read orientation !!!
            int tsoEnd = ((Integer) r.getAttribute(TSOENDTAG) != null) ? (Integer) r.getAttribute(TSOENDTAG) : 0;
            
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
            
            String str = null;
            String readSequence = (String)r.getAttribute(USTAG);

            // detect softclipping starting or ending reads
            if ("H".equals(cigartype[1]) || "S".equals(cigartype[1])){ sizeStartToClip = new Integer(cigarsize[0]).intValue(); }
            if ("H".equals(cigartype[cigartype.length - 1]) || "S".equals(cigartype[cigartype.length - 1])) { sizeEndToClip = new Integer(cigarsize[cigarsize.length - 1]).intValue(); }

            // Strand "+" process
            if(! r.getReadNegativeStrandFlag()){
                if(sizeEndToClip < MAXCLIP){
                    if(load_sequence){
                        // need to trim sizeStartToClip bases at start in case of chimeria
                        if(sizeStartToClip > MAXCLIP){
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
                if(sizeStartToClip < MAXCLIP){
                    if(load_sequence){
                        // need to trim sizeStartToClip bases at start in case of chimeria
                        if(sizeEndToClip > MAXCLIP){
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

            if(sizeStartToClip > MAXCLIP || sizeEndToClip > MAXCLIP)
                record.isChimeria = true;

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

    public int getMapqv() {
        return this.mapqv;
    }
    public float getDe() {
        return this.de;
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
    
    public boolean getIsSecondaryOrSupplementary() {
        return isSecondaryOrSupplementary;
    }

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
