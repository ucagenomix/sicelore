package org.ipmc.sicelore.utils;

import java.util.*;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.CloseableIterator;

public class LongreadRecord implements Comparable<LongreadRecord>
{
    private String name;
    private String barcode;
    private String umi;
    private String chrom;
    private Strand strand;
    private int txStart;
    private int txEnd;
    private int exonCount;
    private int exonBases;
    private int[] exonStarts;
    private int[] exonEnds;
    private String geneId;
    private String[] geneIds;
    private List<int []> exons;
    private double rpkm;
    private String cdna;
    private String orientation;
    private int umiEnd;
    private int tsoEnd;
    private Float dv;
    private boolean isSoftOrHardClipped = false;
    private boolean isSecondaryOrSupplementary = false;
    private int sizeStartToClip;
    private int sizeEndToClip;
    private boolean is_associated = false;

    private LongreadRecord() {}

    public int compareTo(LongreadRecord lr)
    {
    	Float obj1 = new Float(((LongreadRecord)lr).getDv());
        Float obj2 = new Float(this.dv);
        int retval =  obj2.compareTo(obj1);
    	return retval;
    }

    public static LongreadRecord fromSAMRecord(SAMRecord r) throws LongreadParseException
    {
        LongreadRecord record = new LongreadRecord();
        record.name = r.getReadName();
        
        try {     
        	// IG flag is Illumina gene Name (from Rainer)
        	if((String)r.getAttribute("GE") != null)
        		record.geneIds = ((String)r.getAttribute("GE")).split(",");
        	
        	record.geneId = (String)r.getAttribute("IG");
        	record.barcode = (String)r.getAttribute("BC");
        	record.umi = (String)r.getAttribute("U8");
        	record.chrom = r.getReferenceName();
            record.txStart = r.getAlignmentStart();
            record.txEnd = r.getAlignmentEnd();
            record.strand = Strand.fromString((r.getReadNegativeStrandFlag())?"-":"+");
            record.isSecondaryOrSupplementary = r.isSecondaryOrSupplementary();
            record.dv = (Float)r.getAttribute("dv");
			record.orientation = (String)r.getAttribute("AR");  // if exists US is "TSO------------------------AAAA-UMI-BC-ADAPTOR"
			record.umiEnd = ((Integer)r.getAttribute("UE") != null)?(Integer)r.getAttribute("UE"):0; 		// +1 is start of polyA with --------------TTTT read orientation !!!
			record.tsoEnd = ((Integer)r.getAttribute("TE") != null)?(Integer)r.getAttribute("TE"):0;

            String exon_starts = "";
            String exon_ends = "";
			List<AlignmentBlock> blocks = r.getAlignmentBlocks();
		 	String cigar = r.getCigarString();
			String[] cigartype = cigar.split("[0-9]+");
			String[] cigarsize = cigar.split("[A-Z]");
		 	
			// detect large softclipping starting or ending reads
			if("H".equals(cigartype[1]) || "S".equals(cigartype[1])){
				if(new Integer(cigarsize[0]).intValue() > 450){
					record.isSoftOrHardClipped = true;
					record.sizeStartToClip = new Integer(cigarsize[0]).intValue();
					//System.out.println(record.name+"\t"+ r.isSecondaryOrSupplementary() +"\tstart\t"+cigarsize[0]+"\t"+cigartype[1]);
				}
			}
			if("H".equals(cigartype[cigartype.length-1]) || "S".equals(cigartype[cigartype.length-1])){
				if(new Integer(cigarsize[cigarsize.length-1]).intValue() > 450){
					record.isSoftOrHardClipped = true;
					record.sizeEndToClip = new Integer(cigarsize[cigarsize.length-1]).intValue();
					//System.out.println(record.name+"\t"+r.isSecondaryOrSupplementary()+"\tend\t"+cigarsize[cigarsize.length-1]+"\t"+cigartype[cigartype.length-1]);
				}
			}
			
        	// if we don't have a IG / BC / U8 we do not init cDNA
        	if(record.geneId != null && record.barcode != null && record.umi != null)
        	{
        		record.is_associated = true;

	        	// get the read sequence but not as record attribute,will set the cDNA sequence instead.
				//need to be record.is_associated = true !!!
	            String readSequence = (String)r.getAttribute("US");
	            if(readSequence != null && record.is_associated){
					if(record.orientation == null){
						record.cdna = readSequence.substring(record.umiEnd,(record.tsoEnd != 0)?record.tsoEnd:readSequence.length());
						record.cdna = complementWC(record.cdna);
					}
					else{ record.cdna = readSequence.substring((record.tsoEnd != 0)?record.tsoEnd:0,record.umiEnd); }	
					
					// can not be a good record molecule if clipped on both sides !
					//if(record.sizeStartToClip > 0 && record.sizeEndToClip > 0){ return null; }
					
					/*
					// Clip sequence if hard or soft clipped alignment on one side only otherwise return null, 
					else if(record.sizeStartToClip > 0){
						// case of ADAP-BC-UMI-TTTT map the genome....
						if(record.sizeStartToClip >= record.cdna.length()){ return null; }
						record.cdna = record.cdna.substring(record.sizeStartToClip, record.cdna.length());
					}
					else if(record.sizeEndToClip > 0){
						// case of ADAP-BC-UMI-TTTT map the genome....
						if(record.sizeStartToClip >= record.cdna.length()){ return null; }
						record.cdna = record.cdna.substring(record.sizeEndToClip, record.cdna.length());
					}
					*/
	            }
        	}
	            
            
			// init exons
			int block_index=0;
			int s = blocks.get(0).getReferenceStart();
			int e = blocks.get(0).getLength();
			for(int i=0; i<cigartype.length - 2; i++){
				AlignmentBlock currBlock = blocks.get(block_index);
				
				if("M".equals(cigartype[i+1]))
					block_index++;
				else if("N".equals(cigartype[i+1])){
					exon_starts += s + ",";
					exon_ends += e + ",";
					s = currBlock.getReferenceStart();
				}
				e = currBlock.getReferenceStart() + currBlock.getLength();
     		}
			exon_starts += s + ",";
			exon_ends += e + ",";
			
            record.exonStarts = LongreadRecord.toIntArray(exon_starts);
            record.exonEnds = LongreadRecord.toIntArray(exon_ends);
            record.exonBases = 0;
            record.exons = new ArrayList<int []>();
            
            for(int i=0; i<record.exonStarts.length; i++) {
                int start = record.exonStarts[i];
                int end = record.exonEnds[i];
                
                record.exonBases += end-start;
                record.exons.add(new int[]{start,end});
            }
            record.exonCount = record.exons.size();
            
        } catch (Exception e) { throw new LongreadParseException("Invalid Bam file. " + record.name + ", Can't parse: ", e); }

        return record;
    }

    public static String complementWC(String dna)
    {
        StringBuilder builder = new StringBuilder();
        for(int i=0;i<dna.length();i++){
            char c = dna.charAt(i);
            if(dna.charAt(i) == 'T'){
                builder.append('A');
            }
            if(dna.charAt(i) == 'A'){
                builder.append('T');
            }
            if(dna.charAt(i) == 'C'){
                builder.append('G');
            }
            if(dna.charAt(i) == 'G'){
                builder.append('C');
            }   
        }
        return builder.reverse().toString();
    }
        
    private static int[] toIntArray(String str) throws NumberFormatException
    {
        str = StringUtils.stripEnd(str, ",");
        String[] vals = str.split(",");
        int[] numbers = new int[vals.length];
        
        for(int i = 0; i < vals.length; i++) {
            numbers[i] = Integer.valueOf(vals[i]);
        }
        return numbers;
    }
    
    public String toString()
    {
        String str = "[\n"; 
        str += name+": "+chrom + ':' + txStart + '-' + txEnd + " " + strand + "\n";
        str += "Gene: "+geneId+"\n";
        str += "Exon Count: " + exonCount + "\n";
        str += "Exon Bases: " + exonBases + "\n";
        str += "Exon Starts: "+ ArrayUtils.toString(exonStarts) + "\n";
        str += "Exon Ends: "+ArrayUtils.toString(exonEnds) + "\n";
        str += "\n]\n";
        
        return str;
    }

    public String getGeneId() {  return geneId; }
    public String[] getGeneIds() {  return geneIds; }
    public String getCdna(){ return this.cdna; }
    public void setGeneId(String geneId) { this.geneId=geneId; }
    public float getDv() { return this.dv; }
    public boolean getIsSoftOrHardClipped(){ return isSoftOrHardClipped; }
    public boolean getIsSecondaryOrSupplementary(){ return isSecondaryOrSupplementary; }
    public List<int []> getExons() { return this.exons; }
    public boolean getIs_associated(){ return this.is_associated; }

    public String getName() {
        return name;
    }
    
    public String getBarcode() {
        return barcode;
    }
    
    public String getUmi() {
        return umi;
    }
    
    public String getChrom() {
        return chrom;
    }

    public Strand getStrand() {
        return strand;
    }

    public int getTxStart() {
        return txStart;
    }

    public int getTxEnd() {
        return txEnd;
    }

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
}