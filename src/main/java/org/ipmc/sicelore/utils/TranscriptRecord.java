package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.util.*;
import org.apache.commons.lang3.StringUtils;

public class TranscriptRecord implements Comparable<TranscriptRecord> {

    private String transcriptId;
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
    private String geneId;
    private List<int[]> codingExons;
    private List<int[]> exons;
    private double rpkm;

    private TranscriptRecord() {
    }

    public int compareTo(TranscriptRecord tr) {
        // max exons to min
        return ((TranscriptRecord) tr).getExonCount() - this.exonCount;
        // min exons to max
        //return this.exonCount - ((TranscriptRecord)tr).getExonCount();
    }

    public String toString() {
        return "[ transcriptId=" + transcriptId + ", exonCount=" + exonCount + "]";
    }

    public static TranscriptRecord fromRefFlat(String[] fields) throws GTFParseException {
        if (fields.length < 11) {
            throw new RuntimeException(
                    "Invalid RefGene file. records should have at least 11 fields but found only: "
                    + fields.length);
        }

        TranscriptRecord record = new TranscriptRecord();

        record.geneId = fields[0];
        record.transcriptId = fields[1];
        record.chrom = fields[2];
        record.strand = Strand.fromString(fields[3]);

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
                record.exons.add(new int[]{start, end});

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
        } catch (NumberFormatException e) {
            throw new GTFParseException(
                    "Invalid RefGene file. Can't parse integer value: ", e);
        }

        record.rpkm = 0.0;

        return record;
    }

    public void setRPKM(double rpkm) {
        this.rpkm = rpkm;
    }

    public double getRPKM() {
        return this.rpkm;
    }

    public List<int[]> getCodingExons() {
        return this.codingExons;
    }

    public List<int[]> getExons() {
        return this.exons;
    }

    public List<int[]> getExons(boolean cdsExonsOnly) {
        return cdsExonsOnly ? this.codingExons : this.exons;
    }

    private static int[] toIntArray(String str) throws NumberFormatException {
        str = StringUtils.stripEnd(str, ",");
        String[] vals = str.split(",");
        int[] numbers = new int[vals.length];

        for (int i = 0; i < vals.length; i++) {
            numbers[i] = Integer.valueOf(vals[i]);
        }
        return numbers;
    }

    /*
    public String toString() {
        String str = "[\n"; 
        str += transcriptId+": "+chrom + ':' + txStart + '-' + txEnd + " " + strand + "\n";
        str += "Gene: "+geneId+"\n";
        str += "CDS: " +cdsStart+'-'+cdsEnd + "\n";
        str += "Exon Count: " + exonCount + "\n";
        str += "Exon Bases: " + exonBases + "\n";
        str += "CDS Exon Bases: " + cdsExonBases + "\n";
        str += "Exon Starts: "+ ArrayUtils.toString(exonStarts) + "\n";
        str += "Exon Ends: "+ArrayUtils.toString(exonEnds) + "\n";
        //str += "Exon Frames: "+ArrayUtils.toString(exonFrames) + "\n";
        str += "Coding Exons: ";
        for(int[] x: this.getCodingExons()) {
            str += x[0]+"-"+x[1]+",";
        }
       
        str += "\n]\n";
        
        return str;
    }
     */

    public String getTranscriptId() {
        return transcriptId;
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

    public int getCdsStart() {
        return cdsStart;
    }

    public int getCdsEnd() {
        return cdsEnd;
    }

    public int getExonCount() {
        return exonCount;
    }

    public int getExonBases() {
        return exonBases;
    }

    public int getCdsExonBases() {
        return cdsExonBases;
    }

    public int[] getExonStarts() {
        return exonStarts;
    }

    public int[] getExonEnds() {
        return exonEnds;
    }

    public int[] getExonFrames() {
        return exonFrames;
    }

    public String getGeneId() {
        return geneId;
    }

}
