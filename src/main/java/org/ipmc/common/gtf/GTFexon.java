package org.ipmc.common.gtf;

public class GTFexon {

    private String exon_id;
    private String exon_number;
    private String chromosome;
    private int start;
    private int end;
    private String strand;

    public int count;
    public String label;

    public GTFexon() {
    }

    public GTFexon(String exon_id, String exon_number, String chromosome, int start, int end, String strand) {
        this.exon_id = exon_id;
        this.exon_number = exon_number;
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }

    public String getExon_id() {
        return this.exon_id;
    }

    public String getExon_number() {
        return this.exon_number;
    }

    public String getChromosome() {
        return this.chromosome;
    }

    public int getStart() {
        return this.start;
    }

    public int getEnd() {
        return this.end;
    }

    public String getStrand() {
        return this.strand;
    }

    public int getCount() {
        return this.count;
    }

    public void setCount(int count) {
        this.count = count;
    }

    public String getLabel() {
        return this.label;
    }

    public void setLabel(String label) {
        this.label = label;
    }

    public void addCount() {
        this.count++;
    }
}
