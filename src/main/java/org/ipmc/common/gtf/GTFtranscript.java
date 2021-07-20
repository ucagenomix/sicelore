package org.ipmc.common.gtf;

import java.util.*;

public class GTFtranscript {

    private String transcript_id;
    private String transcript_name;
    private String chromosome;
    private int start;
    private int end;
    private String strand;

    private Vector exons;
    private Vector junctions;

    public GTFtranscript() {
    }

    public GTFtranscript(String transcript_id, String transcript_name, String chromosome, int start, int end, String strand) {
        this.junctions = new Vector();
        this.exons = new Vector();

        this.transcript_id = transcript_id;
        this.transcript_name = transcript_name;
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }

    public String getTranscript_id() {
        return this.transcript_id;
    }

    public String getTranscript_name() {
        return this.transcript_name;
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

    public Vector getJunctions() {
        return this.junctions;
    }

    public Vector getExons() {
        return this.exons;
    }

    public void addExon(GTFexon exon) {
        this.exons.add(exon);
    }

    public void addJunction(GTFjunction j) {
        this.junctions.add(j);
    }

    public String getIsoform() {
        String list = "";

        for (int i = 0; i < this.junctions.size(); i++) {
            GTFjunction j = (GTFjunction) this.junctions.get(i);
            list += j.getLabel() + "-";
        }
        list = list.substring(0, list.length() - 1);

        return list;
    }

}
