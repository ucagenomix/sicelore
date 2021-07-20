package org.ipmc.common.gtf;

import java.util.*;

public class GTFgene {

    private String gene_id;
    private String gene_name;
    private String chromosome;
    private int start;
    private int end;
    private String strand;
    private int nb_reads = 0;

    private HashMap<String, GTFtranscript> transcripts;
    private TreeMap<String, GTFexon> exons;
    private TreeMap<String, GTFjunction> junctions;
    private HashMap<String, Integer> isoforms;

    public GTFgene() {
    }

    public GTFgene(String gene_id, String gene_name, String chromosome, int start, int end, String strand) {
        this.gene_id = gene_id;
        this.gene_name = gene_name;
        this.isoforms = new HashMap<String, Integer>();
        this.transcripts = new HashMap<String, GTFtranscript>();
        this.exons = new TreeMap<String, GTFexon>();
        this.junctions = new TreeMap<String, GTFjunction>();
        this.chromosome = chromosome;
        this.start = start;
        this.end = end;
        this.strand = strand;
    }

    public String getGene_id() {
        return this.gene_id;
    }

    public String getGene_name() {
        return this.gene_name;
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

    public int getNb_reads() {
        return this.nb_reads;
    }

    public TreeMap<String, GTFjunction> getJunctions() {
        return this.junctions;
    }

    public TreeMap<String, GTFexon> getExons() {
        return this.exons;
    }

    public GTFtranscript getTranscript(String transcript_id) {
        return this.transcripts.get(transcript_id);
    }

    public void addTranscript(GTFtranscript transcript) {
        this.transcripts.put(transcript.getTranscript_id(), transcript);
    }

    public void compute() {
        GTFjunction j = null;
        Set cles = this.transcripts.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String id = (String) it.next();
            GTFtranscript t = (GTFtranscript) this.transcripts.get(id);
            Vector ve = t.getExons();
            int nb_exons = ve.size();

            if ("+".equals(this.strand)) {
                for (int i = 0; i < nb_exons - 1; i++) {
                    GTFexon e = (GTFexon) ve.get(i);
                    GTFexon enext = (GTFexon) ve.get(i + 1);

                    String erange = new String(e.getStart() + "-" + e.getEnd());
                    String jrange = new String(e.getEnd() + "-" + enext.getStart());

                    if ((j = this.junctions.get(jrange)) == null) {
                        j = new GTFjunction(e.getEnd(), enext.getStart());
                        this.junctions.put(jrange, j);
                    }
                    t.addJunction(j);
                    this.exons.put(erange, e);
                }
            } else {
                for (int i = nb_exons - 2; i >= 0; i--) {
                    GTFexon e = (GTFexon) ve.get(i);
                    GTFexon eprev = (GTFexon) ve.get(i + 1);

                    String erange = new String(e.getStart() + "-" + e.getEnd());
                    String jrange = new String(eprev.getEnd() + "-" + e.getStart());

                    if ((j = this.junctions.get(jrange)) == null) {
                        j = new GTFjunction(eprev.getEnd(), e.getStart());
                        this.junctions.put(jrange, j);
                    }
                    t.addJunction(j);
                    this.exons.put(erange, e);
                }
            }
        }

        int index = 1;
        Set cles2 = this.junctions.keySet();
        Iterator it2 = cles2.iterator();
        while (it2.hasNext()) {
            String id = (String) it2.next();
            j = (GTFjunction) this.junctions.get(id);
            j.setLabel(new String("J" + index));
            j.setCount(0);
            index++;
        }

        index = 1;
        cles2 = this.exons.keySet();
        it2 = cles2.iterator();
        while (it2.hasNext()) {
            String id = (String) it2.next();
            GTFexon e = (GTFexon) this.exons.get(id);
            e.setLabel(new String("E" + index));
            e.setCount(0);
            index++;
        }
    }

    public String isKnownExon(int start, int end) {
        Set cles = this.exons.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String id = (String) it.next();
            GTFexon e = (GTFexon) exons.get(id);

            if (Math.abs(start - e.getStart()) < 20 && Math.abs(end - e.getEnd()) < 20) {
                return e.getExon_number();
            }
        }
        return "#NA";
    }

    public String isKnownJunction(int start, int end) {
        Set cles = this.junctions.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String id = (String) it.next();
            GTFjunction j = (GTFjunction) this.junctions.get(id);

            if (Math.abs(start - j.getStart()) < 5 && Math.abs(end - j.getEnd()) < 5) {
                j.addCount();
                return j.getLabel();
            }
        }
        return "x";
    }

    public void display() {
        Set cles = this.junctions.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String id = (String) it.next();
            GTFjunction j = (GTFjunction) this.junctions.get(id);
            System.out.println(j.getLabel() + "\t" + j.getStart() + "-" + j.getEnd() + " \t" + j.getCount() + " time(s)");
        }

        cles = this.transcripts.keySet();
        it = cles.iterator();
        while (it.hasNext()) {
            String id = (String) it.next();
            GTFtranscript t = (GTFtranscript) this.transcripts.get(id);
            String isoform = t.getIsoform();

            if (this.isoforms.containsKey(isoform)) {
                System.out.println(this.isoforms.get(isoform) + "\t" + t.getTranscript_id() + "\t" + isoform);
            } else {
                System.out.println("0\t" + t.getTranscript_id() + "\t" + isoform);
            }
        }
    }

    public void addIsoform(String list) {
        this.nb_reads++;
        Integer i = null;
        if ((i = this.isoforms.get(list)) != null) {
            this.isoforms.put(list, new Integer(i + 1));
        } else {
            this.isoforms.put(list, new Integer(1));
        }
    }
}
