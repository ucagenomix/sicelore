package org.ipmc.common.gtf;

import java.util.*;
import java.io.*;
import java.util.regex.Pattern;

public class GTFparser {

    private File gtf_file;
    private int nb_genes = 0;
    private int nb_transcripts = 0;
    private int nb_exons = 0;
    private int nb_uniq_exons = 0;
    private int nb_uniq_junctions = 0;

    private HashMap<String, GTFgene> genes;
    private HashMap<String, String> name2id;

    public GTFparser() {
    }

    public GTFparser(File gtf_file) {
        this.gtf_file = gtf_file;
        this.genes = new HashMap<String, GTFgene>();
        this.name2id = new HashMap<String, String>();

        load();
        init();
    }

    public int getNb_genes() {
        return this.nb_genes;
    }

    public int getNb_transcripts() {
        return this.nb_transcripts;
    }

    public int getNb_exons() {
        return this.nb_exons;
    }

    public int getNb_uniq_exons() {
        return this.nb_uniq_exons;
    }

    public int getNb_uniq_junctions() {
        return this.nb_uniq_junctions;
    }

    public GTFgene getGene(String gene_name) {
        GTFgene g = null;
        if (this.name2id.containsKey(gene_name)) {
            g = (GTFgene) this.genes.get((String) this.name2id.get(gene_name));
        }

        return g;
    }

    public void load() {
        BufferedReader gtffile = null;

        try {
            gtffile = new BufferedReader(new FileReader(gtf_file));
            String line;
            String[] tokens;

            while ((line = gtffile.readLine()) != null) {
                tokens = line.split("[\t]");

                if (Pattern.matches("^#.*", line) == false && Pattern.matches("^[\\dmxyeg].*", tokens[0])) {
                    String[] tmp = tokens[8].split(";[ ]*");
                    HashMap<String, String> genetranscriptTokens = new HashMap();
                    for (int i = 0; i < tmp.length; i++) {
                        String[] tmp2 = tmp[i].split(" ");
                        genetranscriptTokens.put(tmp2[0], tmp2[1].replaceAll("\"", ""));
                    }

                    String gene_id = (String) genetranscriptTokens.get("gene_id");
                    String transcript_id = (String) genetranscriptTokens.get("transcript_id");
                    String gene_name = gene_id;
                    String transcript_name = transcript_id;
                    String exon_id = "";
                    String exon_number = "";
                    if (genetranscriptTokens.containsKey("gene_name")) {
                        gene_name = (String) genetranscriptTokens.get("gene_name");
                    }
                    if (genetranscriptTokens.containsKey("transcript_name")) {
                        transcript_name = (String) genetranscriptTokens.get("transcript_name");
                    }
                    if (genetranscriptTokens.containsKey("exon_id")) {
                        exon_id = (String) genetranscriptTokens.get("exon_id");
                    }
                    if (genetranscriptTokens.containsKey("exon_number")) {
                        exon_number = (String) genetranscriptTokens.get("exon_number");
                    }

                    if (tokens[2].equals("gene")) {
                        GTFgene g = null;
                        if ((g = (GTFgene) this.genes.get(gene_id)) == null) {
                            GTFgene gene = new GTFgene(gene_id, gene_name, tokens[0], new Integer(tokens[3]).intValue(), new Integer(tokens[4]).intValue(), tokens[6]);
                            this.genes.put(gene_id, gene);
                            nb_genes++;
                            name2id.put(gene_name, gene_id);
                        }
                    } else if (tokens[2].equals("transcript")) {
                        GTFtranscript t = null;
                        GTFgene g = (GTFgene) this.genes.get(gene_id);
                        if ((t = (GTFtranscript) g.getTranscript(transcript_id)) == null) {
                            GTFtranscript transcript = new GTFtranscript(transcript_id, transcript_name, tokens[0], new Integer(tokens[3]).intValue(), new Integer(tokens[4]).intValue(), tokens[6]);
                            g.addTranscript(transcript);
                            nb_transcripts++;
                            //System.out.println("Push transcript\t"+transcript_id+" for gene "+gene_id);
                        }
                    } else if (tokens[2].equals("exon")) {
                        GTFexon e = null;
                        GTFgene g = (GTFgene) this.genes.get(gene_id);
                        GTFtranscript t = (GTFtranscript) g.getTranscript(transcript_id);
                        GTFexon exon = new GTFexon(exon_id, exon_number, tokens[0], new Integer(tokens[3]).intValue(), new Integer(tokens[4]).intValue(), tokens[6]);
                        t.addExon(exon);
                        nb_exons++;
                    } else if (tokens[2].equals("cds")) {

                    } else if (tokens[2].equals("start_codon")) {

                    } else if (tokens[2].equals("stop_codon")) {

                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void init() {
        int nb_exons = 0;
        int nb_junctions = 0;
        Set cles = this.genes.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String id = (String) it.next();
            GTFgene g = (GTFgene) this.genes.get(id);

            g.compute();

            nb_uniq_exons += g.getExons().size();
            nb_uniq_junctions += g.getJunctions().size();
        }
    }

    public void displayStats() {
        System.out.println("\n# " + nb_genes + " genes #");
        System.out.println("# " + nb_transcripts + " transcripts #");
        System.out.println("# " + nb_exons + " exons #");
        System.out.println("# " + nb_uniq_exons + " exons uniq #");
        System.out.println("# " + nb_uniq_junctions + " junctions uniq #\n");
    }
}
