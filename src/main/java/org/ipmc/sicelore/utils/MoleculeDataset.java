package org.ipmc.sicelore.utils;

import java.util.*;
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.util.concurrent.*;
import org.biojava.nbio.core.util.ConcurrencyTools;
import com.google.common.util.concurrent.*;

public class MoleculeDataset {

    private final Log log;

    HashMap<String, Molecule> mapMolecules;
    HashMap<String, List<Molecule>> mapGenes;

    private int chromstrange = 0;
    private int noexons = 0;
    private int nogeneid = 0;
    private int nobarcode = 0;
    private int noumi = 0;
    private int softclipped = 0;
    private int genegm = 0;

    private ListeningExecutorService oneNanoporeReadexecutor;
    private Deque<Future<String>> future_list;
    private DataOutputStream os;
    private Iterator itglobal;

    public MoleculeDataset(LongreadParser bam, UCSCRefFlatParser model, int DELTA, boolean SOFT) {
        Molecule molecule = null;
        this.mapMolecules = new HashMap<String, Molecule>();
        this.mapGenes = new HashMap<String, List<Molecule>>();

        log = Log.getInstance(MoleculeDataset.class);

        // 1. Load molecule from bam file
        HashMap<String, Longread> mapLongreads = bam.getMapLongreads();
        Set cles = mapLongreads.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String name = (String) it.next();
            Longread lr = (Longread) mapLongreads.get(name);

            String molkey = lr.getBarcode() + ":" + lr.getUmi();
            if ((molecule = (Molecule) this.mapMolecules.get(molkey)) != null) {
                molecule.addLongread(lr);
            } else {
                this.mapMolecules.put(molkey, new Molecule(lr.getBarcode(), lr.getUmi()));
                ((Molecule) this.mapMolecules.get(molkey)).addLongread(lr);
            }
        }
        log.info(new Object[]{"\tTotal Molecules\t\t[" + this.mapMolecules.size() + "]"});

        // 2. Remove molecules without any associated GI/BC/U8 read
        cles = this.mapMolecules.keySet();
        it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            molecule = (Molecule) this.mapMolecules.get(molkey);

            boolean remove_the_molecule = true;
            List<Longread> longreads = molecule.getLongreads();
            for (Longread lr : longreads) {
                if (lr.getIs_associated()) {
                    remove_the_molecule = false;
                }
            }
            if (remove_the_molecule) {
                it.remove();
            }
        }
        log.info(new Object[]{"\tGI/BC/U8 Molecules\t[" + this.mapMolecules.size() + "]"});

        // 3. clean chimeria reads from molecules based on US length deviation
        // eventually remove the molecule if no more supporting reads.
        // TODO !
        // 4. set final geneId and transcriptId without picking one random in the list
        // should be possible if we have here the model file and select all transcripts from the x genes annotatng the molecule
        // TODO ! or NOT, if need the model, we will need it everytime....
        HashMap<String, List<TranscriptRecord>> mapGenesTranscripts = model.getMapGenesTranscripts();

        log.info(new Object[]{"\tSetIsoforms\t\tstart..."});
        // Running thought all molecules to set for isoform
        cles = this.mapMolecules.keySet();
        it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            molecule = (Molecule) this.mapMolecules.get(molkey);

            List<TranscriptRecord> transcripts = model.select(molecule.getGeneIdsArray());
            // sort all selected isoforms on exons number (max to min)
            Collections.sort(transcripts);
            molecule.setIsoform(transcripts, DELTA, SOFT);
        }
        log.info(new Object[]{"\tSetIsoforms\t\tend..."});

        // 5. magGenes initialization for rapid Molecule retrieval
        List<Molecule> l;
        int[] multiviews = new int[12];
        cles = this.mapMolecules.keySet();
        it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            molecule = (Molecule) this.mapMolecules.get(molkey);
            if ((l = mapGenes.get(molecule.getGeneId())) != null) {
                l.add(molecule);
            } else {
                l = new ArrayList<Molecule>();
                l.add(molecule);
                mapGenes.put(molecule.getGeneId(), l);
            }
            // stats
            int nn = molecule.getLongreads().size();
            if (nn > 10) {
                nn = 10;
            }
            multiviews[nn]++;
        }
        log.info(new Object[]{"\tx-reads\t\t[1:" + multiviews[1] + ",2:" + multiviews[2] + ",3:" + multiviews[3] + ",4:" + multiviews[4] + ",5:" + multiviews[5] + ",6:" + multiviews[6] + ",7:" + multiviews[6] + ",8:" + multiviews[8] + ",9:" + multiviews[9] + ",10:" + multiviews[10] + "]"});
    }

    public LongreadRecord parseSAMRecord(SAMRecord r) throws LongreadParseException {
        LongreadRecord record = LongreadRecord.fromSAMRecord(r);
        //if(record == null) { softclipped++; return null; }
        //if(record.getChrom().contains("_")) { chromstrange++; return null; }
        //if(record.getExonBases() == 0) { noexons++; return null; }
        //if(record.getGeneId() == null) { nogeneid++; return null; }
        //if(record.getBarcode() == null) { nobarcode++; return null; }
        //if(record.getUmi() == null) { noumi++; return null; }
        //if(record.getIsSoftClipped()) { softclipped++; return null; }
        //if(record.getGeneId().substring(0,2).equals("Gm")) { genegm++; return null; }
        return record;
    }

    public List<Molecule> select(String mygene) {
        return (List<Molecule>) mapGenes.get(mygene);
        /*
    	List<Molecule> filteredList = new ArrayList<Molecule>();
        Set cles = this.mapMolecules.keySet();
        Iterator it = cles.iterator();
        while(it.hasNext()){
        	String cle = (String)it.next();
        	Molecule molecule =(Molecule)this.mapMolecules.get(cle);
        	
        	if(mygene.equals(molecule.getGeneId()))
        		filteredList.add(molecule);
        }
        return filteredList;
         */
    }

    public void displayMetrics(File OUTPUT) {
        Integer i = null;
        DataOutputStream os = null;

        int[] count = new int[11];
        double[] somme = new double[11];

        Set cles = this.mapMolecules.keySet();
        Iterator<String> it = cles.iterator();
        while (it.hasNext()) {
            String key = (String) it.next();
            Molecule molecule = (Molecule) mapMolecules.get(key);
            List<Longread> longreads = molecule.getLongreads();

            float min_dv = 1;
            int xtimes = longreads.size();
            for (Longread lr : longreads) {
                List<LongreadRecord> longreadrecords = lr.getLongreadrecords();
                for (LongreadRecord lrr : longreadrecords) {
                    float dv = lrr.getDv();
                    if (dv < min_dv) {
                        min_dv = dv;
                    }
                }
            }
            double pctId = 1.0 - min_dv;
            if (xtimes >= 10) {
                xtimes = 10;
            }
            count[xtimes]++;
            somme[xtimes] += pctId;
        }

        try {
            log.info(new Object[]{"\tx_times\tavg_dv\tnb_mol\tsum"});
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            for (int ii = 1; ii < 11; ii++) {
                double avg = 0.0;
                if (count[ii] > 0) {
                    avg = somme[ii] / count[ii];
                }
                String dvres = String.format("%.3f", avg);
                String sumres = String.format("%.2f", somme[ii]);
                os.writeBytes("x" + ii + "\t" + dvres + "\t" + count[ii] + "\t" + sumres + "\n");
                log.info(new Object[]{"\tx" + ii + "\t" + dvres + "\t" + count[ii] + "\t" + sumres});
            }
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                os.close();
            } catch (Exception e) {
                System.err.println("can not close stream");
            }
        }
    }

    public HashMap<String, Molecule> getMapMolecules() {
        return this.mapMolecules;
    }

    public Matrix DTEMatrix(UCSCRefFlatParser model) {
        int nb = 0;
        Matrix matrix = new Matrix();
        HashMap<String, List<TranscriptRecord>> mapGenesTranscripts = model.getMapGenesTranscripts();

        log.info(new Object[]{"\tDTEMatrix\t\tstart...[" + mapGenesTranscripts.size() + "] genes"});
        Set cles = mapGenesTranscripts.keySet();
        Iterator<String> it = cles.iterator();
        while (it.hasNext()) {
            String mygene = (String) it.next();
            List<Molecule> molecules = this.select(mygene);

            // if we have molecule for this gene
            if (molecules != null) {
                for (Molecule molecule : molecules) {
                    matrix.addMolecule(molecule);
                }
            }

            nb++;
            if (nb % 10000 == 0) {
                log.info(new Object[]{"\tDTEMatrix\t\t[" + nb + "/" + mapGenesTranscripts.size() + "] genes processed"});
            }
        }
        log.info(new Object[]{"\tDTEMatrix\t\t[" + nb + "/" + mapGenesTranscripts.size() + "] genes processed"});

        return matrix;
    }

    public void callConsensus(File OUTPUT, int nThreads) {
        log.info(new Object[]{"\tCalling consensus\tstart with [" + nThreads + "] threads"});

        future_list = new ConcurrentLinkedDeque<Future<String>>();
        oneNanoporeReadexecutor = MoreExecutors.listeningDecorator(Executors.newFixedThreadPool(nThreads));

        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));

            Set cles = this.mapMolecules.keySet();
            itglobal = cles.iterator();
            int i = 0;
            while (itglobal.hasNext() && i++ < (5 * nThreads)) {
                String key = (String) itglobal.next();
                Molecule molecule = (Molecule) this.mapMolecules.get(key);
                ListenableFuture<String> submit = oneNanoporeReadexecutor.submit(molecule);
                future_list.add(submit); //adds to end of queue
                //System.out.println("add:"+molecule.getUmi()+"\t"+ future_list.size());
            }
            ConcurrencyTools.setThreadPoolDefault();

            try {
                //wait until all jobs finished
                while (future_list.isEmpty() == false) {

                    //the following line blocks until job is done. Can also recover results here get() returns OneThreadResult
                    String rslt = future_list.remove().get();
                    write(rslt);

                    if (itglobal.hasNext()) {
                        String key = (String) itglobal.next();
                        Molecule molecule = (Molecule) this.mapMolecules.get(key);
                        ListenableFuture<String> submit = oneNanoporeReadexecutor.submit(molecule);
                        future_list.add(submit);//adds to end of queue
                    }
                }
            } catch (Exception e) {
                e.printStackTrace(System.out);
                System.exit(1);
            }

            oneNanoporeReadexecutor.shutdown();
            ConcurrencyTools.shutdown();

            try {
                oneNanoporeReadexecutor.awaitTermination(1, TimeUnit.MINUTES);
            } catch (InterruptedException ex) {
                log.info(new Object[]{"Error:\t[" + ex + "]"});
            }

            os.close();
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                os.close();
            } catch (Exception e) {
                System.err.println("can not close stream");
            }
        }

    }

    private synchronized void write(String rslt) {
        try {
            os.writeBytes(rslt);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
