package org.ipmc.sicelore.utils;

/**
 *
 * @author kevin lebrigand
 *
 */
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

    public MoleculeDataset(LongreadParser bam)
    {
        Molecule molecule = null;
        this.mapMolecules = new HashMap<String, Molecule>();
        this.mapGenes = new HashMap<String, List<Molecule>>();
        log = Log.getInstance(MoleculeDataset.class);

        HashMap<String, Longread> mapLongreads = bam.getMapLongreads();
        Set cles = mapLongreads.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String name = (String) it.next();
            Longread lr = (Longread) mapLongreads.get(name);

            String molkey = lr.getBarcode() + ":" + lr.getUmi();
            if((molecule = (Molecule) this.mapMolecules.get(molkey)) != null) {
                molecule.addLongread(lr);
            } 
            else{
                this.mapMolecules.put(molkey, new Molecule(lr.getBarcode(), lr.getUmi()));
                ((Molecule) this.mapMolecules.get(molkey)).addLongread(lr);
            }
        }
        log.info(new Object[]{"\tTotal Molecules\t\t[" + this.mapMolecules.size() + "]"});
        
        int multiIG = 0;
        cles = this.mapMolecules.keySet();
        it = cles.iterator();
        while (it.hasNext()){
            molecule = (Molecule) this.mapMolecules.get((String) it.next());
            if(molecule.getGeneIds().size() > 1)
                multiIG++;
        }        
        log.info(new Object[]{"\tTotal Molecules multiIG\t[" + multiIG + "]"});
    }
    
    public void setIsoforms(UCSCRefFlatParser model, int DELTA, boolean SOFT)
    {
        Molecule molecule = null;
        
        int compteur=0;
        log.info(new Object[]{"\tSetIsoforms\t\tstart..."});
        Set cles = this.mapMolecules.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            molecule = (Molecule) this.mapMolecules.get(molkey);

            List<TranscriptRecord> transcripts = model.select(molecule.getGeneIdsArray());
            // sort all selected isoforms on exons number (max to min)
            Collections.sort(transcripts);
            molecule.setIsoform(transcripts, DELTA, SOFT);
            
            compteur++;
            if(compteur%200000 == 0)
                log.info(new Object[]{"\tSetIsoforms\t\t" + compteur + "/" + this.mapMolecules.size()});
        }
        log.info(new Object[]{"\tSetIsoforms\t\tend..."});
        
        List<Molecule> l;
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
        }
    }

    public List<Molecule> select(String mygene){ return (List<Molecule>) mapGenes.get(mygene); }
    public HashMap<String, Molecule> getMapMolecules() { return this.mapMolecules; }

    public Matrix produceMatrix(UCSCRefFlatParser model, HashSet<String> authorizedCells)
    {
        int nb = 0;
        Matrix matrix = new Matrix(authorizedCells);
        HashMap<String, List<TranscriptRecord>> mapGenesTranscripts = model.getMapGenesTranscripts();

        log.info(new Object[]{"\tDTEMatrix\t\tstart...[" + mapGenesTranscripts.size() + "] genes"});
        Set cles = mapGenesTranscripts.keySet();
        Iterator<String> it = cles.iterator();
        while (it.hasNext()) {
            String mygene = (String) it.next();
            List<Molecule> molecules = this.select(mygene);

            if (molecules != null) {
                for (Molecule molecule : molecules)
                    matrix.addMolecule(molecule);
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
