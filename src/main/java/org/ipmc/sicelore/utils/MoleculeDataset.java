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
import org.apache.commons.lang3.StringUtils;
import org.biojava.nbio.core.sequence.DNASequence;

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
    
    public List<Molecule> select(String mygene){ return (List<Molecule>) mapGenes.get(mygene); }
    public HashMap<String, Molecule> getMapMolecules() { return this.mapMolecules; }
    
    public Molecule getMolecule(String isokey)
    {
        return this.mapMolecules.get(isokey);
    }
    
    public void setIsoforms(UCSCRefFlatParser model, int DELTA, boolean SOFT, String METHOD)
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
            //Collections.sort(transcripts);
            
            this.setIsoform(molecule, transcripts, DELTA, SOFT, METHOD);
            
            // record womewhere what happened during this setIsoform
            
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

    public void setIsoform(Molecule molecule, List<TranscriptRecord> transcripts, int DELTA, boolean SOFT, String METHOD)
    {
        List list1=null;
        
        HashMap<String, Integer> candidates = new HashMap<String, Integer>();
        List<Longread> longreads = molecule.getLongreads();
        
        for(Longread lr : longreads){
            List<LongreadRecord> records = lr.getLongreadrecords();
            //LongreadRecord lrr = lr.getBestRecord();

            for(LongreadRecord lrr : records){
                List list = junctionsFromExons(lrr.getExons());

                for(TranscriptRecord transcriptrecord : transcripts){
                    list1 = junctionsFromExons(transcriptrecord.getExons());

                    if(map(list, list1, DELTA, SOFT)) {
                        String key = transcriptrecord.getTranscriptId() + "|" + transcriptrecord.getGeneId();

                        // case SOFT pocess put nb exons as the value to get the more complex transcript
                        if(SOFT)
                            candidates.put(key, transcriptrecord.getExons().size());
                        else{
                            // candidates will be ranked on the number of  maximum exons,
                            // nice when 5' degradation or quality loss of long read
                            if("EXONHIGH".equals(METHOD)){
                                candidates.put(key, transcriptrecord.getExons().size());
                            }
                            // candidates will be ranked on the number of read to isoform match
                            // consistent option, but week when last 5p exons is lost because of long read QV loss
                            else if("NBHIGH".equals(METHOD)){
                                if(candidates.containsKey(key))
                                    candidates.put(key, candidates.get(key) + 1);
                                else
                                    candidates.put(key, 1);
                            }
                            //
                            // smart method on collapsing model to determine rules 
                            // of minimal assignation yet to be implemented here !!!
                            //
                        }
                    }
                }
            }
        }
        
        if(candidates.size() > 0){
            HashSet<String> bestCandidates = new HashSet<String>();
            int maxValueInMap=(Collections.max(candidates.values()));  // This will return max value in the Hashmap
            for(Map.Entry<String, Integer> entry : candidates.entrySet()) {  // Iterate through hashmap
                if (entry.getValue() == maxValueInMap) { bestCandidates.add(entry.getKey()); }
            }
            /*
            if(bestCandidates.size() > 1){
                System.out.println(this.barcode + ":" + this.umi + " " + bestCandidates.size());
                Iterator<String> iterator = bestCandidates.iterator();
                while (iterator.hasNext()) {
                    System.out.println((String)iterator.next());
                }
            }
            */
            int index = new Random().nextInt(bestCandidates.size());
            Iterator<String> iter = bestCandidates.iterator();
            for (int i = 0; i < index; i++) { iter.next(); }
            String g = (String) iter.next();
            molecule.setTranscriptId(g.split("\\|")[0]);
            molecule.setGeneId(g.split("\\|")[1]);
        }
        else{
            int index = new Random().nextInt(molecule.getGeneIds().size());
            Iterator<String> iter = molecule.getGeneIds().iterator();
            for (int i = 0; i < index; i++) { iter.next(); }
            molecule.setTranscriptId("undef");
            molecule.setGeneId((String) iter.next());
        }
        
        // only 1 mono-exonic transcript in the model --> we set the isoform
        if(transcripts.size() == 1 && list1.size() == 0){
           molecule.setTranscriptId(transcripts.get(0).getTranscriptId());
           molecule.setGeneId(transcripts.get(0).getGeneId());
        }
    }
    
    public LinkedHashMap<String,Integer> create_rules(List<TranscriptRecord> transcriptrecord,int DELTA)
    {
        LinkedHashMap<String, Integer> map_result = new LinkedHashMap<>();
        List<List<int[]>> list_transcript = new ArrayList<>();
        for (int cpt = 0; cpt < transcriptrecord.size(); cpt++) {
            list_transcript.add(junctionsFromExons(transcriptrecord.get(cpt).getExons()));
        }
        for (int i = 0; i < list_transcript.size() - 1; i++) {
            for (int k = 0; k < list_transcript.get(i).size(); k++) {
                boolean result = false;
                int[] element_test = list_transcript.get(i).get(k);
                for (int ii = 0; ii < list_transcript.size(); ii++) {
                    for (int jj = 0; jj < list_transcript.get(ii).size(); jj++) {
                        if ((Math.abs(element_test[0] - list_transcript.get(ii).get(jj)[0]) <= DELTA) && (Math.abs(element_test[1] - list_transcript.get(ii).get(jj)[1]) <= DELTA) && i != ii) {
                            result = true;
                            list_transcript.get(ii).remove(jj);
                            jj--;
                        }
                    }
                }
                if (result) {
                    list_transcript.get(i).remove(k);
                }
            }
        }
        for (int iii = 0; iii < list_transcript.size(); iii++) {
            for (int j = 0; j < list_transcript.get(iii).size(); j++) {
                String key = list_transcript.get(iii).get(j)[0] + " " + list_transcript.get(iii).get(j)[1];
                map_result.put(key, iii);
            }
        }
        return map_result;
    }
    
    public TranscriptRecord asign_transcript(List<TranscriptRecord> transcriptrecord, LongreadRecord lrr, int DELTA, LinkedHashMap<String,Integer> map_result)
    {
        List<int[]> list_lrr = junctionsFromExons(lrr.getExons());
        for (int j = 0; j < list_lrr.size(); j++) {
            Set<String> keys = map_result.keySet();
            for (String key : keys) {
                String[] split_result = key.split(" ");
                if ((Math.abs(list_lrr.get(j)[0] - Integer.valueOf(split_result[0])) <= DELTA) && (Math.abs(list_lrr.get(j)[1] - Integer.valueOf(split_result[1])) <= DELTA)) {
                    return transcriptrecord.get(map_result.get(key));
                }
            }
        }
        return null;
    }

    public List<int[]> junctionsFromExons(List<int[]> exons) {
        ArrayList localArrayList = new ArrayList();

        for (int i = 1; i < exons.size(); i++) {
            int j = ((int[]) exons.get(i - 1))[1];
            int k = ((int[]) exons.get(i))[0];
            localArrayList.add(new int[]{j, k});
        }

        return localArrayList;
    }

    public void print(List<int[]> junctions) {
        for (int i = 0; i < junctions.size(); i++) {
            System.out.println(junctions.get(i)[0] + "-" + junctions.get(i)[1]);
        }
    }

    public boolean map(List<int[]> lrr_exons, List<int[]> tr_exons, int DELTA, boolean SOFT)
    {
        boolean bool = true;
	if(SOFT){
            if(tr_exons.size() <= lrr_exons.size()){
	   	for(int i=0; i<tr_exons.size(); i++){
                    if(!isIn((int[])tr_exons.get(i), lrr_exons, DELTA))
	        	bool = false;
	       }
            }
	    else{ bool = false; }
	}
        else{
            if(tr_exons.size() > 0 && tr_exons.size() == lrr_exons.size()) {
                for (int i = 0; i < tr_exons.size(); i++) {
                    if (!isIn((int[]) tr_exons.get(i), lrr_exons, DELTA)) {
                        bool = false;
                        //System.out.println(tr_exons.get(i)[0] +"-" + tr_exons.get(i)[1] + " not in read !");
                    }
                }
            } else {
                bool = false;
                //System.out.println("not same exons size");
            }
        }
        
        //System.out.println(tr_exons.size() + "\t" + lrr_exons.size() + "-->"+bool);
        return bool;
    }

    public boolean isAlreadyIn(int[] paramArrayOfInt, List<int[]> paramList) {
        boolean bool = false;
        for (int[] arrayOfInt : paramList) {
            if ((arrayOfInt[0] == paramArrayOfInt[0]) && (arrayOfInt[1] == paramArrayOfInt[1])) {
                bool = true;
            }
        }
        return bool;
    }

    public boolean isIn(int[] paramArrayOfInt, List<int[]> paramList, int paramInt) {
        boolean bool = false;
        for (int[] arrayOfInt : paramList) {
            if ((Math.abs(arrayOfInt[0] - paramArrayOfInt[0]) <= paramInt) && (Math.abs(arrayOfInt[1] - paramArrayOfInt[1]) <= paramInt)) {
                bool = true;
            }
        }
        return bool;
    }

    private static int[] toIntArray(String paramString) throws NumberFormatException {
        paramString = StringUtils.stripEnd(paramString, ",");
        String[] arrayOfString = paramString.split(",");
        int[] arrayOfInt = new int[arrayOfString.length];

        for (int i = 0; i < arrayOfString.length; i++) {
            arrayOfInt[i] = Integer.valueOf(arrayOfString[i]).intValue();
        }
        return arrayOfInt;
    }

    public String getSizesOfSequence(List<DNASequence> lst) {
        String t = "";
        for (int i = 0; i < lst.size(); i++) {
            t += new Integer(((DNASequence) lst.get(i)).getSequenceAsString().length()) + ",";
        }
        return t;
    }
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
    
    public void callConsensus(File OUTPUT, int nThreads)
    {
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
