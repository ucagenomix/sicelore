package org.ipmc.sicelore.utils;

/**
 *
 * @author kevin lebrigand
 *
 */
import java.util.*;
import java.io.*;
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

    private int nomatch = 0;
    private int onematch = 0;
    private int ambiguous = 0;
    private int monoexon = 0;
    private int multimatchset = 0;
    
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
    
    public void setIsoforms(UCSCRefFlatParser model, int DELTA, String METHOD, boolean AMBIGUOUS_ASSIGN)
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
            
            if("SOFT".equals(METHOD))
                this.setIsoformSoft(molecule, transcripts, DELTA, AMBIGUOUS_ASSIGN);
            else if("STRICT".equals(METHOD))
                this.setIsoformStrict(molecule, transcripts, DELTA, AMBIGUOUS_ASSIGN);
            
            // record womewhere what happened during this setIsoform
            
            compteur++;
            if(compteur%200000 == 0)
                log.info(new Object[]{"\tSetIsoforms\t\t" + compteur + "/" + this.mapMolecules.size()});
        }
        
        log.info(new Object[]{"\tSetIsoforms\t\tend..."});
        log.info(new Object[]{"\tSetIsoforms\t\tmonoexon\t[" + this.monoexon + "]"});
        log.info(new Object[]{"\tSetIsoforms\t\tno match\t[" + this.nomatch + "]"});
        log.info(new Object[]{"\tSetIsoforms\t\tone match\t[" + this.onematch + "]"});
        //log.info(new Object[]{"\tSetIsoforms\t\t>1 match set\t[" + this.multimatchset + "]"});
        log.info(new Object[]{"\tSetIsoforms\t\tambiguous\t[" + this.ambiguous + "]"});

        
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

    public void setIsoformStrict(Molecule molecule, List<TranscriptRecord> transcripts, int DELTA, boolean AMBIGUOUS_ASSIGN)
    {
        List list1=null;
        
        HashMap<String, Integer> candidates = new HashMap<String, Integer>();
        List<Longread> longreads = molecule.getLongreads();
        
        for(Longread lr : longreads){
            List<LongreadRecord> records = lr.getLongreadrecords();

            for(LongreadRecord lrr : records){
                List list = junctionsFromExons(lrr.getExons());

                for(TranscriptRecord transcriptrecord : transcripts){
                    list1 = junctionsFromExons(transcriptrecord.getExons());

                    if(map(list, list1, DELTA)) {
                        String key = transcriptrecord.getTranscriptId() + "|" + transcriptrecord.getGeneId();

                        if(candidates.containsKey(key))
                            candidates.put(key, candidates.get(key) + 1);
                        else
                            candidates.put(key, 1);
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
            
            if(bestCandidates.size() == 1){
                this.onematch++;
                String g = (String)bestCandidates.iterator().next();
                molecule.setTranscriptId(g.split("\\|")[0]);
                molecule.setGeneId(g.split("\\|")[1]);
            }
            else{
                this.ambiguous++;
                int index = new Random().nextInt(bestCandidates.size());
                String g = (String)(bestCandidates.toArray()[index]);
                
                if(AMBIGUOUS_ASSIGN){
                    // set isoform randomly if AMBIGUOUS_ASSIGN=true
                    molecule.setTranscriptId(g.split("\\|")[0]);
                    molecule.setGeneId(g.split("\\|")[1]);
                }
                else{
                    // we set the geneId at least
                    molecule.setTranscriptId("undef");
                    molecule.setGeneId(g.split("\\|")[1]);
                }
            }
        }
        else{
            this.nomatch++;
            int index = new Random().nextInt(molecule.getGeneIds().size());
            Iterator<String> iter = molecule.getGeneIds().iterator();
            for (int i = 0; i < index; i++) { iter.next(); }
            molecule.setTranscriptId("undef");
            molecule.setGeneId((String) iter.next());
        }
        
        // only 1 mono-exonic transcript in the model --> we set the isoform
        if(transcripts.size() == 1 && list1.size() == 0){
           this.monoexon++;
           molecule.setTranscriptId(transcripts.get(0).getTranscriptId());
           molecule.setGeneId(transcripts.get(0).getGeneId());
        }
    }
    
    
    public void setIsoformSoft(Molecule molecule, List<TranscriptRecord> transcripts, int DELTA, boolean AMBIGUOUS_ASSIGN)
    {
        HashMap<Junction, HashSet<TranscriptRecord>> mapper = new HashMap<Junction, HashSet<TranscriptRecord>>();
        List<Junction> junc = new ArrayList<Junction>();
        
        for(TranscriptRecord transcriptrecord : transcripts){
            junc = transcriptrecord.getJunctions();
            for(int i=0; i<junc.size(); i++){
                if(! mapper.containsKey(junc.get(i)))
                    mapper.put(junc.get(i), new HashSet<TranscriptRecord>());
                
                mapper.get(junc.get(i)).add(transcriptrecord);
            }
        }
        
        //System.out.println(molecule.getBarcode() + "|" + molecule.getUmi() + " --> " + molecule.getGeneIds());
        
        /*
        Set cles = mapper.keySet();
        Iterator it = cles.iterator();
        while(it.hasNext()){
           Junction j = (Junction)it.next();
           System.out.println(j);
           HashSet<TranscriptRecord> rec = mapper.get(j);
           Iterator<TranscriptRecord> iterator = rec.iterator();
           while(iterator.hasNext()){
                TranscriptRecord tr = iterator.next();
                System.out.println("\t"+tr.getTranscriptId());
            }
        }
        */
        
        HashSet<TranscriptRecord> h = null;
        HashMap<TranscriptRecord, Integer> scoring = new HashMap<TranscriptRecord, Integer>();
        //HashMap<Integer, Integer> nbJunctions = new HashMap<Integer, Integer>();
        List<Longread> longreads = molecule.getLongreads();
        
        for(Longread lr : longreads){
            List<LongreadRecord> records = lr.getLongreadrecords();
            for(LongreadRecord lrr : records){
                List<Junction> lrr_junc = junctionsListFromExon(lrr.getExons());
                
                // record nb junction of SAM record
                //if(lrr_junc.size() > 0){
                //    if(nbJunctions.containsKey(lrr_junc.size()))
                //        nbJunctions.put(lrr_junc.size(), nbJunctions.get(lrr_junc.size())+1);
                //    else
                //        nbJunctions.put(lrr_junc.size(), 1);
                //}
                
                for(int i=0; i<lrr_junc.size(); i++){
                    Set cles = mapper.keySet();
                    Iterator it = cles.iterator();
                    while(it.hasNext()){
                       Junction j = (Junction)it.next();
                       
                       if(isSameJunction(lrr_junc.get(i), j, DELTA)){
                            HashSet<TranscriptRecord> rec = mapper.get(j);
                            Iterator<TranscriptRecord> iterator = rec.iterator();
                            while(iterator.hasNext()){
                                TranscriptRecord tr = iterator.next();
                                if(scoring.get(tr) != null)
                                    scoring.put(tr, scoring.get(tr) + 1);
                                else
                                    scoring.put(tr, 1);
                            }
                       }
                    }
	       }
            }
        }
        
        if(scoring.isEmpty()){
            // only 1 mono-exonic transcript in the model --> we set the isoform
            if(transcripts.size() == 1 && junc.isEmpty()){
                this.monoexon++;
                molecule.setTranscriptId(transcripts.get(0).getTranscriptId());
                molecule.setGeneId(transcripts.get(0).getGeneId());
            }
            else{
                this.nomatch++;
                // we set the geneId at least
                int index = new Random().nextInt(molecule.getGeneIds().size());
                Iterator<String> iter = molecule.getGeneIds().iterator();
                for (int i = 0; i < index; i++) { iter.next(); }
                molecule.setTranscriptId("undef");
                molecule.setGeneId((String) iter.next());
            }
        }
        else{
            HashMap<TranscriptRecord, Integer> sorted = scoring.entrySet().stream().sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
            .collect(java.util.stream.Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2, LinkedHashMap::new));
        
            /*
            int nb=0;
            Set cles = sorted.keySet();
            Iterator it = cles.iterator();
            while(it.hasNext()){
                TranscriptRecord tr = (TranscriptRecord)it.next();
                if(nb++ < 3)
                    System.out.println(tr + " " + sorted.get(tr));
            }
            */

            HashSet<TranscriptRecord> bestCandidates = new HashSet<TranscriptRecord>();
            int maxValueInMap=(Collections.max(sorted.values()));
            for(Map.Entry<TranscriptRecord, Integer> entry : sorted.entrySet()) 
                if (entry.getValue() == maxValueInMap) { bestCandidates.add(entry.getKey()); }

            if(bestCandidates.size() == 1){
                this.onematch++;
                TranscriptRecord tr = (TranscriptRecord)bestCandidates.iterator().next();
                molecule.setTranscriptId(tr.getTranscriptId());
                molecule.setGeneId(tr.getGeneId());
            }
            else{
                // choose on the number of junctions of transcript and number of junctions of reads
                //int nbJunctionsMaxValue=Collections.max(nbJunctions.entrySet(), Comparator.comparingInt(Map.Entry::getValue)).getKey();
                //HashSet<TranscriptRecord> newCandidates = new HashSet<TranscriptRecord>();
                
                //Iterator<TranscriptRecord> iterator = bestCandidates.iterator();
                //while(iterator.hasNext()){
                //    TranscriptRecord tr = iterator.next();
                //    if(tr.getJunctions().size() == nbJunctionsMaxValue)
                //        newCandidates.add(tr);
                //}
                
                //if(newCandidates.size() == 1){
                //    TranscriptRecord tr = (TranscriptRecord)newCandidates.iterator().next();
                //    molecule.setTranscriptId(tr.getTranscriptId());
                //    molecule.setGeneId(tr.getGeneId());
                //    this.onematch++;
                //}
                //else{
                    this.ambiguous++;
                    int index = new Random().nextInt(bestCandidates.size());
                    TranscriptRecord tr = (TranscriptRecord)(bestCandidates.toArray()[index]);
                    
                    if(AMBIGUOUS_ASSIGN){
                        // set isoform randomly if AMBIGUOUS_ASSIGN=true
                        molecule.setTranscriptId(tr.getTranscriptId());
                        molecule.setGeneId(tr.getGeneId());
                    }
                    else{
                        // we set the geneId at least
                        molecule.setTranscriptId("undef");
                        molecule.setGeneId(tr.getGeneId());
                    }
                //}
            }
        }
    }
    
    public boolean isSameJunction(Junction j1, Junction j2, int DELTA)
    {
        boolean bool = false;
        if ((Math.abs(j1.getStart() - j2.getStart()) <= DELTA) && (Math.abs(j1.getEnd() - j2.getEnd()) <= DELTA))
            bool = true;
        
        return bool;
    }

    public List<int[]> junctionsFromExons(List<int[]> exons) {
        ArrayList lst = new ArrayList();

        for (int i = 1; i < exons.size(); i++) {
            int j = ((int[]) exons.get(i - 1))[1];
            int k = ((int[]) exons.get(i))[0];
            lst.add(new int[]{j, k});
        }

        return lst;
    }
    
    public List<Junction> junctionsListFromExon(List<int[]> exons) {
        ArrayList lst = new ArrayList();

        for (int i = 1; i < exons.size(); i++) {
            int j = ((int[]) exons.get(i - 1))[1];
            int k = ((int[]) exons.get(i))[0];
            lst.add(new Junction(j, k));
        }

        return lst;
    }

    public void print(List<int[]> junctions) {
        for (int i = 0; i < junctions.size(); i++) {
            System.out.println(junctions.get(i)[0] + "-" + junctions.get(i)[1]);
        }
    }

    public boolean map(List<int[]> lrr_exons, List<int[]> tr_exons, int DELTA)
    {
        boolean bool = true;
	if(tr_exons.size() > 0 && tr_exons.size() == lrr_exons.size()) {
            for (int i = 0; i < tr_exons.size(); i++) {
                if (!isIn((int[]) tr_exons.get(i), lrr_exons, DELTA)) {
                    bool = false;
                    //System.out.println(tr_exons.get(i)[0] +"-" + tr_exons.get(i)[1] + " not in read !");
                }
            }
        }
        else
            bool = false;
        
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
                oneNanoporeReadexecutor.awaitTermination(1, TimeUnit.SECONDS);
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
