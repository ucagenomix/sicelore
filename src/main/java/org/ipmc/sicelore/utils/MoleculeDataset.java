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
import java.util.Map.Entry;
import java.util.stream.Collectors;
import gnu.trove.THashMap;
import gnu.trove.THashSet;

public class MoleculeDataset {

    private final Log log;

    THashMap<String, Molecule> mapMolecules;
    THashMap<String, List<Molecule>> mapGenes;
    THashMap<String, Consensus> mapConsensus;

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
    
    private UCSCRefFlatParser model;
    private THashMap<String, Integer> candidates;
    private List<TranscriptRecord> transcripts;
    private THashSet<String> bestCandidates;
    
    public MoleculeDataset()
    {
        log = Log.getInstance(MoleculeDataset.class);
    }
    
    public MoleculeDataset(LongreadParser bam)
    {
        log = Log.getInstance(MoleculeDataset.class);
        log.info(new Object[]{"\tMoleculeDataset init start..."});
        
        Molecule molecule = null;
        this.mapMolecules = new THashMap<String, Molecule>();
        this.mapGenes = new THashMap<String, List<Molecule>>();

        THashMap<String, Longread> mapLongreads = bam.getMapLongreads();
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
                this.mapMolecules.put(molkey, new Molecule(lr.getBarcode(), lr.getUmi(), lr.getRn()));
                ((Molecule) this.mapMolecules.get(molkey)).addLongread(lr);
            }
        }
        log.info(new Object[]{"\tTotal molecules\t\t" + this.mapMolecules.size()});
        
        int multiIG = 0;
        int totalReads = 0;
        //int totalRecords = 0;
        cles = this.mapMolecules.keySet();
        it = cles.iterator();
        while (it.hasNext()){
            molecule = (Molecule) this.mapMolecules.get((String) it.next());
            totalReads += molecule.getLongreads().size();
            if(molecule.getGeneIds().size() > 1)
                multiIG++;
        }        
        log.info(new Object[]{"\tTotal molecule reads\t" + totalReads});
        log.info(new Object[]{"\tTotal molecule multiIG\t" + multiIG});
    }
    
    public List<Molecule> select(String mygene){ return (List<Molecule>) mapGenes.get(mygene); }
    public THashMap<String, Molecule> getMapMolecules() { return this.mapMolecules; }
    
    public void initModel(File REFFLAT)
    {
        this.model = new UCSCRefFlatParser(REFFLAT);
    }
    
    public Molecule getMolecule(String isokey)
    {
        return this.mapMolecules.get(isokey);
    }
    
    public UCSCRefFlatParser getModel()
    {
        return this.model;
    }
    
    public void setIsoforms(int DELTA, String METHOD, boolean AMBIGUOUS_ASSIGN)
    {
        Molecule molecule = null;
        
        int compteur=0;
        log.info(new Object[]{"\tSetIsoforms\t\tstart..."});
        Set cles = this.mapMolecules.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            molecule = (Molecule) this.mapMolecules.get(molkey);

            //List<TranscriptRecord> transcripts = this.model.select(molecule.getGeneIdsArray());
            // sort all selected isoforms on exons number (max to min)
            //Collections.sort(transcripts);
            
            //if(transcripts.size() > 0){
            //    if("SCORE".equals(METHOD))
            //        this.setIsoformScore(molecule, transcripts, DELTA, AMBIGUOUS_ASSIGN);
            //    else if("STRICT".equals(METHOD)) 
                    
                    // KL:29/05/2020 -> for RAM optimization purpose
                    this.setIsoformStrictNew(molecule, DELTA);
            //}
            
            compteur++;
            if(compteur%1000000 == 0)
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
    
    // 29/05/2020
    public void setIsoformStrictNew(Molecule molecule, int DELTA)
    {
        boolean debug = false;
        this.candidates = new THashMap<String, Integer>();
        this.transcripts = this.model.select(molecule.getGeneIdsArray());
        
        if(debug) { System.out.println("transcripts:" + transcripts); }
        if(debug) { System.out.println("molecule:" + molecule.getBarcode() + "\t" + molecule.getUmi() + "\t" + molecule.getGeneIds().toString()); }
        
        // only 1 mono-exonic transcript in the model --> we do set the isoform
        if(transcripts.size() == 1 && transcripts.get(0).getJunctions().isEmpty()){
           this.monoexon++;
           molecule.setTranscriptId(transcripts.get(0).getTranscriptId());
           molecule.setGeneId(transcripts.get(0).getGeneId());
           molecule.setSupporting_reads(1);
           if(debug) { System.out.println("mono-exonic --> transcript_id/gene_id:" + transcripts.get(0).getTranscriptId()+"|"+transcripts.get(0).getGeneId()); } 
        }
        // we need to define which is the best transcriptrecord
        else{
            List<Longread> longreads = molecule.getLongreads();
            for(Longread lr : longreads){
                List<LongreadRecord> records = lr.getLongreadrecords();

                for(LongreadRecord lrr : records){
                    List<Junction> list = lrr.getJunctions();

                    if(debug) { System.out.println("o "+lrr.getName() + "("+list.size()+" junctions)"); }

                    for(TranscriptRecord transcriptrecord : transcripts){
                        List<Junction> list1 = transcriptrecord.getJunctions();

                        if(debug) { System.out.println("\t"+transcriptrecord.getTranscriptId() + "|" + transcriptrecord.getGeneId() + "("+list1.size()+" junctions)"); }

                        if(map(list, list1, DELTA, molecule)){

                            if(debug) { System.out.println("\tmatch"); }

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

                if(debug) { System.out.println("we have candidate(s):" + candidates.size()); }

                this.bestCandidates = new THashSet<String>();
                int maxValueInMap=(Collections.max(candidates.values()));  // This will return max value in the Hashmap
                for(Map.Entry<String, Integer> entry : candidates.entrySet()) {  // Iterate through hashmap
                    if (entry.getValue() == maxValueInMap) { bestCandidates.add(entry.getKey()); }
                }

                if(bestCandidates.size() == 1){
                    this.onematch++;
                    String g = (String)bestCandidates.iterator().next();
                    molecule.setTranscriptId(g.split("\\|")[0]);
                    molecule.setGeneId(g.split("\\|")[1]);
                    molecule.setSupporting_reads(candidates.get(g));

                    if(debug) { System.out.println("only one best candidate --> transcript_id/gene_id:" + g); }
                }
                else{
                    // ambiguous is true, mutiple isoforms are valid
                    // but in STICT mode we have all exons and we set to one of the possible isoform
                    // need to solve Gapsh case where competing with pseudogenes
                    // get the transcripts with the more exons for instance
                    // but need to records the transcriptRecords before
                    this.ambiguous++;
                    int index = new Random().nextInt(bestCandidates.size());
                    String g = (String)(bestCandidates.toArray()[index]);
                    molecule.setTranscriptId(g.split("\\|")[0]);
                    molecule.setGeneId(g.split("\\|")[1]);
                    molecule.setSupporting_reads(candidates.get(g));
                    // we set the geneId at least
                    //molecule.setTranscriptId("undef");
                    //molecule.setGeneId(g.split("\\|")[1]);
                }
            }
            // no condaidates and several transcripts
            else if(transcripts.size() > 0){
                if(debug) { System.out.println("no candidate: choose beteen " + transcripts.size() + " transcripts -> get the most complex one"); }
                // KL 22/04/2020
                // case multiGene here, like Pkm / RP23-320D23.6
                // set to the more complex gene, the one with the more isoforms in List<TranscriptRecord> transcripts
                // need some bug fix here in the future verson !!!!
                this.nomatch++;
                //int index = new Random().nextInt(molecule.getGeneIds().size());
                //Iterator<String> iter = molecule.getGeneIds().iterator();
                //for (int i = 0; i < index; i++) { iter.next(); }
                molecule.setTranscriptId("undef");
                //molecule.setGeneId((String) iter.next());
                molecule.setGeneId(getGeneIdForMostComplexTranscript());
                if(debug) { System.out.println("no candidate --> gene_id:" + molecule.getGeneId()); }
            }
            //else{ // we have a problem here no transcripts for genes (GE:tag) in model

                //System.out.println(transcripts.size() + "\t" + molecule.getBarcode()+ ":"  + molecule.getUmi() + "\t" + molecule.getGeneIdsArray());

            //}
        }
    }
    
    public String getGeneIdForMostComplexTranscript()
    {
        String mostRepeatedWord = "";
        try{
            ArrayList<String> list = new ArrayList<>();
            for(TranscriptRecord transcriptrecord : transcripts)
                list.add(transcriptrecord.getGeneId());

            mostRepeatedWord  = list.stream()
              .collect(Collectors.groupingBy(w -> w, Collectors.counting()))
              .entrySet()
              .stream()
              .max(Comparator.comparing(Entry::getValue))
              .get()
              .getKey();
            
            list = null;
        }
        catch(Exception e){ e.printStackTrace(); }
        
        return mostRepeatedWord;
    }
    /*
    public void setIsoformStrict(Molecule molecule, List<TranscriptRecord> transcripts, int DELTA)
    {
        List list1=null;
        boolean debug = false;
        
        if(debug) { System.out.println("transcripts:" + transcripts); }
        if(debug) { System.out.println("molecule:" + molecule.getBarcode() + "\t" + molecule.getUmi() + "\t" + molecule.getGeneIds().toString()); }
        
        THashMap<String, Integer> candidates = new THashMap<String, Integer>();
        List<Longread> longreads = molecule.getLongreads();

        for(Longread lr : longreads){
            List<LongreadRecord> records = lr.getLongreadrecords();
            
            for(LongreadRecord lrr : records){
                List list = junctionsFromExons(lrr.getExons());
                
                if(debug) { System.out.println("o "+lrr.getName() + "("+list.size()+" junctions)"); }
                
                for(TranscriptRecord transcriptrecord : transcripts){
                    list1 = junctionsFromExons(transcriptrecord.getExons());
                    
                    if(debug) { System.out.println("\t"+transcriptrecord.getTranscriptId() + "|" + transcriptrecord.getGeneId() + "("+list1.size()+" junctions)"); }
                    
                    if(map(list, list1, DELTA, molecule)) {
                        
                        if(debug) { System.out.println("\tmatch"); }
                        
                        String key = transcriptrecord.getTranscriptId() + "|" + transcriptrecord.getGeneId();

                        if(candidates.containsKey(key))
                            candidates.put(key, candidates.get(key) + 1);
                        else
                            candidates.put(key, 1);
                    }
                }
            }
        }
        
        if(debug) { System.out.println("candidates:" + candidates); }
        
        if(candidates.size() > 0){
            
            if(debug) { System.out.println("we have candidate(s):" + candidates.size()); }
            
            THashSet<String> bestCandidates = new THashSet<String>();
            int maxValueInMap=(Collections.max(candidates.values()));  // This will return max value in the Hashmap
            for(Map.Entry<String, Integer> entry : candidates.entrySet()) {  // Iterate through hashmap
                if (entry.getValue() == maxValueInMap) { bestCandidates.add(entry.getKey()); }
            }
            
            if(bestCandidates.size() == 1){
                this.onematch++;
                String g = (String)bestCandidates.iterator().next();
                molecule.setTranscriptId(g.split("\\|")[0]);
                molecule.setGeneId(g.split("\\|")[1]);
                molecule.setSupporting_reads(candidates.get(g));
                
                if(debug) { System.out.println("only one best candidate --> transcript_id/gene_id:" + g); }
            }
            else{
                // ambiguous is true, mutiple isoforms are valid
                // but in STICT mode we have all exons and we set to one of the possible isoform
                // need to solve Gapsh case where competing with pseudogenes
                // get the transcripts with the more exons for instance
                // but need to records the transcriptRecords before
                this.ambiguous++;
                int index = new Random().nextInt(bestCandidates.size());
                String g = (String)(bestCandidates.toArray()[index]);
                molecule.setTranscriptId(g.split("\\|")[0]);
                molecule.setGeneId(g.split("\\|")[1]);
                molecule.setSupporting_reads(candidates.get(g));
                // we set the geneId at least
                //molecule.setTranscriptId("undef");
                //molecule.setGeneId(g.split("\\|")[1]);
            }
        }
        else{
            // KL 22/04/2020
            // case multiGene here, like Pkm / RP23-320D23.6
            // set to the more complex gene, the one with the more isoforms in List<TranscriptRecord> transcripts
            // need some bug fix here in the future verson !!!!
            this.nomatch++;
            String geneIdChoice = this.getMostComplexGene(transcripts);
            //int index = new Random().nextInt(molecule.getGeneIds().size());
            //Iterator<String> iter = molecule.getGeneIds().iterator();
            //for (int i = 0; i < index; i++) { iter.next(); }
            molecule.setTranscriptId("undef");
            //molecule.setGeneId((String) iter.next());
            molecule.setGeneId(geneIdChoice);
            if(debug) { System.out.println("no candidate --> gene_id:" + molecule.getGeneId()); }
        }
        
        // only 1 mono-exonic transcript in the model --> we do set the isoform
        if(transcripts.size() == 1 && list1.isEmpty()){
           this.monoexon++;
           molecule.setTranscriptId(transcripts.get(0).getTranscriptId());
           molecule.setGeneId(transcripts.get(0).getGeneId());
           molecule.setSupporting_reads(1);
           if(debug) { System.out.println("mono-exonic --> transcript_id/gene_id:" + transcripts.get(0).getTranscriptId()+"|"+transcripts.get(0).getGeneId()); } 
        }
    }

    public void setIsoformScore(Molecule molecule, List<TranscriptRecord> transcripts, int DELTA, boolean AMBIGUOUS_ASSIGN)
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
        
        HashSet<TranscriptRecord> h = null;
        HashMap<TranscriptRecord, Integer> scoring = new HashMap<TranscriptRecord, Integer>();
        //HashMap<Integer, Integer> nbJunctions = new HashMap<Integer, Integer>();
        List<Longread> longreads = molecule.getLongreads();
        
        for(Longread lr : longreads){
            List<LongreadRecord> records = lr.getLongreadrecords();
            for(LongreadRecord lrr : records){
                List<Junction> lrr_junc = junctionsListFromExon(lrr.getExons());
                
                for(int i=0; i<lrr_junc.size(); i++){
                    Set cles = mapper.keySet();
                    Iterator it = cles.iterator();
                    while(it.hasNext()){
                       Junction j = (Junction)it.next();
                       
                       if(isSameJunction(lrr_junc.get(i), j, DELTA)){
                            
                            // add junction to molecule
                            molecule.addJunction(j.getIntArray());
                            
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
            }
        }
    }
    */
    
    public void setFusions()
    {
        Molecule molecule = null;
        Map<String, Integer> count = new HashMap<String, Integer>();
        
        log.info(new Object[]{"\tSetFusions\t\tstart..."});
        Set cles = this.mapMolecules.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            molecule = (Molecule) this.mapMolecules.get(molkey);
            
            if(molecule.getUmi() != null && molecule.getGeneIds().size() > 1){
                //System.out.println(molecule.toString());
                String key = molecule.getGeneIds().toString();
                int c = count.containsKey(key) ? count.get(key) : 0;
                count.put(key, c + 1);
            }
        }
        
        HashMap<String, Integer> sorted = count.entrySet().stream().sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
            .collect(java.util.stream.Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2, LinkedHashMap::new));
        
        cles = sorted.keySet();
        it = cles.iterator();
        while (it.hasNext()) {
            String key = (String) it.next();
            int c = (Integer)sorted.get(key);
            if(c >= 5)
                log.info(new Object[]{"\t" + c + " distincts molecules support " + key + " genes fusion"});
        }
    }
    
    public boolean isSameJunction(Junction j1, Junction j2, int DELTA)
    {
        boolean bool = false;
        if ((Math.abs(j1.getStart() - j2.getStart()) <= DELTA) && (Math.abs(j1.getEnd() - j2.getEnd()) <= DELTA))
            bool = true;
        
        return bool;
    }
    /*
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
    */
    public void print(List<int[]> junctions) {
        for (int i = 0; i < junctions.size(); i++) {
            System.out.println(junctions.get(i)[0] + "-" + junctions.get(i)[1]);
        }
    }

    public boolean map(List<Junction> juncRead, List<Junction> juncRef, int DELTA, Molecule molecule)
    {   
        boolean bool = true;
	if(juncRef.size() > 0 && juncRef.size() == juncRead.size()) {
            for (int i = 0; i < juncRef.size(); i++) {
                if (!isIn((Junction) juncRef.get(i), juncRead, DELTA))
                    bool = false;
            }
        }
        else
            bool = false;
        
        // get all junctions of all reads/molecules
        for (int i = 0; i < juncRef.size(); i++) {
            if (isIn((Junction) juncRef.get(i), juncRead, DELTA))
                 molecule.addJunction(juncRef.get(i));
        }
        
        return bool;
    }

    public boolean isIn(Junction j, List<Junction> lst, int DELTA)
    {
        boolean bool = false;
        for(Junction jlist : lst) {
            if ((Math.abs(jlist.getStart() - j.getStart()) <= DELTA) && (Math.abs(jlist.getEnd() - j.getEnd()) <= DELTA)) {
                bool = true;
            }
        }
        return bool;
    }
    
    public Matrix produceMatrix(HashSet<String> authorizedCells)
    {
        int nb = 0;
        Matrix matrix = new Matrix(authorizedCells);
        THashMap<String, List<TranscriptRecord>> mapGenesTranscripts = model.getMapGenesTranscripts();

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
        log.info(new Object[]{"\tCalling consensus start with " + nThreads + " threads"});
        
        this.mapConsensus = new THashMap<String, Consensus>();
        for(String key : this.mapMolecules.keySet()) {
            Molecule m = (Molecule) this.mapMolecules.get(key);
            String mykey = m.getBarcode() + "-" + m.getUmi() + "-" + m.getLongreads().size();
            this.mapConsensus.put(mykey, new Consensus(mykey, m.getLongreads(), true));
        }
        
        log.info(new Object[]{"\tTotal consensus to compute\t" + this.mapConsensus.size()});
        
        future_list = new ConcurrentLinkedDeque<Future<String>>();
        oneNanoporeReadexecutor = MoreExecutors.listeningDecorator(Executors.newFixedThreadPool(nThreads));

        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));

            Set cles = this.mapConsensus.keySet();
            itglobal = cles.iterator();
            int i = 0;
            while (itglobal.hasNext() && i++ < (5 * nThreads)) {
                String key = (String) itglobal.next();
                Consensus consensus = (Consensus) this.mapConsensus.get(key);
                ListenableFuture<String> submit = oneNanoporeReadexecutor.submit(consensus);
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
                        Consensus consensus = (Consensus) this.mapConsensus.get(key);
                        ListenableFuture<String> submit = oneNanoporeReadexecutor.submit(consensus);
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

    private synchronized void write(String rslt)
    {
        try {
            os.writeBytes(rslt);
            // thats a mess.... but it works, need to change ket of 
            // this.mapMolecules ligne 68, does this will break something else ?
            String[] tmp = rslt.split("\n");
            tmp[0] = tmp[0].replaceAll("@","");
            String[] tmp2 = tmp[0].split("-");
            this.mapMolecules.get(tmp2[0]+":"+tmp2[1]).setConsensus(tmp[1].getBytes());
        } catch (Exception e) { e.printStackTrace(); }
    }
}
