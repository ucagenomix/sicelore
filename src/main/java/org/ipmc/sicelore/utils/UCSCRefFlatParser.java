package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */ 
import com.google.common.util.concurrent.ListenableFuture;
import com.google.common.util.concurrent.ListeningExecutorService;
import com.google.common.util.concurrent.MoreExecutors;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import java.util.*;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.annotation.Strand;
import java.util.concurrent.ConcurrentLinkedDeque;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;
import org.biojava.nbio.core.util.ConcurrencyTools;
import gnu.trove.THashMap;

public class UCSCRefFlatParser implements GeneModelParser {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;
    
    public THashMap<String, List<TranscriptRecord>> mapGenesTranscripts;
    public THashMap<String, Consensus> mapConsensus;

    public UCSCRefFlatParser refmodel;
    public int DELTA = 2;
    public int MINEVIDENCE = 5;
    public int RNMIN = 1;
    public int NOVELINDEX=1;
               
    private ListeningExecutorService oneNanoporeReadexecutor;
    private Deque<Future<String>> future_list;
    private DataOutputStream os;
    private Iterator itglobal;
    
    public UCSCRefFlatParser(File refFlat)
    {
        log = Log.getInstance(UCSCRefFlatParser.class);
        log.info(new Object[]{"UCSCRefFlatParser\tstart..."});

        this.mapGenesTranscripts = new THashMap<String, List<TranscriptRecord>>();
        int nb = 0;

        try {
            FileInputStream in = new FileInputStream(refFlat);
            BufferedReader reader = new BufferedReader(new InputStreamReader(in));
            String line = null;
            while ((line = reader.readLine()) != null) {
                TranscriptRecord f = parseLine(line);
                if (f != null) {
                    nb++;
                    if (this.mapGenesTranscripts.containsKey(f.getGeneId())) {
                        ((ArrayList<TranscriptRecord>) this.mapGenesTranscripts.get(f.getGeneId())).add(f);
                    } else {
                        List<TranscriptRecord> l = new ArrayList<TranscriptRecord>();
                        l.add(f);
                        this.mapGenesTranscripts.put(f.getGeneId(), l);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        log.info(new Object[]{"UCSCRefFlatParser\tend..."});
        log.info(new Object[]{"Number of Genes Symbols\t[" + this.mapGenesTranscripts.size() + "]"});
        log.info(new Object[]{"Number of Transcripts\t[" + nb + "]"});
    }
    
    public UCSCRefFlatParser(int DELTA, int MINEVIDENCE, int RNMIN, UCSCRefFlatParser refmodel)
    {
        log = Log.getInstance(UCSCRefFlatParser.class);
        this.mapGenesTranscripts = new THashMap<String, List<TranscriptRecord>>();
        
        this.refmodel=refmodel;
        this.DELTA =DELTA;
        this.MINEVIDENCE=MINEVIDENCE;
        this.RNMIN=RNMIN;
        this.NOVELINDEX=1;
    }

    public TranscriptRecord parseLine(String line) throws GTFParseException {
        String[] fields = line.split("\t");
        TranscriptRecord record = TranscriptRecord.fromRefFlat(fields);
        
        //if (record.getChrom().contains("_")) {
        //    return null;
        //}
        if (record.getExonBases() == 0) {
            return null;
        }

        return record;
    }

    public List<TranscriptRecord> select(String[] filter) {
        List<TranscriptRecord> filteredList = new ArrayList<TranscriptRecord>();
        for (String f : filter) {
            if (this.mapGenesTranscripts.containsKey(f)) {
                for (TranscriptRecord t : this.mapGenesTranscripts.get(f)) {
                    filteredList.add(t);
                }
            }
        }
        return filteredList;
    }

    public TranscriptRecord select(String geneId, String transcriptId)
    {
        TranscriptRecord tr = null;
        List<TranscriptRecord> list = select(new String[] {geneId});
        ListIterator<TranscriptRecord> it = list.listIterator();
        while(it.hasNext()){
            TranscriptRecord t = (TranscriptRecord)it.next();
            if(t.getTranscriptId().equals(transcriptId))
                tr = t;
        }
        return tr;
    }

    public THashMap<String, List<TranscriptRecord>> getMapGenesTranscripts() {
        return mapGenesTranscripts;
    }
    
    public void loader(File INPUT, CellList cellList, String CELLTAG, String UMITAG, String GENETAG, String ISOFORMTAG, String RNTAG)
    {
        pl = new htsjdk.samtools.util.ProgressLogger(log, 1000000, "Processed\t", "Records");
        
        log.info(new Object[]{String.format("Loader Bam Start...")});
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
        htsjdk.samtools.SAMSequenceDictionary dictionnary = samFileHeader.getSequenceDictionary();
        
        try{
            for(SAMSequenceRecord x : dictionnary.getSequences()){
                SAMRecordIterator iter = samReader.query(x.getSequenceName(), 1, x.getSequenceLength(), false);
                //log.info(new Object[]{"\tProcessing ref. " + x.getSequenceName() + "\t[" + mymodel.getMapGenesTranscripts().size() +" genes]"});
                
                while(iter.hasNext()){
                    SAMRecord r = iter.next();
                    pl.record(r);
                    
                    String BC = (String)r.getAttribute(CELLTAG);
                    String U8 = (String)r.getAttribute(UMITAG);
                    String IG = (String)r.getAttribute(GENETAG);
                    String IT = (String)r.getAttribute(ISOFORMTAG);
                    int    RN = ((Integer) r.getAttribute(RNTAG) != null) ? (Integer) r.getAttribute(RNTAG) : 0;
                    
                    LongreadRecord lrr = LongreadRecord.fromSAMRecord(r, true);
                    
                    // never null case if umifound or isobam from isoformMatrix pipeline used
                    // but we filter out some not reliable reads just as in isoformMatrix
                    if(lrr != null && lrr.getMapqv() > 0 && !lrr.getIsChimeria() && !lrr.getIsReversed() && RN >= this.RNMIN && cellList.contains(BC)){
                        
                        // we have a geneId
                        if(! "undef".equals(IG)){
                            TranscriptRecord tr = new TranscriptRecord(IG, IT);
                                
                            // never seen this gene, create a new List<TranscriptRecord> including IT of lrr
                            if(! mapGenesTranscripts.containsKey(IG)) {
                                if(! "undef".equals(IT))
                                    tr = refmodel.select(IG, IT);
                                
                                List<TranscriptRecord> lst = new ArrayList<TranscriptRecord>();
                                lst.add(tr);
                                mapGenesTranscripts.put(IG, lst);
                            }
                            // we know this gene
                            else {
                                // but we don't know this transcript
                                if(!(mapGenesTranscripts.get(IG)).contains(tr)){
                                    if(! "undef".equals(IT))
                                        tr = refmodel.select(IG, IT);
                                    
                                    ((ArrayList<TranscriptRecord>) mapGenesTranscripts.get(IG)).add(tr);
                                }
                            }
                            
                            // finally add the molecule to the TranscriptRecord
                            mapGenesTranscripts.get(IG).get(mapGenesTranscripts.get(IG).indexOf(tr)).add(lrr);
                        }
                        // intergenic region, new genes ?
                        else{ }
                    }
                }
                iter.close();
            }
            samReader.close();
        } catch (Exception e) { e.printStackTrace(); }
    }
    
    // parcours de la map de genes and collapse undef categorie
    public void collapser()
    {
        log.info(new Object[]{String.format("Collapser Start...[" + mapGenesTranscripts.size() + " total genes]")});
        
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            
            if(tAll.contains(new TranscriptRecord(geneId, "undef"))){
                
                TranscriptRecord tUndef = tAll.get(tAll.indexOf(new TranscriptRecord(geneId, "undef")));
                List<TranscriptRecord> tNovel = this.collapse(tUndef.getEvidenceList(), geneId);
                
                tAll.remove(tUndef);
                for(int i=0; i<tNovel.size(); i++){
                    if(tNovel.get(i).getEvidenceList().size() >= MINEVIDENCE)
                        tAll.add(tNovel.get(i));
                }
            }
        }
    }
    
    // parcours de la map de genes and intilialize TranscriptRecord attributes based on LongreadRecord stored
    public void initialize()
    {    
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            for(int i=0; i<tAll.size(); i++)
                tAll.get(i).initialize();
        }
    }
    
    // filtration of sub-part degradated isoforms
    public void filter()
    {
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            Collections.sort(tAll);
            
            List<TranscriptRecord> keep = new ArrayList<TranscriptRecord>();
            for(int i=0; i<tAll.size(); i++){
                TranscriptRecord t = tAll.get(i);
                
                // we keep all Ensembl transcripts
                if(t.getIs_known())
                    keep.add(t);
                else{
                    if(! isPartOfLonger(t, keep, refmodel.select(new String[]{t.getGeneId()})))
                        keep.add(t);
                }
            }
            mapGenesTranscripts.put(geneId, keep);
        }
    }
    
    // classification
    public void classifier()
    {
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            for(int i=0; i<tAll.size(); i++){
                TranscriptRecord t = tAll.get(i);
                if(t.getIs_novel())
                    this.noveltyDetector(t, refmodel.select(new String[]{t.getGeneId()}));
            }
        }
    }
    
    // validator
    public void validator(BEDParser cage, BEDParser polyA, File SHORT, int cageCo, int polyaCo, int juncCo)
    {
        int nb=0;
        SamReader samReaderShort = SamReaderFactory.makeDefault().open(SHORT);
        HashMap<String, Integer> isDone = new HashMap<String, Integer>();
        int e_prev=0;
        
        log.info(new Object[]{String.format("Validator Start...[" + mapGenesTranscripts.size() + " total genes]")});
        
        try {
            for(String geneId : mapGenesTranscripts.keySet()) {
                List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
                nb++;
                if(nb%2500 == 0)
                    log.info(new Object[]{String.format(nb+" genes processed")});
                
                for(int i=0; i<tAll.size(); i++){
                    TranscriptRecord t = tAll.get(i);

                    // get distance of isoforms to cage-peak
                    int dist_cage_peak = cage.getDistanceCage(t.getChrom(), t.getStrand(), (t.getStrand() == Strand.POSITIVE)?t.getTxStart():t.getTxEnd());
                    //boolean within_cage_peak = cage.isWithin(t.getChrom(), t.getStrand(), (t.getStrand() == Strand.POSITIVE)?t.getTxStart():t.getTxEnd());
                    // get distance of isoforms to polyA
                    int dist_polya = polyA.getDistancePolyA(t.getChrom(), t.getStrand(), (t.getStrand() == Strand.POSITIVE)?t.getTxEnd():t.getTxStart());
                    //boolean within_polya = polyA.isWithin(t.getChrom(), t.getStrand(), (t.getStrand() == Strand.POSITIVE)?t.getTxEnd():t.getTxStart());

                    t.setDist_cage(dist_cage_peak);
                    t.setIs_valid_cage((Math.abs(dist_cage_peak) <= cageCo)? true: false);
                    t.setDist_polya(dist_polya);
                    t.setIs_valid_polya((Math.abs(dist_polya) <= polyaCo)? true: false);
                    
                    List<int[]> lst = t.getNovelJunctions();
                    boolean is_junction_valid=true;
                    int totalSupportReads=0;
                    for(int index=0; index<lst.size(); index++){
                        int donor = lst.get(index)[0];
                        int acceptor = lst.get(index)[1];
                        String jkey = t.getChrom()+":"+donor+"-"+acceptor;
                        
                        int supportReads=0;
                        int totalReads=0;
                        if(! isDone.containsKey(jkey)){
                            SAMRecordIterator iter = samReaderShort.query(t.getChrom(), donor, acceptor, false);
                            while(iter.hasNext()){
                                SAMRecord r = iter.next();

                                List<int[]> junctions = new ArrayList<int[]>();
                                List<AlignmentBlock> blocks = r.getAlignmentBlocks();
                                for (int b=0; b<blocks.size(); b++) {
                                    AlignmentBlock currBlock = blocks.get(b);
                                    int s = currBlock.getReferenceStart();
                                    int e = s + currBlock.getLength();
                                    if(b>0){ junctions.add(new int[]{e_prev-1, s}); }
                                    e_prev = e;
                                }
                                
                                // for junction validation no DELTA
                                if(this.isIn(new int[]{donor, acceptor}, junctions, 0)){
                                    supportReads++;
                                    totalSupportReads++;
                                }

                                totalReads++;
                            }
                            iter.close();
                            
                            isDone.put(jkey, supportReads);
                        }
                        else
                            totalSupportReads += isDone.get(jkey);
                        
                        if(isDone.get(jkey) < juncCo){
                            is_junction_valid = false;
                        }
                    }
                    
                    t.setIs_valid_junction(is_junction_valid);
                    t.setJunctionReads(totalSupportReads);
                    
                    //System.out.println(t.getIs_valid_cage() + "," + t.getIs_valid_polya() + "," + t.getIs_valid_junction());
                    
                    if(t.getIs_valid_cage() && t.getIs_valid_polya() && t.getIs_valid_junction())
                        t.setIs_valid(true);
                }
            }
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { samReaderShort.close();  } catch (Exception e) { System.err.println("can not close stream"); } }
    }
        
    public boolean isIn(int[] junc, List<int[]> list, int delta)
    {
        boolean bool = false;
        for (int[] a : list) {
            //allcmp.add(junc[0] + "-" + junc[1] + "/"+ a[0] +"-" + a[1]);
            if ((Math.abs(a[0] - junc[0]) <= delta) && (Math.abs(a[1] - junc[1]) <= delta))
                bool = true;
        }
        return bool;
    }
    
    public void noveltyDetector(TranscriptRecord t, List<TranscriptRecord> modelLst)
    {
        // novel_in_catalog
        //  - combination_of_known_junctions
        //  - combination_of_known_splicesites
        // novel_not_in_catalog
        //  - at_least_one_novel_splicesite

        List<int[]> junctions = junctionsFromExons(t.getExons());
        List<int[]> modelJunctions = new ArrayList<int[]>();
        List<Integer> modelSplice = new ArrayList<Integer>();
        
        for(int i=0; i<modelLst.size(); i++){
            TranscriptRecord tMod = modelLst.get(i);
            List<int[]> jMod = junctionsFromExons(tMod.getExons());
            for(int j=0; j<jMod.size(); j++){
                modelJunctions.add(jMod.get(j));
                modelSplice.add(jMod.get(j)[0]);
                modelSplice.add(jMod.get(j)[1]);
            }
        }
        
        for(int i=0; i<junctions.size(); i++){
            if(!isIn(junctions.get(i), modelJunctions, this.DELTA)){
                
                // cette junctions n'est pas connue
                // est-elle formées sur la base des sites de splices connus ?
                // meaning start and stop are in the modelJunctions start/stop
                if(modelSplice.contains(junctions.get(i)[0]) && modelSplice.contains(junctions.get(i)[1])){
                    if("undef".equals(t.getCategorie())){
                        t.setCategorie("novel_in_catalog");
                        t.setSubcategorie("combination_of_known_splicesites");
                    }
                    t.addNovelJunction(junctions.get(i));
                }
                else{
                    t.setCategorie("novel_not_in_catalog");
                    t.setSubcategorie("at_least_one_novel_splicesite");
                    t.addNovelJunction(junctions.get(i));
                }
            }
        }
        
        // si toutes les junctions sont connues
        if("undef".equals(t.getCategorie())){
            t.setCategorie("novel_in_catalog");
            t.setSubcategorie("combination_of_known_junctions");
        }
    }
    
    public boolean isPartOfLonger(TranscriptRecord t, List<TranscriptRecord> lst, List<TranscriptRecord> modelLst)
    {
        // -------------
        // Here we remove degradated but also all intron-retention isoforms 
        // which have all their junctions in full gencode isoforms
        // -------------
        boolean bool = false;
        List<int[]> junctions = junctionsFromExons(t.getExons());
        
        for(int i=0; i<lst.size(); i++){
            TranscriptRecord tkept = lst.get(i);
            List<int[]> jkept = junctionsFromExons(tkept.getExons());
            
            if(isAllInclude(junctions,jkept)){
                bool = true;
                //log.info(new Object[]{t.getTranscriptId() + " is part of " + tkept.getTranscriptId()});
            }
        }
        if(t.getIs_novel()){
            for(int i=0; i<modelLst.size(); i++){
                TranscriptRecord tkept = modelLst.get(i);
                List<int[]> jkept = junctionsFromExons(tkept.getExons());
            
                if(isAllInclude(junctions,jkept)){
                    bool = true;
                    //log.info(new Object[]{t.getTranscriptId() + " is part of " + tkept.getTranscriptId()});
                }
            }
        }
        
        return bool;
    }
    
    public void callConsensus(File OUTPUT, int nThreads)
    {
        log.info(new Object[]{"Calling consensus start with " + nThreads + " threads"});
        
        this.mapConsensus = new THashMap<String, Consensus>();
        for(String geneId : this.mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            for(int i=0; i<tAll.size(); i++){
                TranscriptRecord t = tAll.get(i);
                
                this.mapConsensus.put(t.getTranscriptId(), new Consensus(t.getTranscriptId(), t.getEvidenceList()));
            }
        }
        
        log.info(new Object[]{"Total consensus to compute\t" + this.mapConsensus.size()});
        
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
            } catch (InterruptedException ex) { log.info(new Object[]{"Error:\t[" + ex + "]"}); }

            os.close();
        }
        catch (Exception e) {  e.printStackTrace(); } 
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
    }

    private synchronized void write(String rslt)
    {
        try { os.writeBytes(rslt); } catch (Exception e) { e.printStackTrace(); }
    }
    
    // statistiques
    public void statistics()
    {
        int sum_isoforms = 0;
        int sum_evidences = 0;
        int sum_valid_isoforms = 0;
        int sum_valid_evidences = 0;
        
        log.info(new Object[]{String.format("Printing statistics...")});
        
        HashMap<String, Integer> stats = new HashMap<String, Integer>();
        String[] allkeys = {"undef","undef2","full_splice_match","gencode","novel_in_catalog","novel_not_in_catalog","combination_of_known_junctions","combination_of_known_splicesites","at_least_one_novel_splicesite"};
        for(int i=0; i<allkeys.length; i++){
            stats.put(allkeys[i]+"_count", 0);
            stats.put(allkeys[i]+"_evidences", 0);
            stats.put(allkeys[i]+"_count_valid", 0);
            stats.put(allkeys[i]+"_evidences_valid", 0);
        }
        
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            for(int i=0; i<tAll.size(); i++){
                TranscriptRecord t = tAll.get(i);
                
                //log.info(new Object[]{t.getCategorie() + "," + t.getSubcategorie()});
                stats.put(t.getCategorie()+"_count", stats.get(t.getCategorie()+"_count")+1);
                stats.put(t.getCategorie()+"_evidences", stats.get(t.getCategorie()+"_evidences")+t.getEvidenceList().size());
                stats.put(t.getSubcategorie()+"_count", stats.get(t.getSubcategorie()+"_count")+1);
                stats.put(t.getSubcategorie()+"_evidences", stats.get(t.getSubcategorie()+"_evidences")+t.getEvidenceList().size());
                
                sum_isoforms++;
                sum_evidences += t.getEvidenceList().size();
                
                // les stats des known + valid novel
                if(t.getIs_known() || (t.getIs_novel() && t.getIs_valid())){
                    sum_valid_isoforms++;
                    sum_valid_evidences += t.getEvidenceList().size();
                    
                    stats.put(t.getCategorie()+"_count_valid", stats.get(t.getCategorie()+"_count_valid")+1);
                    stats.put(t.getCategorie()+"_evidences_valid", stats.get(t.getCategorie()+"_evidences_valid")+t.getEvidenceList().size());
                    stats.put(t.getSubcategorie()+"_count_valid", stats.get(t.getSubcategorie()+"_count_valid")+1);
                    stats.put(t.getSubcategorie()+"_evidences_valid", stats.get(t.getSubcategorie()+"_evidences_valid")+t.getEvidenceList().size());
                }
            }
        }
        
        log.info(new Object[]{"-----------------------------------------------------------------------"});
        log.info(new Object[]{String.format("\t\t\t\t\tall_set (UMI)\tvalid_set (UMI)")});
        log.info(new Object[]{String.format("total_genes\t\t\t\t%d",mapGenesTranscripts.size())});
        log.info(new Object[]{String.format("total_isoforms\t\t\t\t%d (%d)\t%d (%d)",sum_isoforms,sum_evidences,sum_valid_isoforms,sum_valid_evidences)});
        log.info(new Object[]{String.format("full_splice_match")});
        log.info(new Object[]{String.format(" o gencode\t\t\t\t%d (%d)\t%d (%d)",stats.get("gencode_count"),stats.get("gencode_evidences"),stats.get("gencode_count_valid"),stats.get("gencode_evidences_valid"))});
        log.info(new Object[]{String.format("novel_in_catalog")});
        log.info(new Object[]{String.format(" o combination_of_known_junctions\t%d (%d)\t%d (%d)",stats.get("combination_of_known_junctions_count"),stats.get("combination_of_known_junctions_evidences"),stats.get("combination_of_known_junctions_count_valid"),stats.get("combination_of_known_junctions_evidences_valid"))});
        log.info(new Object[]{String.format(" o combination_of_known_splicesites\t%d (%d)\t%d (%d)",stats.get("combination_of_known_splicesites_count"),stats.get("combination_of_known_splicesites_evidences"),stats.get("combination_of_known_splicesites_count_valid"),stats.get("combination_of_known_splicesites_evidences_valid"))});
        log.info(new Object[]{String.format("novel_not_in_catalog")});
        log.info(new Object[]{String.format(" o at_least_one_novel_splicesite\t%d (%d)\t%d (%d)",stats.get("at_least_one_novel_splicesite_count"),stats.get("at_least_one_novel_splicesite_evidences"),stats.get("at_least_one_novel_splicesite_count_valid"),stats.get("at_least_one_novel_splicesite_evidences_valid"))});
        log.info(new Object[]{"------------------------------------------------------------------------"});
    }
    
    //export files
    public void exportFiles(File TXT, File FLAT, File FLATVALID,File GFF, File GFFVALID)
    {
        BufferedOutputStream ostxt = null;
        BufferedOutputStream osrefflat = null;
        BufferedOutputStream osrefflatvalid = null;
        BufferedOutputStream osgff = null;
        BufferedOutputStream osgffvalid = null;
        
        try {
            ostxt = new BufferedOutputStream(new java.io.FileOutputStream(TXT));
            osrefflat = new BufferedOutputStream(new java.io.FileOutputStream(FLAT));
            osrefflatvalid = new BufferedOutputStream(new java.io.FileOutputStream(FLATVALID));
            ostxt.write(new TranscriptRecord().printLegendTxt().getBytes());
            osgff = new BufferedOutputStream(new java.io.FileOutputStream(GFF));
            osgffvalid = new BufferedOutputStream(new java.io.FileOutputStream(GFFVALID));
            
            for(String geneId : mapGenesTranscripts.keySet()) {
                // running through all ENST + NOVEL
                List<TranscriptRecord> tList = mapGenesTranscripts.get(geneId);

                for(int i=0; i<tList.size(); i++) {
                    TranscriptRecord t = tList.get(i);
                    
                    osrefflat.write(t.printRefflat().getBytes());
                    ostxt.write(t.printTxt().getBytes());
                    osgff.write(t.printGff().getBytes());
                    
                    // valif.gff is all gencode + the validated novels isoforms
                    if(t.getIs_known() || (t.getIs_novel() && t.getIs_valid())){
                        osgffvalid.write(t.printGff().getBytes());
                        osrefflatvalid.write(t.printRefflat().getBytes());
                    }
                }
            }
            
            ostxt.close();
            osrefflat.close();
            osrefflatvalid.close();
            osgff.close();
            osgffvalid.close();
        } catch (Exception e) { e.printStackTrace(); try { ostxt.close(); osgff.close(); osgffvalid.close(); osrefflat.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { ostxt.close(); osgff.close(); osgffvalid.close(); osrefflat.close();} catch (Exception e3) { System.err.println("can not close stream");  } }
    }
    
    public List<TranscriptRecord> collapse(List<LongreadRecord> lrrList, String geneId)
    {
        List<TranscriptRecord> novel = new ArrayList<TranscriptRecord>();
        
        // run tought all undef LongReadRecord for this gene and create collpased TranscriptRecord 
        for(int i=0; i<lrrList.size(); i++) {
            LongreadRecord lrr = (LongreadRecord)lrrList.get(i);
            List<int[]> lrrJunctions = junctionsFromExons(lrr.getExons());
            
            boolean hasBeenSeenBefore = false;
            for(int j=0; j<novel.size(); j++) {
                TranscriptRecord t = (TranscriptRecord)novel.get(j);
                
                if(this.isExactSameStructure(lrrJunctions, this.junctionsFromExons(t.getExons()))){
                    if(!hasBeenSeenBefore)
                        t.add(lrr);
                    hasBeenSeenBefore = true;
                }
            }
            
            // if not in the pool of collapsed TranscriptRecord, create a new one
            if(!hasBeenSeenBefore && lrrJunctions.size()>0){
                TranscriptRecord novelT = new TranscriptRecord(geneId, "Novel." + NOVELINDEX++);
                novelT.setIs_novel(true);
                novelT.setIs_known(false);
                novelT.setExons(lrr.getExons());
                novelT.add(lrr);
                novel.add(novelT);
            }
        }
        //log.info(new Object[]{"undef collapse " + geneId + "\t" + lrrList.size() + " --> " + novel.size()});
        return novel;
    }
    
    public boolean isExactSameStructure(List<int[]> juncLrr, List<int[]> juncTr)
    {   
        boolean bool = true;
        
        // geneId  nbUmis  nbIsoformSet    nbIsoformNotSet
        // Clta    2031     1667           364
        // ENST         --> 1667
        // ENST + Novel --> 2018 only not 2031 because we do not map monexonic here !

	if(juncTr.size() > 0 && juncTr.size() == juncLrr.size()) {
            for (int i = 0; i < juncLrr.size(); i++) {
                if (!isIn((int[]) juncTr.get(i), juncLrr, this.DELTA))
                    bool = false;
            }
        }
        else
            bool = false;
        
        return bool;
    }
        
    public boolean isAllInclude(List<int[]> j1, List<int[]> j2)
    {
        // check if all junctions from j1 are in j2
        boolean bool = true;
        for (int[] a : j1) {
            if(!isIn(a,j2, this.DELTA))
                bool = false;
        }
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
}
