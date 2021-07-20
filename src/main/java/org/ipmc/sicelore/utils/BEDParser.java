package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */ 
import java.io.*;
import java.util.*;
import htsjdk.samtools.util.Log;
import htsjdk.tribble.annotation.Strand;

import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;
import htsjdk.tribble.bed.BEDCodec.StartOffset;
import htsjdk.tribble.readers.AsciiLineReader;
import htsjdk.tribble.readers.LineIteratorImpl;

public class BEDParser {

    private final Log log;
    List<BEDFeature> allFeatures;
    int entries;
    
    HashMap<String, List<BEDFeature>> chrToBF;
    
    public BEDParser(File bedFile)
    {
        log = Log.getInstance(UCSCRefFlatParser.class);
        
        
        chrToBF = new HashMap<String, List<BEDFeature>>();
        allFeatures = new ArrayList<BEDFeature>();
        
        try {
            InputStream is = new java.io.FileInputStream(bedFile);
            AsciiLineReader lr = new AsciiLineReader(is);
            LineIteratorImpl li = new LineIteratorImpl(lr);
            BEDCodec bc = new BEDCodec(StartOffset.ZERO);
            while(!bc.isDone(li))
            {
                BEDFeature bf = bc.decode(li);
                if(bf != null)
                {
                    allFeatures.add(bf);
                    if(chrToBF.containsKey(bf.getChr())){
                        chrToBF.get(bf.getChr()).add(bf);
                    }
                    else{
                        List<BEDFeature> a = new ArrayList<BEDFeature>();
                        a.add(bf);
                        chrToBF.put(bf.getChr(), a);
                    }
                }
            }
            this.entries = allFeatures.size();
        } catch (Exception e) { e.printStackTrace(); }
        
        log.info(new Object[]{"BEDParser\t" + bedFile + "\t[references="+this.chrToBF.size()+",entries="+this.entries+"]"});
    }
    public HashMap<String, List<BEDFeature>> getChrToBF(){ return chrToBF; }
    public int getEntries(){ return entries; }
    
    // distance relative to start of CAGE peak
    // negative value mean CAGE peak is upstream transcript start
    // 0 mean transcript start at the start of CAGE signal
    // positive value mean transcript cover the CAGE peak start
    public int getDistanceCage(String chromosome, Strand strand, int pos)
    {
        BEDFeature evidence = null;
        int min=Integer.MAX_VALUE;
        int minabs=Integer.MAX_VALUE;
        if(chrToBF.containsKey(chromosome)){
            List<BEDFeature> a = chrToBF.get(chromosome);
            for (BEDFeature f : a) {
                if(strand.equals(f.getStrand())){
                    int pp = (strand.equals(Strand.POSITIVE))?f.getStart():f.getEnd();
                    //int pp = f.getStart();
                    if(Math.abs(pos-pp) < minabs){
                        min = pos-pp;
                        minabs = Math.abs(pos-pp);
                        evidence = f;
                    }
                }
            }
        }
        if(strand.equals(Strand.POSITIVE)){ min = -min; }
        //log.info(new Object[]{evidence.getContig() + ":" + evidence.getStart() + "-" + evidence.getEnd() + "," + evidence.getStrand()});
        
        return min;
    }
    
    // distance relative to end of polyA signal
    // negative value mean transcript cover the polyA
    // 0 mean transcript end on the end of the polyA
    // positive value mean polyA is further than end of transcript
    public int getDistancePolyA(String chromosome, Strand strand, int pos)
    {
        BEDFeature evidence = null;
        int min=Integer.MAX_VALUE;
        int minabs=Integer.MAX_VALUE;
        if(chrToBF.containsKey(chromosome)){
            List<BEDFeature> a = chrToBF.get(chromosome);
            for (BEDFeature f : a) {
                if(strand.equals(f.getStrand())){
                    int pp = (strand.equals(Strand.POSITIVE))?f.getStart():f.getEnd();
                    if(Math.abs(pos-pp) < minabs){
                        min = pos-pp;
                        minabs = Math.abs(pos-pp);
                        evidence = f;
                    }
                }
            }
        }
        if(strand.equals(Strand.POSITIVE)){ min = -min; }
        //log.info(new Object[]{evidence.getContig() + ":" + evidence.getStart() + "-" + evidence.getEnd() + "," + evidence.getStrand()});
        
        return min;
    }
    
    
    public boolean isWithin(String chromosome, Strand strand, int pos)
    {
        boolean within = false;
        if(chrToBF.containsKey(chromosome)){
            List<BEDFeature> a = chrToBF.get(chromosome);
            for (BEDFeature f : a) {
                if(strand.equals(f.getStrand())){
                    if(pos < f.getEnd() && pos > f.getStart())
                        within=true;
                }
            }
        }
        return within;
    }
}






/*
public class UCSCRefFlatParser implements GeneModelParser {

    private final Log log;
    private HashMap<String, List<TranscriptRecord>> mapGenesTranscripts;

    public int DELTA = 2;
    public int MINEVIDENCE = 5;
    public int NOVELINDEX=1;
    
    HashSet<String> knownJunctions;
    HashSet<String> knownSpliceSites;
    HashSet<String> novelSpliceSites;
               
    public UCSCRefFlatParser(File refFlat)
    {
        log = Log.getInstance(UCSCRefFlatParser.class);
        log.info(new Object[]{"UCSCRefFlatParser\tstart..."});

        this.mapGenesTranscripts = new HashMap<String, List<TranscriptRecord>>();
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
    
    public UCSCRefFlatParser(int DELTA, int MINEVIDENCE)
    {
        log = Log.getInstance(UCSCRefFlatParser.class);
        this.mapGenesTranscripts = new HashMap<String, List<TranscriptRecord>>();
        
        this.DELTA =DELTA;
        this.MINEVIDENCE=MINEVIDENCE;
        this.NOVELINDEX=1;
        knownJunctions = new HashSet<String>();
        knownSpliceSites = new HashSet<String>();
        novelSpliceSites = new HashSet<String>();
    }

    public TranscriptRecord parseLine(String line) throws GTFParseException {
        String[] fields = line.split("\t");
        TranscriptRecord record = TranscriptRecord.fromRefFlat(fields);

        if (record.getChrom().contains("_")) {
            return null;
        }
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

    public HashMap<String, List<TranscriptRecord>> getMapGenesTranscripts() {
        return mapGenesTranscripts;
    }
    
    // parcours de la map de genes and collapse undef categorie
    public void collapser()
    {
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
    public void setter()
    {    
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            for(int i=0; i<tAll.size(); i++)
                tAll.get(i).initialize();
        }
    }
    
    // filtration of 3p-part degradated isoforms
    public void filter(UCSCRefFlatParser model)
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
                    if(! is3pPart(t, keep, model.select(new String[]{t.getGeneId()})))
                        keep.add(t);
                }
            }
            mapGenesTranscripts.put(geneId, keep);
        }
    }
    
    // classification
    public void classifier(UCSCRefFlatParser model)
    {
        for(String geneId : mapGenesTranscripts.keySet()) {
            List<TranscriptRecord> tAll = mapGenesTranscripts.get(geneId);
            for(int i=0; i<tAll.size(); i++){
                TranscriptRecord t = tAll.get(i);
                if(t.getIs_novel())
                    this.noveltyDetector(t, model.select(new String[]{t.getGeneId()}));
            }
        }
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
            if(!isIn(junctions.get(i), modelJunctions)){
                
                // cette junctions n'est pas connue
                // est-elle formées sur la base des sites de splices connus ?
                // meaning start and stop are in the modelJunctions start/stop
                if(modelSplice.contains(junctions.get(i)[0]) && modelSplice.contains(junctions.get(i)[1])){
                    if("undef".equals(t.getCategorie())){
                        t.setCategorie("novel_in_catalog");
                        t.setSubcategorie("combination_of_known_splicesites");
                    }
                    t.setWhynovel(t.getWhynovel() + new String("kss=" + t.getChrom() + ":" + junctions.get(i)[0] + "-" + junctions.get(i)[1] + ","));
                    this.knownSpliceSites.add(new String(t.getChrom() + ":" + junctions.get(i)[0] + "-" + junctions.get(i)[1]));
                }
                else{
                    t.setCategorie("novel_not_in_catalog");
                    t.setSubcategorie("at_least_one_novel_splicesite");
                    t.setWhynovel(t.getWhynovel() + new String("nss=" + t.getChrom() + ":" + junctions.get(i)[0] + "-" + junctions.get(i)[1] + ","));
                    this.novelSpliceSites.add(new String(t.getChrom() + ":" + junctions.get(i)[0] + "-" + junctions.get(i)[1]));
                    //log.info(new Object[]{t.getTranscriptId() + "\t" + junctions.get(i)[0] + "-" + junctions.get(i)[1] + " is novel"});
                }
            }
            else{
                this.knownJunctions.add(new String("kj=" + t.getChrom() + ":" + junctions.get(i)[0] + "-" + junctions.get(i)[1]));
            }
        }
        
        // si toutes les junctions sont connues
        if("undef".equals(t.getCategorie())){
            t.setCategorie("novel_in_catalog");
            t.setSubcategorie("combination_of_known_junctions");
        }
    }
    
    public boolean is3pPart(TranscriptRecord t, List<TranscriptRecord> lst, List<TranscriptRecord> modelLst)
    {
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
    
    // statistiques
    public void statistics()
    {
        int sum_isoforms = 0;
        int sum_evidences = 0;
        HashMap<String, Integer> stats = new HashMap<String, Integer>();
        String[] allkeys = {"undef","undef2","full_splice_match","gencode","novel_in_catalog","novel_not_in_catalog","combination_of_known_junctions","combination_of_known_splicesites","at_least_one_novel_splicesite"};
        for(int i=0; i<allkeys.length; i++){
            stats.put(allkeys[i]+"_count", 0);
            stats.put(allkeys[i]+"_evidences", 0);
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
            }
        }
        
        log.info(new Object[]{"-STATS-----------------------------------------------------------"});
        log.info(new Object[]{String.format("|\t\t\t\t\t| nb\t| evidences\t|")});
        log.info(new Object[]{String.format("|total_genes\t\t\t\t|%d\t|\t\t|",mapGenesTranscripts.size())});
        log.info(new Object[]{String.format("|total_isoforms\t\t\t\t|%d\t| %d\t\t|",sum_isoforms,sum_evidences)});
        log.info(new Object[]{String.format("|full_splice_match\t\t\t|\t|\t\t|")});
        log.info(new Object[]{String.format("| o gencode\t\t\t\t|%d\t| %d\t\t|",stats.get("gencode_count"),stats.get("gencode_evidences"))});
        log.info(new Object[]{String.format("|novel_in_catalog\t\t\t|\t|\t\t|")});
        log.info(new Object[]{String.format("| o combination_of_known_junctions\t|%d\t| %d\t\t|",stats.get("combination_of_known_junctions_count"),stats.get("combination_of_known_junctions_evidences"))});
        log.info(new Object[]{String.format("| o combination_of_known_splicesites\t|%d\t| %d\t\t|",stats.get("combination_of_known_splicesites_count"),stats.get("combination_of_known_splicesites_evidences"))});
        log.info(new Object[]{String.format("|novel_not_in_catalog\t\t\t|\t|\t\t|")});
        log.info(new Object[]{String.format("| o at_least_one_novel_splicesite\t|%d\t| %d\t\t|",stats.get("at_least_one_novel_splicesite_count"),stats.get("at_least_one_novel_splicesite_evidences"))});
        log.info(new Object[]{"-----------------------------------------------------------------"});
    }
    
    //export files
    public void exportFiles(File TXT, File GFF, File FAS)
    {
        BufferedOutputStream ostxt = null;
        BufferedOutputStream osgff = null;
        BufferedOutputStream osfas = null;
        
        try {
            ostxt = new BufferedOutputStream(new java.io.FileOutputStream(TXT));
            osgff = new BufferedOutputStream(new java.io.FileOutputStream(GFF));
            osfas = new BufferedOutputStream(new java.io.FileOutputStream(FAS));
            
            for(String geneId : mapGenesTranscripts.keySet()) {
                // running through all ENST + NOVEL
                List<TranscriptRecord> tList = mapGenesTranscripts.get(geneId);

                for(int i=0; i<tList.size(); i++) {
                    TranscriptRecord t = tList.get(i);
                    
                    ostxt.write(t.printTxt().getBytes());
                    osgff.write(t.printGff().getBytes());
                    osfas.write(t.printFas().getBytes());
                }
            }
            
            ostxt.close();
            osgff.close();
            osfas.close();
        } catch (Exception e) { e.printStackTrace(); try { ostxt.close(); osgff.close(); osfas.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { ostxt.close(); osgff.close(); osfas.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
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
            
            // if not in the pool of collpased TranscriptRecord, create a new one
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
                if (!isIn((int[]) juncTr.get(i), juncLrr))
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
            if(!isIn(a,j2))
                bool = false;
        }
        return bool;
    }
    
    public boolean isIn(int[] junction, List<int[]> list)
    {
        boolean bool = false;
        for (int[] a : list) {
            if ((Math.abs(a[0] - junction[0]) <= DELTA) && (Math.abs(a[1] - junction[1]) <= DELTA))
                bool = true;
        }
        return bool;
    }
    
    public boolean isSameJunction(int[] j1, int[] j2)
    {
        boolean bool = false;
        if ((Math.abs(j1[0] - j2[0]) <= DELTA) && (Math.abs(j1[1] - j2[1]) <= DELTA))
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
}
*/