package org.ipmc.sicelore.utils;

import java.io.*;
import java.util.*;
import htsjdk.samtools.util.Log;

public class UCSCRefFlatParser implements GeneModelParser
{
	private final Log log; 
    private HashMap<String, List<TranscriptRecord>> mapGenesTranscripts;
    
    public UCSCRefFlatParser(File refFlat)
    {
    	log = Log.getInstance(UCSCRefFlatParser.class);
		log.info(new Object[] { "UCSCRefFlatParser\tstart..." });
		
		this.mapGenesTranscripts = new HashMap<String, List<TranscriptRecord>>();
        int nb=0;
        
        try{
	        FileInputStream in = new FileInputStream(refFlat);
	        BufferedReader reader = new BufferedReader(new InputStreamReader(in));
	        String line = null;
	        while((line = reader.readLine()) != null){
	            TranscriptRecord f = parseLine(line);
	            if(f!=null){
	            	nb++;
	            	if(this.mapGenesTranscripts.containsKey(f.getGeneId())){
	            		((ArrayList<TranscriptRecord>)this.mapGenesTranscripts.get(f.getGeneId())).add(f);
	            	}
	            	else{
	            		List<TranscriptRecord> l = new ArrayList<TranscriptRecord>();
	            		l.add(f);
	            		this.mapGenesTranscripts.put(f.getGeneId(), l);
	            	}
	            }
	        }
        }catch(Exception e){ e.printStackTrace(); }
        
        log.info(new Object[] { "UCSCRefFlatParser\tend..." });
		log.info(new Object[] { "Number of Genes Symbols\t[" + this.mapGenesTranscripts.size() + "]"});
		log.info(new Object[] { "Number of Transcripts\t[" + nb + "]"});
    }

    public TranscriptRecord parseLine(String line) throws GTFParseException
    {
        String[] fields = line.split("\t");
        TranscriptRecord record = TranscriptRecord.fromRefFlat(fields);
        
        if(record.getChrom().contains("_")) return null;
        if(record.getExonBases() == 0) return null;

        return record;
    }
    
    public List<TranscriptRecord> select(String[] filter)
    {
    	List<TranscriptRecord> filteredList = new ArrayList<TranscriptRecord>();
        for(String f : filter){
            if(this.mapGenesTranscripts.containsKey(f)){
            	for (TranscriptRecord t : this.mapGenesTranscripts.get(f))
            		filteredList.add(t);
            }
        }
        return filteredList;
    }
    
    public HashMap<String, List<TranscriptRecord>> getMapGenesTranscripts() { return mapGenesTranscripts; }
 }