package org.ipmc.sicelore.utils;

import java.util.*;
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;

public class LongreadParser implements LongreadModelParser
{
	private final Log log; 
	private htsjdk.samtools.util.ProgressLogger pl;
	
	HashMap<String, Longread> mapLongreads;
	
    public LongreadParser(File bam)
    {
    	Longread longread = null;
    	log = Log.getInstance(LongreadParser.class);
		log.info(new Object[] { "\tstart..." });
    	pl = new htsjdk.samtools.util.ProgressLogger(log, 200000, "\tProcessed\t", "Records");
    	
    	this.mapLongreads = new HashMap<String, Longread>();
    	htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(bam);
    	
    	int total_records = 0;
    	int valid_records = 0;
    	
        try{
			for(SAMRecord r : inputSam){
				pl.record(r);
				LongreadRecord lrr = parseSAMRecord(r);
				total_records++;
				
				if(lrr != null){
					valid_records++;
						
					if((longread = (Longread)this.mapLongreads.get(lrr.getName())) != null)
						longread.addRecord(lrr);
					else{
						this.mapLongreads.put(lrr.getName(), new Longread(lrr.getName()));
						((Longread)this.mapLongreads.get(lrr.getName())).addRecord(lrr);
					}
				}
			}
        }catch(Exception e){ e.printStackTrace(); }
        
        log.info(new Object[] { "\tend..." });
        log.info(new Object[] { "\tTotal SAMrecords\t[" + total_records + "]"});
        log.info(new Object[] { "\tTotal SAMrecords valid\t[" + valid_records + "]"});
        log.info(new Object[] { "\tTotal Longreads\t\t[" + mapLongreads.size() + "]"});
   }
    
   public LongreadRecord parseSAMRecord(SAMRecord r) throws LongreadParseException
   {
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
   
   public HashMap<String, Longread> getMapLongreads() { return this.mapLongreads; }
}
