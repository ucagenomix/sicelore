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

public class LongreadParser implements LongreadModelParser {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;
    private int unvalid_records = 0;
    private int chimeria_records = 0;
    private int reversed_records = 0;
    private boolean load_sequence = false;
    
    HashMap<String, Longread> mapLongreads;

    public LongreadParser(File bam, boolean load_sequence)
    {
        this.load_sequence = load_sequence;
        
        Longread longread = null;
        log = Log.getInstance(LongreadParser.class);
        log.info(new Object[]{"\tstart..."});
        pl = new htsjdk.samtools.util.ProgressLogger(log, 500000, "\tProcessed\t", "Records");

        this.mapLongreads = new HashMap<String, Longread>();
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(bam);

        int total_records = 0;
        int valid_records = 0;
        //HashSet<String> multiRec = new HashSet<String>();
        
        try {
            for(SAMRecord r : inputSam) {
                pl.record(r);
                LongreadRecord lrr = parseSAMRecord(r);
                total_records++;
                
                //if(total_records%1000000 == 0){
                //    System.gc();System.gc();System.gc();System.gc();System.gc();System.gc();
                //    log.info(new Object[]{"\tSystem.gc() done..."});
                //}
                
                // IG/BC/U8
                if (lrr != null) {
                    valid_records++;

                    if ((longread = (Longread) this.mapLongreads.get(lrr.getName())) != null) {
                        longread.addRecord(lrr);
                        //multiRec.add(longread.getName());
                    } 
                    else {
                        this.mapLongreads.put(lrr.getName(), new Longread(lrr.getName()));
                        ((Longread)this.mapLongreads.get(lrr.getName())).addRecord(lrr);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        log.info(new Object[]{"\tend..."});
        log.info(new Object[]{"\tTotal SAMrecords\t[" + total_records + "]"});
        log.info(new Object[]{"\tSAMrecords valid\t[" + valid_records + "]"});
        log.info(new Object[]{"\tSAMrecords unvalid\t[" + unvalid_records + "]"});
        log.info(new Object[]{"\tSAMrecords chimeria\t[" + chimeria_records + "]"});
        log.info(new Object[]{"\tSAMrecords reversed\t[" + reversed_records + "]"});
        log.info(new Object[]{"\tTotal reads\t\t[" + mapLongreads.size() + "]"});
        //log.info(new Object[]{"\tTotal reads multiSAM\t[" + multiRec.size() + "]"});
    }

    public LongreadRecord parseSAMRecord(SAMRecord r) throws LongreadParseException
    {
        LongreadRecord record = LongreadRecord.fromSAMRecord(r, this.load_sequence);
        if(record == null) { unvalid_records++; return null; }
        if(record.getIsChimeria()) { chimeria_records++; return null; }
        if(record.getIsReversed()) { reversed_records++; return null; }
        return record;
    }

    public HashMap<String, Longread> getMapLongreads() {
        return this.mapLongreads;
    }
}
