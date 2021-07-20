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
import gnu.trove.THashMap;

public class LongreadParser implements LongreadModelParser {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;
    private int unvalid_records = 0;
    private int chimeria_records = 0;
    private int reversed_records = 0;
    private int not_primary_records = 0;
    private int mapqv0_records = 0;
    private int null_records = 0;
    private int gene_unset = 0;
    private int umi_unset = 0;
    
    private boolean keep_mapqv0 = false;
    private boolean load_sequence = false;
    private boolean is_gene_mandatory = true;
    private boolean is_umi_mandatory = true;
    
    THashMap<String, Longread> mapLongreads;
    
    public LongreadParser()
    {
        log = Log.getInstance(LongreadParser.class);
        this.mapLongreads = new THashMap<String, Longread>();
    }
    
    public LongreadParser(File bam, boolean keep_mapqv0, boolean load_sequence, boolean is_gene_mandatory, boolean is_umi_mandatory)
    {
        this.keep_mapqv0 = keep_mapqv0;
        this.load_sequence = load_sequence;
        this.is_gene_mandatory = is_gene_mandatory;
        this.is_umi_mandatory = is_umi_mandatory;
        
        Longread longread = null;
        log = Log.getInstance(LongreadParser.class);
        log.info(new Object[]{"\tstart..."});
        pl = new htsjdk.samtools.util.ProgressLogger(log, 1000000, "\tProcessed\t", "Records");

        this.mapLongreads = new THashMap<String, Longread>();
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(bam);

        int total_records = 0;
        int valid_records = 0;
        HashSet<String> multiRec = new HashSet<String>();
        
        try {
            for(SAMRecord r : inputSam) {
                pl.record(r);
                LongreadRecord lrr = parseSAMRecord(r);
                total_records++;

                // IG/BC/U8
                if (lrr != null) {
                    valid_records++;

                    if ((longread = (Longread) this.mapLongreads.get(lrr.getName())) != null){
                        longread.addRecord(lrr);
                        multiRec.add(longread.getName());
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
        log.info(new Object[]{"\tTotal SAMrecords\t" + total_records});
        log.info(new Object[]{"\tSAMrecords valid\t" + valid_records});
        log.info(new Object[]{"\tSAMrecords unvalid\t" + unvalid_records});
        log.info(new Object[]{"\tSAMrecords mapqv=0\t" + mapqv0_records});
        log.info(new Object[]{"\tSAMrecords no gene\t" + gene_unset});
        log.info(new Object[]{"\tSAMrecords no UMI\t" + umi_unset});
        log.info(new Object[]{"\tSAMrecords chimeria\t" + chimeria_records});
        log.info(new Object[]{"\tSAMrecords reversed\t" + reversed_records});
        //log.info(new Object[]{"\tSAMrecords null\t[" + null_records + "]"});
        //log.info(new Object[]{"\tSAMrecords not primary\t[" + not_primary_records + "]"});
        log.info(new Object[]{"\tTotal reads\t\t" + mapLongreads.size()});
        log.info(new Object[]{"\tTotal reads multiSAM\t" + multiRec.size()});
    }

    public LongreadRecord parseSAMRecord(SAMRecord r) throws LongreadParseException
    {
        LongreadRecord record = LongreadRecord.fromSAMRecord(r, this.load_sequence);
        
        if(record == null) { unvalid_records++; null_records++; return null; }
        if(record.getIsChimeria()) { unvalid_records++; chimeria_records++; return null; }
        if(record.getIsReversed()) { unvalid_records++; reversed_records++; return null; }
        if(this.is_gene_mandatory && (record.getGeneId() == null || "undef".equals(record.getGeneId()))) { unvalid_records++; gene_unset++; return null; }
        if(this.is_umi_mandatory && record.getUmi() == null) { unvalid_records++; umi_unset++; return null; }
        if(!this.keep_mapqv0 && record.getMapqv() == 0) { unvalid_records++; mapqv0_records++; return null; }
        
        return record;
    }

    public THashMap<String, Longread> getMapLongreads() {
        return this.mapLongreads;
    }
}
