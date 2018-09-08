package org.ipmc.sicelore.utils;

import java.util.*;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.CloseableIterator;

public class Longread implements Comparable<Longread> {

    private List<LongreadRecord> longreadrecords;
    private String name;
    private String geneId;
    private String barcode;
    private String umi;
    private boolean is_associated = false;

    public Longread(String name) {
        this.name = name;
        this.longreadrecords = new ArrayList<LongreadRecord>();
    }
    
    public int compareTo(Longread lr)
    {
        Float obj1 = new Float(lr.getAssociatedRecord().getDv());
        Float obj2 = new Float(this.getAssociatedRecord().getDv());
        int retval = obj2.compareTo(obj1);
        return retval;

        //return ((Longread)lr).getAssociatedRecord().getExonBases() - this.getAssociatedRecord().getExonBases();
    }

    public String toString() {
        String str = "[" + name + "|" + geneId + "|" + barcode + "|" + umi + "|" + longreadrecords.size() + "]";
        return str;
    }

    public void addRecord(LongreadRecord lrr)
    {
        this.longreadrecords.add(lrr);

        // is associated record
        if (lrr.getIs_associated()) {
            this.geneId = lrr.getGeneId();
            this.barcode = lrr.getBarcode();
            this.umi = lrr.getUmi();
            this.is_associated = true;
        }
    }

    public LongreadRecord getAssociatedRecord() {
        LongreadRecord rec = null;
        for (LongreadRecord lrr : this.longreadrecords) {
            if (lrr.getIs_associated()) {
                rec = lrr;
            }
        }
        //System.out.println(name+"\t"+rec);
        return rec;
    }

    public List<LongreadRecord> getLongreadrecords() {
        return longreadrecords;
    }

    public String getName() {
        return name;
    }

    public String getGeneId() {
        return geneId;
    }

    public String getBarcode() {
        return barcode;
    }

    public String getUmi() {
        return umi;
    }

    public boolean getIs_associated() {
        return is_associated;
    }
}
