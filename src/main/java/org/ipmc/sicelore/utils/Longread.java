package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.util.*;
import gnu.trove.THashSet;

public class Longread implements Comparable<Longread> {
 
    private List<LongreadRecord> longreadrecords;
    private String name;
    //private String geneId;
    private THashSet<String> geneIds;
    private String barcode;
    private String umi;
    private int rn;
    //private boolean is_associated = false;

    public Longread(String name) {
        this.name = name;
        this.longreadrecords = new ArrayList<LongreadRecord>();
        this.geneIds = new THashSet<String>();
    }
    
    public int compareTo(Longread lr){
        Float obj1 = new Float(lr.getBestRecord().getDe());
        Float obj2 = new Float(this.getBestRecord().getDe());
        int retval = obj2.compareTo(obj1);
        return retval;
    }

    public String toString() {
        String str = "lr[" + name + "|" + barcode + "|" + umi + "|" + longreadrecords.size() + "|" + this.geneIds + "]";
        return str;
    }

    public void addRecord(LongreadRecord lrr)
    {
        if(lrr.getGeneId() != null){
            String[] allgenes = lrr.getGeneId().split(",");
            for(int i=0;i<allgenes.length;i++)
                this.geneIds.add(allgenes[i]);
        }
        if(lrr.getBarcode() != null)
            this.barcode = lrr.getBarcode();
        if(lrr.getUmi() != null)
            this.umi = lrr.getUmi();
        
        this.rn = lrr.getRn();
        this.longreadrecords.add(lrr);
    }
    
    public LongreadRecord getBestRecord()
    {
        Collections.sort(this.longreadrecords);
        return this.longreadrecords.get(0);
    }

    public List<LongreadRecord> getLongreadrecords() { return longreadrecords; }

    public String getName() { return name; }
    public void setName(String name) { this.name=name; }

    public String getBarcode() { return barcode; }
    public void setBarcode(String barcode) { this.barcode=barcode; }

    public String getUmi() {return umi; }
    public void setUmi(String umi) { this.umi=umi; }
    
    public int getRn() {return rn; }
    public void setRn(int rn) { this.rn=rn; }
    
    public THashSet<String> getGeneIds() { return geneIds; }
}
