package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */

public class FastqRecord {

    String name;
    String seq;
    String qual;
    
    public FastqRecord(String name,String seq,String qual)
    {
        this.name=name;
        this.seq=seq;
        this.qual=qual;
    }

    public String getName() { return name; }
    public String getSeq() { return seq; }
    public String getQual() { return qual; }
    
    public void setName(String name) { this.name=name; }
    public void setSeq(String seq) { this.seq=seq; }
    public void setQual(String qual) { this.qual=qual; }
}
