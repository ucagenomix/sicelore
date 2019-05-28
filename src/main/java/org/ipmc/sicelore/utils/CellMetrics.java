package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class CellMetrics
{
    private int isoform_known_count=0;
    private int isoform_undef_count=0;
    private int nb_reads=0;
    private int nb_umis=0;
    
    public CellMetrics(){ }
    
    public void addCount(String geneId, String transcriptId, int nb_reads){
        
        this.nb_umis++;
        this.nb_reads+=nb_reads;
        
        if("undef".equals(transcriptId))
            this.isoform_undef_count++;
        else
            this.isoform_known_count++;
    }
    public int getIsoform_known_count() {
        return isoform_known_count;
    }
    public int getIsoform_undef_count() {
        return isoform_undef_count;
    }
    public int getNb_reads() {
        return nb_reads;
    }
    public int getNb_umis() {
        return nb_umis;
    }

}
