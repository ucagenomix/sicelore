package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class JunctionMetrics
{
    private int isoform_known_count=0;
    private int isoform_undef_count=0;
    private int nb_reads=0;
    private int nb_umis=0;
    
    public JunctionMetrics(){ }
    
    public void addCount(String geneId, int[] pos){
        
        this.nb_umis++;
        this.nb_reads+=nb_reads;
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
