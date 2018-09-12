package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class GeneMetrics
{
    private int isoform_known_count=0;
    private int isoform_undef_count=0;
    
    public GeneMetrics(){ }
    
    public void addCount(String geneId, String transcriptId) {
        if("undef".equals(transcriptId)){
            this.isoform_undef_count++;
        }
        else{
            this.isoform_known_count++;
        }
    }
    public int getIsoform_known_count() {
        return isoform_known_count;
    }
    public int getIsoform_undef_count() {
        return isoform_undef_count;
    }

}
