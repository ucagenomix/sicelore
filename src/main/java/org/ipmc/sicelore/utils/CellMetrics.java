/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ipmc.sicelore.utils;

/**
 *
 * @author kevin
 */
public class CellMetrics
{
    private int isoform_known_count=0;
    private int isoform_undef_count=0;
    
    public CellMetrics(){ }
    
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
