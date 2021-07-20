package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class MoleculeMetrics
{
    private int xReads = 0;
    private double pctId = 0.0;
    
    public MoleculeMetrics(int xReads, double pctId, int zero)
    {
        this.xReads = xReads;
        this.pctId = pctId;
    }
    
    public int getXReads() {
        return xReads;
    }
    public double getPctId() {
        return pctId;
    }
}
