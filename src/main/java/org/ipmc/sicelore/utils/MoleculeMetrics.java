package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class MoleculeMetrics
{
    private int xReads = 0;
    private int xCleanReads = 0;
    private int xConsensusReads = 0;
    private double pctId = 0.0;
    
    public MoleculeMetrics(int xReads, int xCleanReads, int xConsensusReads, double pctId)
    {
        this.xReads = xReads;
        this.xCleanReads = xCleanReads;
        this.xConsensusReads = xConsensusReads;
        this.pctId = pctId;
    }
    
    public int getXReads() {
        return xReads;
    }
    public int getXCleanReads() {
        return xCleanReads;
    }
    public int getXConsensusReads() {
        return xConsensusReads;
    }
    public double getPctId() {
        return pctId;
    }


}
