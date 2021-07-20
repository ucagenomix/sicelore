package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 *  
 */

public class SubConsensus 
{
    private String name;
    private String sequence;
    
    public SubConsensus() {}
    
    public SubConsensus(String name, String sequence) {
        this.name = name;
        this.sequence = sequence;
    }
    
    public String getName() { return this.name; }
    public String getSequence() { return this.sequence; }
}
