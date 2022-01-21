package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class ReadInfo
{
    private boolean is_mapped=false;
    private boolean is_incell=false;
    private boolean is_ingene=false;
    private boolean is_umi=false;
    
    public ReadInfo(){ }
    
    public boolean getIs_mapped() { return is_mapped; }
    public boolean getIs_incell() { return is_incell; }
    public boolean getIs_ingene() { return is_ingene; }
    public boolean getIs_umi() { return is_umi; }

    public void setIs_mapped(boolean is_mapped) { this.is_mapped = is_mapped; }
    public void setIs_incell(boolean is_incell) { this.is_incell = is_incell; }
    public void setIs_ingene(boolean is_ingene) { this.is_ingene = is_ingene; }
    public void setIs_umi(boolean is_umi) { this.is_umi = is_umi; }
}
