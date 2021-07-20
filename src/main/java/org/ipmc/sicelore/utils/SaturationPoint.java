package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class SaturationPoint{
    
    private int order;
    private int read;
    private double proba;
    private int count;
    private int totalRead;

    public SaturationPoint(int order, int read, double proba) {
        this.order = order;
        this.read = read;
        this.proba = proba;
        this.count = 0;
        this.totalRead = 0;
    }

    public int getOrder() { return order; }
    public void setOrder(int order) { this.order=order; }

    public int getRead() { return read; }
    public void setRead(int read) { this.read=read; }

    public double getProba() { return proba; }
    public void setProba(double proba) { this.proba=proba; }

    public int getcount() { return count; }
    public void setCount(int count) { this.count=count; }
    
    public void addCount(){ count++; }
    public void addCountRead(){ totalRead++; }
    
    @Override
    public String toString()
    {
        double saturation = 100.0 * (1.0-(new Double(count).doubleValue()/new Double(totalRead).doubleValue()));
        return order + "," + read + "," + proba + "," + count + "," + totalRead + "," + saturation;
    }
}
