package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class Junction{
    
    private int start;
    private int end;

    public Junction(int start, int end) {
        this.start = start;
        this.end = end;
    }

    public String toString() {
        return "[" + this.start + "-" + this.end + "]";
    }
    public int[] getIntArray() {
        return new int[] {this.start,this.end};
    }

    public int getStart() { return start; }
    public void setStart(int start) { this.start=start; }

    public int getEnd() { return end; }
    public void setEnd(int end) { this.end=end; }

        
    @Override
    public int hashCode(){
        int result=17;
        result=31*result+start+end;
        return result;
    }
    
    @Override
    public boolean equals(final Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        
        final Junction other = (Junction)obj;
        return other.getStart() == this.getStart() && other.getEnd() == this.getEnd();
    }
}
