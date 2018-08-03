package org.ipmc.common.gtf;

import java.util.*;
import java.io.*;
import java.util.concurrent.Callable;
import org.biojava.nbio.core.sequence.io.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.sequence.DNASequence;
import org.biojava.nbio.core.sequence.io.FastaReaderHelper;
import org.biojava.nbio.core.util.ConcurrencyTools;
import org.biojava.nbio.core.sequence.compound.DNACompoundSet; 
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.biojava.nbio.core.alignment.SimpleAlignedSequence;
import org.biojava.nbio.core.alignment.template.AlignedSequence;

public class GTFjunction
{
   private String chromosome;
   private int start;
   private int end;
   private String strand;
   
   public int count;
   public String label;
   
   public GTFjunction(){}
   
   public GTFjunction(int start, int end)
   {
       this.start = start;
       this.end = end;
   }
   
   public int getStart(){ return this.start; }
   public int getEnd(){ return this.end; }
   public int getCount(){ return this.count; }
   public void setCount(int count ){ this.count=count; }
   public String getLabel(){ return this.label; }
   public void setLabel(String label){ this.label=label; }
   
   public void addCount()
   {
	   this.count++;
   }
}