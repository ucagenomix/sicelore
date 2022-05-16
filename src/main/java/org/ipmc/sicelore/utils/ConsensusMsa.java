package org.ipmc.sicelore.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashSet;
import java.util.regex.Pattern;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class ConsensusMsa
{
    private String cons;
    private String qv;
    private HashSet<char[]> msa;
            
    protected static int MAXPS;
    
    public ConsensusMsa() {}
    
    public ConsensusMsa(String file)
    {
        msa = new HashSet<char[]>();
        
        try{
            BufferedReader input = new BufferedReader(new FileReader(file));
            String line = input.readLine();
            while(line != null && !"".equals(line)){
                if(Pattern.matches("^>Consensus.*", line))
                    this.cons = input.readLine();
                else if(Pattern.matches("^Consensus.*", line))
                    this.cons = input.readLine();
                else
                    msa.add(input.readLine().toCharArray());
                
                line = input.readLine();
            }                
            input.close();
        }catch(Exception e){ e.printStackTrace(); }
        
        process();
    }
    
    public void setStaticParams(int MAXPS)
    {
        this.MAXPS = MAXPS;
    }

    public void process()
    {
        String seq = "";
        char[] consArray = cons.toCharArray();
        int size = cons.length();
        double[] qvs = new double[size];
        for(int i=0; i<consArray.length; i++){
            double same = 0;
            for(char[] s: msa){
               if(s[i] == consArray[i])
                    same++;
            }
            qvs[i] = same/msa.size();
            //System.out.println(same + "/" + msa.size() + "=" + qvs[i]);
        }
        
        this.cons = this.cons.replaceAll("-", "");
        int consSize = this.cons.length();
        int[] final_qvs = new int[consSize];
        int index=0;
        for(int i=0; i<consArray.length; i++){
            if(consArray[i] != '-'){
                seq += consArray[i];
                if(qvs[i] == 1.0)
                    final_qvs[index] = 33 + MAXPS;
                else
                    final_qvs[index] = 33 + (int)Math.round(-10*Math.log10(1.0-qvs[i]));
                
                //System.out.println(i + " = " + qvs[i] + "\t" + final_qvs[index]);
                
                index++;
            }
        }
        this.cons = seq;
        this.qv = "";
        for(int i: final_qvs){
            String s = Character.toString((char)i);
            this.qv += s;
            //System.out.println(i + " => " + s);
        }
    }
    
    public String getCons() { return cons; }
    public String getQv() { return qv; }
}
