package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.util.*;
import java.io.*;
import htsjdk.samtools.util.Log;

public class Matrix
{
    private final Log log;
    private HashMap<String, HashMap<String, Integer>> matrice;
    private HashMap<String, HashMap<String, Integer>> matriceGene;
    private HashMap<String, CellMetrics> cellMetrics;
    private HashMap<String, GeneMetrics> geneMetrics;

    private int nb_count = 0;
    private int total_count = 0;
    private int total_remove = 0;
    private int total_isoform_def = 0;
    private int total_isoform_undef = 0;
    
    public Matrix(HashSet<String> cells)
    {
        log = Log.getInstance(Matrix.class);
        this.matrice = new HashMap<String, HashMap<String, Integer>>();
        this.matriceGene = new HashMap<String, HashMap<String, Integer>>();
        this.geneMetrics = new HashMap<String, GeneMetrics>();
        this.cellMetrics = new HashMap<String, CellMetrics>();
        Iterator it = cells.iterator();
        while (it.hasNext())
            cellMetrics.put((String)it.next(), new CellMetrics());
    }

    public HashMap<String, HashMap<String, Integer>> getMatrice(){ return matrice; }
    public HashMap<String, HashMap<String, Integer>> getMatriceGene(){ return matriceGene; }
    public HashMap<String, CellMetrics> getCellMetrics(){ return cellMetrics; }
    public HashMap<String, GeneMetrics> getGeneMetrics(){ return geneMetrics; }
    
    public int getNb_count(){ return nb_count; }
    public int getTotal_count(){ return total_count; }
    public int getTotal_remove(){ return total_remove; }
    public int getTotal_isoform_def(){ return total_isoform_def; }
    public int getTotal_isoform_undef(){ return total_isoform_undef; } 
       
    public void addMolecule(Molecule molecule)
    {
        HashMap mapCell = null;
        HashSet setUmi = null;
        
        if(! geneMetrics.containsKey(molecule.getGeneId()))
            geneMetrics.put(molecule.getGeneId(), new GeneMetrics());
        
        //System.out.println(molecule.getBarcode()+","+molecule.getGeneId()+","+molecule.getTranscriptId());
        
        ((CellMetrics)cellMetrics.get(molecule.getBarcode())).addCount(molecule.getGeneId(), molecule.getTranscriptId());
        ((GeneMetrics)geneMetrics.get(molecule.getGeneId())).addCount(molecule.getGeneId(), molecule.getTranscriptId());
        
        if ("undef".equals(molecule.getTranscriptId()))
            total_isoform_undef++;
        else
            total_isoform_def++;
        
        // we already have this gene/transcript isokey
        String isokey = molecule.getGeneId() + "\t" + molecule.getTranscriptId();
        if ((mapCell = (HashMap) matrice.get(isokey)) != null) {
            // we already have this cell
            if ((setUmi = (HashSet) mapCell.get(molecule.getBarcode())) != null) {
                setUmi.add(molecule.getUmi());
            }
            else {
                setUmi = new HashSet();
                setUmi.add(molecule.getUmi());
                mapCell.put(molecule.getBarcode(), setUmi);
            }
        }
        else {
            matrice.put(isokey, new HashMap());
            setUmi = new HashSet();
            setUmi.add(molecule.getUmi());
            ((HashMap) matrice.get(isokey)).put(molecule.getBarcode(), setUmi);
        }
        
        // now work on matrixGene
        isokey = molecule.getGeneId();
        if ((mapCell = (HashMap) matriceGene.get(isokey)) != null) {
            // we already have this cell
            if ((setUmi = (HashSet) mapCell.get(molecule.getBarcode())) != null) {
                setUmi.add(molecule.getUmi());
            }
            else {
                setUmi = new HashSet();
                setUmi.add(molecule.getUmi());
                mapCell.put(molecule.getBarcode(), setUmi);
            }
        }
        else {
            matriceGene.put(isokey, new HashMap());
            setUmi = new HashSet();
            setUmi.add(molecule.getUmi());
            ((HashMap) matriceGene.get(isokey)).put(molecule.getBarcode(), setUmi);
        }
    }

    public void writeIsoformMatrix(java.io.File isomatrix,java.io.File isometrics)
    {
        DataOutputStream os = null;
        DataOutputStream os2 = null;
        HashSet setUmi = null;
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(isomatrix));
            os2 = new DataOutputStream(new java.io.FileOutputStream(isometrics));
            
            os.writeBytes("geneId\ttranscriptId");
            for(String key : cellMetrics.keySet())
                os.writeBytes("\t" + key);
            os.writeBytes("\n");

            for(String isokey : matrice.keySet()){
                os.writeBytes(isokey);
                os2.writeBytes(isokey);
                int total = 0;
                for(String cell_barcode : cellMetrics.keySet()){
                    if ((setUmi = (HashSet) ((HashMap) matrice.get(isokey)).get(cell_barcode)) != null) {
                        os.writeBytes("\t" + setUmi.size());
                        this.total_count += setUmi.size();
                        total += setUmi.size();
                    }
                    else
                        os.writeBytes("\t0");
                }
                os.writeBytes("\n");
                os2.writeBytes("\t" + total + "\n");
            }
            os.close();
            os2.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close();  os2.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close();  os2.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
    
    public void writeGeneMatrix(java.io.File paramFile)
    {
        DataOutputStream os = null;
        HashSet setUmi = null;
        
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(paramFile));
            os.writeBytes("geneId");
            for(String cell_barcode : cellMetrics.keySet())
                os.writeBytes("\t" + cell_barcode);
            os.writeBytes("\n");

            for(String isokey : matriceGene.keySet()){
                os.writeBytes(isokey);
                for(String cell_barcode : cellMetrics.keySet()){
                    if ((setUmi = (HashSet) ((HashMap) matriceGene.get(isokey)).get(cell_barcode)) != null)
                        os.writeBytes("\t" + setUmi.size());
                    else
                        os.writeBytes("\t0");
                }
                os.writeBytes("\n");
            }
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }

    public void writeCellMetrics(java.io.File paramFile)
    {
        DataOutputStream os = null;
        
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(paramFile));
            os.writeBytes("barcode\ttotal\tisoform_def\tisoform_undef\n");
            
            for(String cell_barcode : cellMetrics.keySet()){
                CellMetrics cm = cellMetrics.get(cell_barcode);
                if(cm != null){
                    int total = cm.getIsoform_known_count() + cm.getIsoform_undef_count();
                    os.writeBytes(cell_barcode+"\t"+total+"\t"+cm.getIsoform_known_count()+"\t" +cm.getIsoform_undef_count()+"\n");
                }
                else{
                    os.writeBytes(cell_barcode+"\t0\t0\t0\n");
                }
            }
            
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
    
    public void writeGeneMetrics(java.io.File paramFile)
    {
        DataOutputStream os = null;
        
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(paramFile));
            os.writeBytes("geneId\ttotal\tisoform_def\tisoform_undef\n");
            for(String geneId : geneMetrics.keySet())
            {
                GeneMetrics gm = this.geneMetrics.get(geneId);
                if(gm != null){
                    int total = gm.getIsoform_known_count() + gm.getIsoform_undef_count();
                    os.writeBytes(geneId+"\t"+total+"\t"+gm.getIsoform_known_count()+"\t" +gm.getIsoform_undef_count()+"\n");
                }
                else{
                    os.writeBytes(geneId+"\t0\t0\t0\n");
                }
            }
            
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
}
