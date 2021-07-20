package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */ 
import java.util.*;
import java.io.*;
import htsjdk.samtools.util.Log;
import gnu.trove.THashMap;
import gnu.trove.THashSet;

public class Matrix
{
    public final Log log;
    public THashMap<String, THashMap<String, Integer>> matrice;
    public THashMap<String, THashMap<String, Integer>> matriceGene;
    public THashMap<String, THashMap<String, Integer>> matriceJunction;
    public THashMap<String, CellMetrics> cellMetrics;
    public THashMap<String, GeneMetrics> geneMetrics;
    public THashSet<Molecule> molecules;

    public int nb_count = 0;
    public int total_count = 0;
    public int total_remove = 0;
    public int total_isoform_def = 0;
    public int total_isoform_undef = 0;
    
    public Matrix(HashSet<String> cells)
    {
        log = Log.getInstance(Matrix.class);
        this.matrice = new THashMap<String, THashMap<String, Integer>>();
        this.matriceGene = new THashMap<String, THashMap<String, Integer>>();
        this.matriceJunction = new THashMap<String, THashMap<String, Integer>>();
        this.geneMetrics = new THashMap<String, GeneMetrics>();
        this.molecules = new THashSet<Molecule>();
        this.cellMetrics = new THashMap<String, CellMetrics>();
        Iterator it = cells.iterator();
        while (it.hasNext())
            cellMetrics.put((String)it.next(), new CellMetrics());
    }

    public THashMap<String, THashMap<String, Integer>> getMatrice(){ return matrice; }
    public THashMap<String, THashMap<String, Integer>> getMatriceGene(){ return matriceGene; }
    public THashMap<String, THashMap<String, Integer>> getMatriceJunction(){ return matriceJunction; }
    public THashMap<String, CellMetrics> getCellMetrics(){ return cellMetrics; }
    public THashMap<String, GeneMetrics> getGeneMetrics(){ return geneMetrics; }
    public THashSet<Molecule> getMolecules(){ return molecules; }
    
    public int getNb_count(){ return nb_count; }
    public int getTotal_count(){ return total_count; }
    public int getTotal_remove(){ return total_remove; }
    public int getTotal_isoform_def(){ return total_isoform_def; }
    public int getTotal_isoform_undef(){ return total_isoform_undef; } 
       
    public void addMolecule(Molecule molecule)
    {
        THashMap mapCell = null;
        THashSet umiSet = null;
        
        // only produce matrix for cells that we need to consider
        // we may have some molecules from other cell barcodes
        if(cellMetrics.containsKey(molecule.getBarcode())){
            this.molecules.add(molecule);
            //System.out.println(molecule);
        
            if(! geneMetrics.containsKey(molecule.getGeneId()))
                geneMetrics.put(molecule.getGeneId(), new GeneMetrics());
        
            //System.out.println(molecule.getBarcode() + "|" + molecule.getUmi()+ "|" + molecule.getGeneId()+ "|" + molecule.getTranscriptId()+ "|" + molecule.getLongreads());
        
            ((CellMetrics)cellMetrics.get(molecule.getBarcode())).addCount(molecule.getGeneId(), molecule.getTranscriptId(), molecule.getLongreads().size());
            ((GeneMetrics)geneMetrics.get(molecule.getGeneId())).addCount(molecule.getGeneId(), molecule.getTranscriptId());
            //((JunctionMetrics)junctionMetrics.get(molecule.getGeneId())).addCount(molecule.getGeneId(), molecule.getTranscriptId());

            if ("undef".equals(molecule.getTranscriptId()))
                total_isoform_undef++;
            else
                total_isoform_def++;

            // we already have this gene/transcript isokey
            String isokey = molecule.getGeneId() + "\t" + molecule.getTranscriptId();
            if ((mapCell = (THashMap) matrice.get(isokey)) != null) {
                // we already have this cell
                if ((umiSet = (THashSet) mapCell.get(molecule.getBarcode())) != null) {
                    umiSet.add(molecule.getUmi());
                }
                else {
                    umiSet = new THashSet();
                    umiSet.add(molecule.getUmi());
                    mapCell.put(molecule.getBarcode(), umiSet);
                }
            }
            else {
                matrice.put(isokey, new THashMap());
                umiSet = new THashSet();
                umiSet.add(molecule.getUmi());
                ((THashMap) matrice.get(isokey)).put(molecule.getBarcode(), umiSet);
            }

            // now work on matrixGene
            isokey = molecule.getGeneId();
            if ((mapCell = (THashMap) matriceGene.get(isokey)) != null) {
                // we already have this cell
                if ((umiSet = (THashSet) mapCell.get(molecule.getBarcode())) != null) {
                    umiSet.add(molecule.getUmi());
                }
                else {
                    umiSet = new THashSet();
                    umiSet.add(molecule.getUmi());
                    mapCell.put(molecule.getBarcode(), umiSet);
                }
            }
            else {
                matriceGene.put(isokey, new THashMap());
                umiSet = new THashSet();
                umiSet.add(molecule.getUmi());
                ((THashMap) matriceGene.get(isokey)).put(molecule.getBarcode(), umiSet);
            }

            // now work on matrixJunction
            HashSet<Junction> junctionSet = molecule.getJunctionSet();
            Iterator<Junction> it = junctionSet.iterator();
            while(it.hasNext()) {
                Junction j = it.next();
                String junckey = molecule.getGeneId() + ":" + j.getStart() + "-" + j.getEnd();

                // we already have this junction
                if ((mapCell = (THashMap) matriceJunction.get(junckey)) != null) {
                    if ((umiSet = (THashSet) mapCell.get(molecule.getBarcode())) != null) {
                        umiSet.add(molecule.getUmi());
                    }
                    else {
                        umiSet = new THashSet();
                        umiSet.add(molecule.getUmi());
                        mapCell.put(molecule.getBarcode(), umiSet);
                    }
                }
                else {
                    matriceJunction.put(junckey, new THashMap());
                    umiSet = new THashSet();
                    umiSet.add(molecule.getUmi());
                    ((THashMap) matriceJunction.get(junckey)).put(molecule.getBarcode(), umiSet);
                }
            }
        }
        //else{
        //    log.info(new Object[]{"cellMetrics does not contains cellBC\t"+molecule.getBarcode()});
        //}
    }

    public void writeIsoformMatrix(java.io.File isomatrix, java.io.File isometrics, java.io.File molmetrics, UCSCRefFlatParser model)
    {
        BufferedOutputStream os = null;
        BufferedOutputStream os2 = null;
        BufferedOutputStream os3 = null;
        THashSet setUmi = null;
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(isomatrix));
            os2 = new BufferedOutputStream(new java.io.FileOutputStream(isometrics));
            os3 = new BufferedOutputStream(new java.io.FileOutputStream(molmetrics));
            
            os.write(new String("geneId\ttranscriptId\tnbExons").getBytes());
            os2.write(new String("geneId\ttranscriptId\tnbExons\tnbUmis\n").getBytes());
            
            for(String key : cellMetrics.keySet())
                os.write(new String("\t" + key).getBytes());
               
            os.write(new String("\n").getBytes());

            for(String isokey : matrice.keySet()){
                String[] tmp = isokey.split("\t");
                
                if(model != null){
                    TranscriptRecord tr = (TranscriptRecord)model.select(tmp[0],tmp[1]);
                    int nb_exon = 0;
                    if(tr != null)
                        nb_exon = tr.getExons().size();
                    
                    os.write(new String(isokey + "\t" + nb_exon).getBytes());
                    os2.write(new String(isokey + "\t" + nb_exon).getBytes());
                }
                else{
                    os.write(new String(isokey+"\tna").getBytes());
                    os2.write(new String(isokey+"\tna").getBytes());
                }
                    
                int total = 0;
                for(String cell_barcode : cellMetrics.keySet()){
                    if ((setUmi = (THashSet) ((THashMap) matrice.get(isokey)).get(cell_barcode)) != null) {
                        os.write(new String("\t" + setUmi.size()).getBytes());
                        this.total_count += setUmi.size();
                        total += setUmi.size();
                    }
                    else
                        os.write(new String("\t0").getBytes());
                }
                os.write(new String("\n").getBytes());
                os2.write(new String("\t" + total + "\n").getBytes());
            }
            
            os3.write(new String("cellBC\tUMI\tnbReads\tnbSupportingReads\tmappingPctId\tsnpPhredScore\tgeneId\ttranscriptId\n").getBytes());
            Iterator<Molecule> it = this.molecules.iterator();
            while(it.hasNext()) {
                Molecule m = it.next();
                os3.write(new String(m.getBarcode() + "\t" + m.getUmi() + "\t" + m.getNumberOfReads() + "\t" + m.getSupporting_reads() + "\t" + m.getPctId() + "\t" + m.getSnpPhredScore() + "\t" + m.getGeneId() + "\t" + m.getTranscriptId() + "\n").getBytes());
            }
            
            os.close();
            os2.close();
            os3.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close();  os2.close(); os3.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close();  os2.close(); os3.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }

    public void writeJunctionMatrix(java.io.File juncmatrix, java.io.File juncmetrics)
    {
        BufferedOutputStream os = null;
        BufferedOutputStream os2 = null;
        THashSet umiSet = null;
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(juncmatrix));
            os2 = new BufferedOutputStream(new java.io.FileOutputStream(juncmetrics));
            
            os.write(new String("junctionId").getBytes());
            os2.write(new String("junctionId\tnbUmis\n").getBytes());
            for(String key : cellMetrics.keySet())
                os.write(new String("\t" + key).getBytes());
            os.write(new String("\n").getBytes());

            for(String junckey : matriceJunction.keySet()){
                os.write(junckey.getBytes());
                os2.write(junckey.getBytes());
                int total = 0;
                for(String cell_barcode : cellMetrics.keySet()){
                    if ((umiSet = (THashSet) ((THashMap) matriceJunction.get(junckey)).get(cell_barcode)) != null) {
                        os.write(new String("\t" + umiSet.size()).getBytes());
                        //this.total_count += umiSet.size();
                        total += umiSet.size();
                    }
                    else
                        os.write(new String("\t0").getBytes());
                }
                os.write(new String("\n").getBytes());
                os2.write(new String("\t" + total + "\n").getBytes());
            }
            
            os.close();
            os2.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close();  os2.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close();  os2.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
    
    public void writeGeneMatrix(java.io.File genematrix, java.io.File genemetrics)
    {
        BufferedOutputStream os = null;
        BufferedOutputStream os2 = null;
        THashSet setUmi = null;
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(genematrix));
            os.write(new String("geneId").getBytes());
            for(String cell_barcode : cellMetrics.keySet())
                os.write(new String("\t" + cell_barcode).getBytes());
            os.write(new String("\n").getBytes());

            for(String isokey : matriceGene.keySet()){
                os.write(isokey.getBytes());
                for(String cell_barcode : cellMetrics.keySet()){
                    if ((setUmi = (THashSet) ((THashMap) matriceGene.get(isokey)).get(cell_barcode)) != null)
                        os.write(new String("\t" + setUmi.size()).getBytes());
                    else
                        os.write(new String("\t0").getBytes());
                }
                os.write(new String("\n").getBytes());
            }
            
            os2 = new BufferedOutputStream(new java.io.FileOutputStream(genemetrics));
            os2.write(new String("geneId\tnbUmis\tnbIsoformSet\tnbIsoformNotSet\n").getBytes());
            for(String geneId : geneMetrics.keySet())
            {
                GeneMetrics gm = this.geneMetrics.get(geneId);
                if(gm != null){
                    int total = gm.getIsoform_known_count() + gm.getIsoform_undef_count();
                    os2.write(new String(geneId+"\t"+total+"\t"+gm.getIsoform_known_count()+"\t" +gm.getIsoform_undef_count()+"\n").getBytes());
                }
                else{
                    os2.write(new String(geneId+"\t0\t0\t0\n").getBytes());
                }
            }
            
            os.close();
            os2.close();
        } catch (Exception e) { e.printStackTrace(); try { os.close();  os2.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close();  os2.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }

    public void writeCellMetrics(java.io.File paramFile)
    {
        BufferedOutputStream os = null;
        
        try {
            os = new BufferedOutputStream(new java.io.FileOutputStream(paramFile));
            os.write(new String("cellBC\tnbReads\tnbUmis\tnbIsoformSet\tnbIsoformNotSet\n").getBytes());
            
            for(String cell_barcode : cellMetrics.keySet()){
                CellMetrics cm = cellMetrics.get(cell_barcode);
                if(cm != null){
                    //int total = cm.getIsoform_known_count() + cm.getIsoform_undef_count();
                    os.write(new String(cell_barcode+"\t"+cm.getNb_reads()+"\t"+cm.getNb_umis()+"\t"+cm.getIsoform_known_count()+"\t" +cm.getIsoform_undef_count()+"\n").getBytes());
                }
                else{
                    os.write(new String(cell_barcode+"\t0\t0\t0\t0\t0\n").getBytes());
                }
            }
            
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try { os.close(); } catch (Exception e2) { System.err.println("can not close stream"); }
        } finally { try { os.close(); } catch (Exception e3) { System.err.println("can not close stream");  } }
    }
}