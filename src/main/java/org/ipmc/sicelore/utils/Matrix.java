package org.ipmc.sicelore.utils;

import java.util.*;
import java.io.*;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import htsjdk.samtools.util.Log;

public class Matrix {

    private final Log log;
    private HashMap<String, HashMap<String, Integer>> matrice;
    private HashMap<String, CellMetrics> cellMetrics;

    private HashSet<String> molecules;
    private HashSet<String> cells;
    private HashSet<String> bcumi;

    private int nb_count = 0;
    private int total_count = 0;
    private int total_remove = 0;
    private int total_isoform_def = 0;
    private int total_isoform_undef = 0;
    
    private boolean is_gene_level = false;
    
    public Matrix(boolean is_gene_level)
    {
        log = Log.getInstance(Matrix.class);
        this.matrice = new HashMap<String, HashMap<String, Integer>>();
        this.cellMetrics = new HashMap<String, CellMetrics>();
        
        this.cells = new HashSet<String>();
        this.bcumi = new HashSet<String>();
        this.is_gene_level = is_gene_level;
    }

    public void addMolecule(Molecule molecule)
    {
        HashMap mapCell = null;
        HashSet setUmi = null;
        
        String molkey = molecule.getBarcode() + ":" + molecule.getUmi();
        String isokey = molecule.getGeneId() + "\t" + molecule.getTranscriptId();
        if(this.is_gene_level)
            isokey = molecule.getGeneId() + "\tundef";
            
        if (this.bcumi.contains(molkey)){
            System.out.println(molkey + " already in matrix !!!");
        }

        this.bcumi.add(molkey);
        
        if(! cellMetrics.containsKey(molecule.getBarcode())){    
            cells.add(molecule.getBarcode());
            cellMetrics.put(molecule.getBarcode(), new CellMetrics());
        }
        ((CellMetrics)cellMetrics.get(molecule.getBarcode())).addCount(molecule.getGeneId(),molecule.getTranscriptId());
        
        if ("undef".equals(molecule.getTranscriptId())) {
            total_isoform_undef++;
        } else {
            total_isoform_def++;
        }

        // we already have this gene/transcript isokey
        if ((mapCell = (HashMap) matrice.get(isokey)) != null) {
            // we already have this cell
            if ((setUmi = (HashSet) mapCell.get(molecule.getBarcode())) != null) {
                setUmi.add(molecule.getUmi());
            } else {
                setUmi = new HashSet();
                setUmi.add(molecule.getUmi());
                mapCell.put(molecule.getBarcode(), setUmi);
            }
        } else {
            matrice.put(isokey, new HashMap());
            setUmi = new HashSet();
            setUmi.add(molecule.getUmi());
            ((HashMap) matrice.get(isokey)).put(molecule.getBarcode(), setUmi);
        }
    }

    public void writeMatrix(java.io.File paramFile, HashSet<String> authorizedCells)
    {
        DataOutputStream os = null;
        HashSet setUmi = null;
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(paramFile));
            os.writeBytes("geneId\ttranscriptId");
            Iterator<String> itcell = authorizedCells.iterator();
            while (itcell.hasNext()) {
                String cell_barcode = (String) itcell.next();
                os.writeBytes("\t" + cell_barcode);
            }
            os.writeBytes("\n");

            Set cles = this.matrice.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String isokey = (String) it.next();
                os.writeBytes(isokey);

                Iterator itcell2 = authorizedCells.iterator();
                while (itcell2.hasNext()) {
                    String cell_barcode = (String) itcell2.next();

                    if ((setUmi = (HashSet) ((HashMap) matrice.get(isokey)).get(cell_barcode)) != null) {
                        os.writeBytes("\t" + setUmi.size());
                        this.total_count += setUmi.size();
                    } else {
                        os.writeBytes("\t0");
                    }
                }
                os.writeBytes("\n");
            }
            os.close();
        } catch (Exception e) {
            e.printStackTrace();
            try {
                os.close();
            } catch (Exception e2) {
                System.err.println("can not close stream");
            }
        } finally {
            try {
                os.close();
            } catch (Exception e3) {
                System.err.println("can not close stream");
            }
        }
    }

    public void writeMetrics(java.io.File paramFile, HashSet<String> authorizedCells)
    {
        DataOutputStream os = null;
        
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(paramFile));
            os.writeBytes("cell\tisoform_total\tisoform_known\tisoform_undef\n");
            Iterator<String> itcell = authorizedCells.iterator();
            while (itcell.hasNext()) {
                String cell_barcode = (String) itcell.next();
                CellMetrics m = this.cellMetrics.get(cell_barcode);
                if(m != null){
                    int total = m.getIsoform_known_count() + m.getIsoform_undef_count();
                    os.writeBytes(cell_barcode+"\t"+total+"\t"+m.getIsoform_known_count()+"\t" +m.getIsoform_undef_count()+"\n");
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


        log.info(new Object[]{"\t\tMatrix cells\t\t[" + cells.size() + "]"});
        log.info(new Object[]{"\t\tMatrix gene isoforms\t[" + matrice.size() + "]"});
        log.info(new Object[]{"\t\tMatrix total counts\t[" + total_count + "]"});
        log.info(new Object[]{"\t\tMatrix isoform found\t[" + total_isoform_def + "]"});
        log.info(new Object[]{"\t\tMatrix isoform unknown\t[" + total_isoform_undef + "]"});
    }
}
