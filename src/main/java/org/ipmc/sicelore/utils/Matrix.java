package org.ipmc.sicelore.utils;

import java.util.*;
import java.io.*;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.StringUtils;
import htsjdk.samtools.util.Log;

public class Matrix {

    private final Log log;
    private HashMap<String, HashMap<String, Integer>> matrice;

    private HashSet<String> molecules;
    private HashSet<String> cells;
    private HashSet<String> bcumi;

    private int nb_count = 0;
    private int total_count = 0;
    private int total_remove = 0;
    private int total_isoform_def = 0;
    private int total_isoform_undef = 0;

    public Matrix() {
        log = Log.getInstance(Matrix.class);
        this.matrice = new HashMap<String, HashMap<String, Integer>>();
        this.cells = new HashSet<String>();
        this.bcumi = new HashSet<String>();
    }

    public void addMolecule(Molecule molecule) {
        HashMap mapCell = null;
        HashSet setUmi = null;
        String isokey = molecule.getGeneId() + "\t" + molecule.getTranscriptId();
        String molkey = molecule.getBarcode() + ":" + molecule.getUmi();

        if (this.bcumi.contains(molkey)) {
            System.out.println(molkey + " already in matrix !!!");
        }

        this.bcumi.add(molkey);

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
                cells.add(molecule.getBarcode());
                mapCell.put(molecule.getBarcode(), setUmi);
            }
        } else {
            matrice.put(isokey, new HashMap());
            setUmi = new HashSet();
            setUmi.add(molecule.getUmi());
            ((HashMap) matrice.get(isokey)).put(molecule.getBarcode(), setUmi);
        }
    }

    public void write(java.io.File paramFile, HashSet<String> authorizedCells) {
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

        this.displayMetrics();
    }

    public void displayMetrics() {
        log.info(new Object[]{"\t\tMatrix cells\t\t[" + cells.size() + "]"});
        log.info(new Object[]{"\t\tMatrix gene isoforms\t[" + matrice.size() + "]"});
        log.info(new Object[]{"\t\tMatrix total counts\t[" + total_count + "]"});
        log.info(new Object[]{"\t\tMatrix isoform found\t[" + total_isoform_def + "]"});
        log.info(new Object[]{"\t\tMatrix isoform unknown\t[" + total_isoform_undef + "]"});
        /*
    	
    	log.info(new Object[] { "\tMatrix molecules\t[" + molecules.size() + "]"});
    	log.info(new Object[] { "\tMatrix molecule2isoform\t[" + molecule2isoform.size() + "]"});
    	
    	int nb=0;
    	int ok=0;
        Iterator it = molecule2isoform.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            String molkey = (String)pair.getKey();
            String isokey = (String)pair.getValue();
        	String[] dat = molkey.split("_");
        	
        	
        	if(! this.features.contains(isokey)){
        		System.out.println(isokey + " not found"); 
        	}
        	if(! this.cells.contains(dat[0])){
        		System.out.println(dat[0] + " not found"); 
        	}
        	
        	if(((HashSet)((HashMap)matrice.get(isokey)).get(dat[0])).contains(dat[1])){
        		ok++;
        		System.out.println(ok + "\tin\t" + isokey + "\t" + dat[0] + "\t" + dat[1]);
        	}
        	else{
        		nb++;
        		System.out.println(nb + "\tout\t" + isokey + "\t" + dat[0] + "\t" + dat[1]);
        	}
        }
         */
    }
}
