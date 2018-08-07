package org.ipmc.sicelore.programs;

import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Collapse Isoforms Expression Matrix to Genes Expression Matrix", oneLineSummary = "Collapse Isoforms Expression Matrix to Genes Expression Matrix", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class CollapseToGeneExpressionMatrix extends CommandLineProgram {

    private final Log log;
    private htsjdk.samtools.util.ProgressLogger pl;
    @Argument(shortName = "I", doc = "The iput Isoforms Expression Matrix")
    public java.io.File INPUT;
    @Argument(shortName = "MATRIX", doc = "The output collapsed Gene Expression Matrix")
    public java.io.File MATRIX;

    public CollapseToGeneExpressionMatrix() {
        log = Log.getInstance(CollapseToGeneExpressionMatrix.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
    }

    protected int doWork() {
        Integer value;
        Integer k;
        String str = null;
        DataOutputStream os = null;

        HashMap<String, HashMap<String, Integer>> matrix = new HashMap<String, HashMap<String, Integer>>();
        HashSet<String> genes = new HashSet<String>();
        HashSet<String> cells = new HashSet<String>();
        Vector vectorGenes = new Vector();
        Vector vectorCells = new Vector();

        IOUtil.assertFileIsReadable(INPUT);

        try {
            BufferedReader fichier = new BufferedReader(new FileReader(INPUT));
            String[] cellbc = fichier.readLine().split("\t");
            String line = fichier.readLine();
            while (line != null) {
                String[] dat = line.split("\t");

                if (!genes.contains(dat[0])) {
                    genes.add(dat[0]);
                    vectorGenes.addElement(dat[0]);
                    matrix.put(dat[0], new HashMap());
                }

                for (int i = 2; i < dat.length; i++) {
                    if ((value = (Integer) ((HashMap) matrix.get(dat[0])).get(cellbc[i])) == null) {
                        ((HashMap) matrix.get(dat[0])).put(cellbc[i], new Integer(dat[i]));
                    } else {
                        ((HashMap) matrix.get(dat[0])).put(cellbc[i], new Integer(value.intValue() + new Integer(dat[i]).intValue()));
                    }
                }
                line = fichier.readLine();
            }
            fichier.close();

            os = new DataOutputStream(new java.io.FileOutputStream(MATRIX));
            os.writeBytes("geneId");
            for (int i = 2; i < cellbc.length; i++) {
                os.writeBytes("\t" + cellbc[i]);
            }
            os.writeBytes("\n");

            for (int i = 0; i < vectorGenes.size(); i++) {
                String g = (String) vectorGenes.get(i);
                os.writeBytes(g);
                for (int j = 2; j < cellbc.length; j++) {
                    if ((k = ((HashMap<String, Integer>) matrix.get(g)).get(cellbc[j])) != null) {
                        os.writeBytes("\t" + k.intValue());
                    } else {
                        os.writeBytes("\t0");
                    }
                }
                os.writeBytes("\n");
            }
            os.close();

        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            try {
                os.close();
            } catch (Exception e) {
                System.err.println("can not close stream");
            }
        }

        return 0;

    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new CollapseToGeneExpressionMatrix().instanceMain(paramArrayOfString));
    }
}
