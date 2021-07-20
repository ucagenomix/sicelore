package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.CellList;
import java.io.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.Set;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.LongreadParser;
import org.ipmc.sicelore.utils.LongreadRecord;
import org.ipmc.sicelore.utils.Matrix;
import org.ipmc.sicelore.utils.Molecule;
import org.ipmc.sicelore.utils.MoleculeDataset;
import picard.cmdline.CommandLineProgram;
import gnu.trove.THashMap;

@CommandLineProgramProperties(summary = "Detect fusion genes.", oneLineSummary = "Detect fusion genes.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class FusionDetector extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "PREFIX", doc = "The output file prefix (default=fusion)")
    public String PREFIX = "fusion";
   
    private ProgressLogger pl;
    private final Log log;
    
    public CellList cellList;

    public FusionDetector() {
        log = Log.getInstance(FusionDetector.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        
        this.cellList = new CellList(CSV);
        log.info(new Object[]{"\tCells detected\t[" + this.cellList.size() + "]"});
        Matrix matrix = new Matrix(cellList);
        
	LongreadRecord lrr = new LongreadRecord();
	lrr.setStaticParams("BC","U8","GE","TE","UE","PE","US",10000,"RN");

        // umi not mandatory cos' 2ndary fusion SAM records without UMI !
        LongreadParser bam = new LongreadParser(INPUT, false, false, true, false);
        MoleculeDataset dataset = new MoleculeDataset(bam);
        
        log.info(new Object[]{"\tSetFusions\t\tstart..."});
        
        Map<String, Integer> count = new HashMap<String, Integer>();
        THashMap<String, Molecule> mapMolecules = dataset.getMapMolecules();
        Set cles = mapMolecules.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            Molecule molecule = (Molecule) mapMolecules.get(molkey);
            HashSet<String> list = molecule.getGeneIds();
            
            if(this.cellList.contains(molecule.getBarcode()) && molecule.getUmi() != null && list.size() == 2){
                String key = list.toString();
                key = key.replace(", ", "|");
                key = key.replace("[", "");
                key = key.replace("]", "");
                int c = count.containsKey(key) ? count.get(key) : 0;
                count.put(key, c + 1);
                molecule.setGeneId(key);
                molecule.setTranscriptId(key);
                matrix.addMolecule(molecule);
            }
        }
        
        HashMap<String, Integer> sorted = count.entrySet().stream().sorted(Collections.reverseOrder(Map.Entry.comparingByValue()))
            .collect(java.util.stream.Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue, (e1, e2) -> e2, LinkedHashMap::new));
        
        cles = sorted.keySet();
        it = cles.iterator();
        while (it.hasNext()) {
            String key = (String) it.next();
            int c = (Integer)sorted.get(key);
            if(c >= 10){
                log.info(new Object[]{"\t" + c + " distincts molecules support fusion [" + key + "]"});
            }
        }

        File MATRIX = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_fusmatrix.txt");
        File METRICS = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_fusmetrics.txt");
        File molinfos = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_fusmolinfos.txt");
        matrix.writeIsoformMatrix(MATRIX, METRICS, molinfos, null);
        
        return 0;
    }

    public static void main(String[] args) {
        System.exit(new FusionDetector().instanceMain(args));
    }
}
