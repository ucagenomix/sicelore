package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import java.util.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Compute UMI depth histogram for shortread, longread or molecule bam file.", oneLineSummary = "Compute UMI depth histogram for shortread, longread or molecule bam file.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class HistoUMIDepth extends CommandLineProgram
{
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "O", doc = "The metrics output .csv file")
    public File OUTPUT;
    @Argument(shortName = "T", doc = "Type of input SAM or BAM file (shortread, longread or molecule)")
    public String TYPE = "molecule";

    private final Log log;
    private ProgressLogger pl;
    HashMap<Integer, Integer> counts;
    public CellList cellList;

    public HistoUMIDepth() {
        log = Log.getInstance(HistoUMIDepth.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
        this.counts = new HashMap<Integer, Integer>();
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        process();

        return 0;
    }

    protected void process()
    {
        this.cellList = new CellList(CSV);
        log.info(new Object[]{"Cells detected\t\t\t[" + this.cellList.size() + "]"});

        HashMap<String, Integer> m = new HashMap<String, Integer>();
        
        if("longread".equals(TYPE)){
            log.info(new Object[]{"\tloadLongreadsBam"});
            
            LongreadParser bam = new LongreadParser(INPUT, false, false, true, true);
            MoleculeDataset dataset = new MoleculeDataset(bam);
            
            Set cles = dataset.getMapMolecules().keySet();
            Iterator it = cles.iterator();
            while(it.hasNext()){
                String molkey = (String) it.next();
                Molecule molecule = (Molecule)dataset.getMapMolecules().get(molkey);
                Integer nbreads = molecule.getLongreads().size();
                
                if(this.cellList.contains(molecule.getBarcode()))
                    m.put(molecule.getBarcode()+"|"+molecule.getUmi(), nbreads);
            }
        }
        else if("molecule".equals(TYPE)){
            m = loadMoleculesBam();
        }
        else if("shortread".equals(TYPE)){
            m = loadShortreadsBam();
        }
        
        log.info(new Object[]{"\tMolecules\t\t[" + m.size() + "]"});
        this.writeUMIDepthHisto(m);
    }
    
    public HashMap<String, Integer> loadMoleculesBam()
    {
        log.info(new Object[]{"\tloadMoleculesBam"});
        
        HashMap<String, Integer> m = new HashMap<String, Integer>();
        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String name = r.getReadName();
                String BC = (String)r.getAttribute("BC");
                if(BC != null)
                    BC = BC.replace("-1", "");
                
                String U8 = (String)r.getAttribute("U8");
                int rn = (Integer) r.getAttribute("RN");
                
                if(this.cellList.contains(BC))
                    m.put(BC+"|"+U8, rn);
            }
            localSamReader.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        return m;
    }
    
    public HashMap<String, Integer> loadShortreadsBam()
    {
        log.info(new Object[]{"\tloadShortreadsBam"});
        HashMap<String, Integer> m = new HashMap<String, Integer>();
        HashMap<String, HashSet<String>> mr = new HashMap<String, HashSet<String>>();
        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        HashSet<String> value;
        
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String name = r.getReadName();
                String UB = (String)r.getAttribute("UB");
                String CB = (String)r.getAttribute("CB");
                if(CB != null)
                    CB = CB.replace("-1", "");
                
                if(CB != null && UB != null && this.cellList.contains(CB)){
                    if((value = (HashSet<String>)mr.get(CB+"|"+UB)) != null)
                        value.add(name);
                    else{
                        value = new HashSet<String>();
                        value.add(name);
                    }
                    mr.put(CB+"|"+UB, value);
                }
            }
            localSamReader.close();
            
            Set cles = mr.keySet();
            Iterator it = cles.iterator();
            while (it.hasNext()) {
                String molkey = (String)it.next();
                Integer nbreads = ((HashSet<String>)mr.get(molkey)).size();
                m.put(molkey, nbreads);
            }
        } catch (Exception e) { e.printStackTrace(); }
        
        return m;
    }    
    
    public void writeUMIDepthHisto(HashMap<String, Integer> hash)
    {
        DataOutputStream os = null;
        int max=0;
        Integer value;
        
        log.info(new Object[]{"\twriteUMIDepthHisto"});
        Set cles = hash.keySet();
        Iterator it = cles.iterator();
        while (it.hasNext()) {
            String molkey = (String) it.next();
            Integer nbreads = (Integer)hash.get(molkey);
            
            if((value = (Integer)counts.get(nbreads)) != null)
                value += 1;
            else
                value = 1;
            
            counts.put(nbreads, value);

            if(max<nbreads)
                max=nbreads;
        }
        
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(OUTPUT));
            os.writeBytes("depth\tnumber\n");
            for(Integer depth : counts.keySet())
                os.writeBytes(depth + "\t" + counts.get(depth) + "\n");

            os.close();
        } catch (Exception e) { e.printStackTrace(); }
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream");  } }
    }

    public static void main(String[] args) {
        System.exit(new HistoUMIDepth().instanceMain(args));
    }
}
