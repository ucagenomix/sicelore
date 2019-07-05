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

@CommandLineProgramProperties(summary = "UMI depth histogram calculation", oneLineSummary = "UMI depth histogram calculation", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class UMIDepthHisto extends CommandLineProgram
{
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "T", doc = "The type of input SAM or BAM file (shortread, longread or molecule)")
    public String TYPE = "molecule";
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "OUT", doc = "The metrics output .csv file")
    public File OUT;

    private final Log log;
    private ProgressLogger pl;
    HashMap<Integer, Integer> counts;
    public HashSet<String> DTEcells;

    public UMIDepthHisto() {
        log = Log.getInstance(UMIDepthHisto.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log);
        this.counts = new HashMap<Integer, Integer>();
        this.DTEcells = new HashSet<String>();
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        process();

        return 0;
    }

    protected void process()
    {
        loadDTEcells();
        HashMap<String, Integer> m = new HashMap<String, Integer>();
        
        if("longread".equals(TYPE)){
            log.info(new Object[]{"loadLongreadsBam"});
            
            LongreadParser bam = new LongreadParser(INPUT, false, false);
            MoleculeDataset dataset = new MoleculeDataset(bam);
            
            Set cles = dataset.getMapMolecules().keySet();
            Iterator it = cles.iterator();
            while(it.hasNext()){
                String molkey = (String) it.next();
                Molecule molecule = (Molecule)dataset.getMapMolecules().get(molkey);
                Integer nbreads = molecule.getLongreads().size();
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
        log.info(new Object[]{"loadMoleculesBam"});
        
        HashMap<String, Integer> m = new HashMap<String, Integer>();
        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String name = r.getReadName();
                String[] tab = name.split("|");
                
                m.put(tab[2]+"|"+tab[3], new Integer(tab[4]));
            }
            localSamReader.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        return m;
    }
    
    public HashMap<String, Integer> loadShortreadsBam()
    {
        log.info(new Object[]{"loadShortreadsBam"});
        HashMap<String, Integer> m = new HashMap<String, Integer>();
        HashMap<String, HashSet<String>> mr = new HashMap<String, HashSet<String>>();
        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        HashSet<String> value;
        
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String name = r.getReadName();
                String CB = (String)r.getAttribute("CB");
                if(CB != null)
                    CB = CB.replace("-1", "");
                
                String UB = (String)r.getAttribute("UB");
                String GN = (String)r.getAttribute("GN");
                
                if(GN != null && CB != null && UB != null && this.DTEcells.contains(CB)){
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
    
    public void loadDTEcells()
    {
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(CSV));
            String line = fichier.readLine();
            while (line != null) {
                DTEcells.add(line);
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
     public void writeUMIDepthHisto(HashMap<String, Integer> hash)
     {
        DataOutputStream os = null;
        int max=0;
        Integer value;
        
        log.info(new Object[]{"\tRunning through molecules..."});
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
        log.info(new Object[]{"\tSetIsoforms\t\tend..."});
        
        try {
            os = new DataOutputStream(new java.io.FileOutputStream(OUT));
            os.writeBytes("UMIdepth\tNbMolecules\n");
            for(Integer depth : counts.keySet())
                os.writeBytes(depth + "\t" + counts.get(depth) + "\n");

            os.close();
        } catch (Exception e) { e.printStackTrace(); }
        finally { try { os.close(); } catch (Exception e) { System.err.println("can not close stream");  } }
    }

    public static void main(String[] args) {
        System.exit(new UMIDepthHisto().instanceMain(args));
    }
}
