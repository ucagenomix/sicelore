package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.Matrix;
import org.ipmc.sicelore.utils.Molecule;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Get Isoforms-2-Editing association for Gria2 molecules", oneLineSummary = "Get Isoforms-2-Editing association for Gria2 molecules", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class GetEditingLabel extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    //@Argument(shortName = "O", doc = "The output directory")
    //public File OUTPUT;
    //@Argument(shortName = "CSV", doc = "The .csv cluster file")
    //public File CSV;
    //@Argument(shortName = "CELL_FLAG", doc = "The cell barcode flag (default CB)")
    //public String CELL_FLAG = "CB";

    public HashSet<String> DTEcells;
    
    public GetEditingLabel() {
        log = Log.getInstance(GetEditingLabel.class);
        pl = new ProgressLogger(log);
        this.DTEcells = new HashSet<String>();
    }

    protected int doWork() {
        String str1 = null;

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);

        loadDTEcells();        
        Matrix matrix = new Matrix(DTEcells);
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        //htsjdk.samtools.SAMFileHeader localSAMFileHeader = samReader.getFileHeader();
        //SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            for(SAMRecord r : samReader) {
                pl.record(r);
                String readName = r.getReadName();
                String readString = r.getReadString();
                //String US = (String)r.getAttribute("US");
                String IG = (String)r.getAttribute("IG");
                String BC = (String)r.getAttribute("BC");
                String U8 = (String)r.getAttribute("U8");
                
                int RGpos = r.getReadPositionAtReferencePosition(80692286);
                int QRpos = r.getReadPositionAtReferencePosition(80706912);
                
                boolean is_short = ((int)r.getReadPositionAtReferencePosition(80692000)>0)?true:false;
                boolean is_flop  = ((int)r.getReadPositionAtReferencePosition(80691420)>0)?true:false;
                boolean is_flip  = ((int)r.getReadPositionAtReferencePosition(80690470)>0)?true:false;
                boolean is_exon  = ((int)r.getReadPositionAtReferencePosition(80691720)>0)?true:false;
                
                if(RGpos > 0 && readString.length() > RGpos){
                    String RGbase = readString.substring(RGpos - 1, RGpos);
                    String QRbase = "NA";
                    if(QRpos > 0){ QRbase = readString.substring(QRpos - 1, QRpos); }
                    
                    String transcriptId = "undef";
                    if(is_short){ transcriptId = "short"; }
                    else if(is_flop && !is_flip && !is_exon){ transcriptId = "flop"; }
                    else if(!is_flop && is_flip && !is_exon){ transcriptId = "flip";}
                    else if(is_flop && is_flip && !is_exon){ transcriptId = "flipflop"; }
                    else if(is_flop && is_flip && is_exon){ transcriptId = "flipflopexon"; }
                    
                    transcriptId = transcriptId + "." + RGbase + "." + QRbase;
                    
                    //System.out.println(readName+","+IG+"\t"+transcriptId);
                    
                    Molecule molecule = new Molecule(BC,U8);
                    molecule.setGeneId("Gria2");
                    molecule.setTranscriptId(transcriptId);
                    matrix.addMolecule(molecule);
                }
            }
            samReader.close();

        } catch (Exception e) { e.printStackTrace(); }

        matrix.writeIsoformMatrix(new File("./gria2_matrix.txt"), new File("./gria2_metrics.txt"));

        return 0;
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

    public static void main(String[] paramArrayOfString) {
        System.exit(new SplitBamPerCluster().instanceMain(paramArrayOfString));
    }
}
