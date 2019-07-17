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
import org.ipmc.sicelore.utils.LongreadRecord;
import org.ipmc.sicelore.utils.Longread;
import org.ipmc.sicelore.utils.Matrix;
import org.ipmc.sicelore.utils.Molecule;
import org.ipmc.sicelore.utils.TranscriptRecord;
import org.ipmc.sicelore.utils.UCSCRefFlatParser;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SNP/editing events detection/quantification (cellBC/UMI count matrix).", oneLineSummary = "SNP/editing event detection/quantification (cellBC/UMI count matrix).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SNPMatrix extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input molecule SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "SNP", doc = "The SNP/editing events .csv file \n#-----\nname,gene,chromosome,start,end,strand,pos,ref,alt\nRG,Gria2,3,80681450,80802835,-,80692286,T,C\nQR,Gria2,3,80681450,80802835,-,80706912,T,C\n#-----")
    public File SNP;
    //@Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    //public File REFFLAT;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "PREFIX", doc = "The output file prefix (default=snp)")
    public String PREFIX = "snp";
    //@Argument(shortName = "CELL_FLAG", doc = "The cell barcode flag (default CB)")
    //public String CELL_FLAG = "CB";

    public HashSet<String> DTEcells;
    
    public SNPMatrix() {
        log = Log.getInstance(SNPMatrix.class);
        pl = new ProgressLogger(log);
        this.DTEcells = new HashSet<String>();
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        IOUtil.assertFileIsReadable(SNP);
        //IOUtil.assertFileIsReadable(REFFLAT);

        loadDTEcells();        
        Matrix matrix = new Matrix(DTEcells);
        
        //UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
        //String[] gria2 = {"Gria2"};
        //List<TranscriptRecord> transcripts = model.select(gria2);
        
        /* 
        
        editing.csv event list to quantify description file
        name,gene,chromosome,start,end,strand,pos,ref,alt
        RG,Gria2,3,80681450,80802835,-,80692286,T,C
        QR,Gria2,3,80681450,80802835,-,80706912,T,C
        
        */
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        
        try {
            int nb=0;
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(SNP));
            String line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split(",");
                
                // header line
                if("name".equals(tok[0])){
                    line = fichier.readLine(); 
                    tok = line.split(",");
                }
                
                log.info(new Object[]{"\tprocessing...\t\t" + line});
                
                nb++;
                
                String name = tok[0];
                String gene = tok[1];
                String chromosome = tok[2];
                int start = Integer.parseInt(tok[3]);
                int end = Integer.parseInt(tok[4]);
                boolean strand = ("-".equals(tok[5]))?true:false;
                int position = Integer.parseInt(tok[6]);
                String ref = tok[7];
                String alt = tok[8];
                
                // query the whole gene molecules
                SAMRecordIterator iter = samReader.query(chromosome, start, end, false);
                
         	while(iter.hasNext()){
                    SAMRecord r = iter.next();
                    
                    if(strand == r.getReadNegativeStrandFlag()){
                        pl.record(r);

                        // Molecule SAMrecord
                        LongreadRecord lrr = LongreadRecord.fromSAMRecord(r, false);
                        if(lrr != null){
                            Longread lr = new Longread(r.getReadName());
                            lr.addRecord(lrr);
                            String readString = r.getReadString();
                            int posInt = r.getReadPositionAtReferencePosition(position);

                            if(posInt > 0 && readString.length() > posInt){
                                String nuc = readString.substring(posInt - 1, posInt);
                                
                                Molecule molecule = new Molecule(lrr.getBarcode(), lrr.getUmi());
                                molecule.addLongread(lr);
                                molecule.setGeneId(gene);
                                molecule.setTranscriptId(name + "-" + nuc);
                                matrix.addMolecule(molecule);
                            }
                        }
                    }
                }
                
                iter.close();
                line = fichier.readLine();
            }
            samReader.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        File MATRIX = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_matrix.txt");
        File METRICS = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_metrics.txt");
        matrix.writeIsoformMatrix(MATRIX, METRICS);

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
        System.exit(new SNPMatrix().instanceMain(paramArrayOfString));
    }
}
