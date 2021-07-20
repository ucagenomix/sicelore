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
import java.util.Arrays;
import java.util.Collections;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.*;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SNP/editing events detection/quantification (cellBC/UMI count matrix).", oneLineSummary = "SNP/editing event detection/quantification (cellBC/UMI count matrix).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class SNPMatrix extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input molecule SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CSV", doc = "The cell barcodes .csv file")
    public File CSV;
    @Argument(shortName = "SNP", doc = "The SNP/editing events .csv file \n#-----\nchromosome,position,strand,name\n3,80692286,-,Gria2_RG\n3,80706912,-,Gria2_QR\n3,80692286|80706912,-,Gria2_RGQR\n#-----")
    public File SNP;
    //@Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    //public File REFFLAT;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "PREFIX", doc = "The output file prefix (default=snp)")
    public String PREFIX = "snp";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC; illumina=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8, illumina=UB)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "RNTAG", doc = "Read number tag (default=RN)", optional=true)
    public String RNTAG = "RN";
    @Argument(shortName = "MINRN", doc = "Minimum read number to keep the molecule for SNP calling (default=0, means all)")
    public int MINRN = 0;
    @Argument(shortName = "MINQV", doc = "Minimum QV score at position to keep the molecule for SNP calling (default=0, means all)")
    public int MINQV = 0;

    public CellList cellList;
    
    public SNPMatrix() {
        log = Log.getInstance(SNPMatrix.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);
        IOUtil.assertFileIsReadable(SNP);
        //IOUtil.assertFileIsReadable(REFFLAT);
        
        int total_hits=0;
        int total_lowrn=0;
        int total_lowqv=0;
        
        LongreadRecord lrr = new LongreadRecord();
	lrr.setStaticParams(CELLTAG,UMITAG,"GN","TE","UE","PE","US",150, RNTAG);

        this.cellList = new CellList(CSV); 
        log.info(new Object[]{"Cells detected\t\t[" + this.cellList.size() + "]"});
        Matrix matrix = new Matrix(this.cellList);
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
        htsjdk.samtools.SAMSequenceDictionary dictionnary = samFileHeader.getSequenceDictionary();
        
        //editing.csv event list to quantify description file
        //chromosome,position,strand,name
        //11,75300577,-,Rpa1
        
        try {
            int nb=0;
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(SNP));
            String line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split(",");
                
                //System.out.println(tok[0] + "\t" + dictionnary.getSequence(tok[0]));
                
                if(dictionnary.getSequence(tok[0]) != null){
                    String chromosome = tok[0];
                    boolean strand = ("-".equals(tok[2]))?true:false;
                    String gene = tok[3];
                    int total=0;
                    int lowQV=0;
                    int lowRN=0;
                    nb++;
                    
                    String[] pos = { tok[1] };
                    if(java.util.regex.Pattern.matches(".*\\|.*", tok[1]))
                        pos = tok[1].split("\\|");
                        
                    int[] arr = stringToIntArray(pos);
                    SAMRecordIterator iter = samReader.query(chromosome, arr[0], arr[arr.length - 1], false);
                    while (iter.hasNext()) {
                        SAMRecord r = iter.next();

                        if (strand == r.getReadNegativeStrandFlag()) {
                            pl.record(r);
                            lrr = LongreadRecord.fromSAMRecord(r, false);
                            
                            //System.out.println(lrr + "\t" + lrr.getRn() + "\t" + r.getReadString());
                            
                            if (lrr != null){
                                Longread lr = new Longread(r.getReadName());
                                lr.addRecord(lrr);
                                String readString = r.getReadString();
                                String qvString = r.getBaseQualityString();

                                Integer[] posInt = new Integer[arr.length];
                                for (int i = 0; i < arr.length; i++) {
                                    posInt[i] = r.getReadPositionAtReferencePosition(arr[i]);
                                }

                                int min = Collections.min(Arrays.asList(posInt));
                                int max = Collections.max(Arrays.asList(posInt));

                                if(min > 0 && readString.length() > max) {
                                    String[] nuc = new String[posInt.length];
                                    int[] qv = new int[posInt.length];
                                    int min_qv = 100;

                                    for(int i = 0; i < arr.length; i++) {
                                        nuc[i] = readString.substring(posInt[i] - 1, posInt[i]);
                                        //System.out.println(r.getReadName()+"\t"+lrr.getRn()+"\t"+posInt[i]+"\t"+readString.length()+"\t"+qvString.length());
                                        qv[i] = (int)(qvString.substring(posInt[i] - 1, posInt[i]).charAt(0)) - 33;
                                        if(qv[i]<min_qv)
                                            min_qv = qv[i];
                                    }

                                    if (r.getReadNegativeStrandFlag()) {
                                        for (int i = 0; i < nuc.length; i++) {
                                            nuc[i] = complementBase(nuc[i].charAt(0));
                                        }
                                    }

                                    if(lrr.getRn() >= MINRN) {
                                        // keep high QV molecules only
                                        if(min_qv >= MINQV){
                                            total++;
                                            Molecule molecule = new Molecule(lrr.getBarcode(), lrr.getUmi(), lrr.getRn());
                                            molecule.addLongread(lr);
                                            molecule.setGeneId(gene);
                                            molecule.setTranscriptId(chromosome + ":" + String.join("|", pos) + ".." + String.join("", nuc));
                                            molecule.setSnpPhredScore(StringUtils.join(ArrayUtils.toObject(qv), ","));
                                            matrix.addMolecule(molecule);
                                        }
                                        else
                                            lowQV++;
                                    }
                                    else
                                        lowRN++;
                                }
                            }
                        }
                    }
                    iter.close();
                    log.info(new Object[]{"processing...\t\t" + line + "\t" + total + " hits, " + lowRN + " lowRN, " + lowQV + " lowQV"});
                    total_hits+=total;
                    total_lowrn+=lowRN;
                    total_lowqv+=lowQV;
                }
                line = fichier.readLine();
            }
            samReader.close();
        } catch (Exception e) { e.printStackTrace(); }

        if(matrix.getMatrice().size() > 0){
            File MATRIX = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_snpmatrix.txt");
            File METRICS = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_snpmetrics.txt");
            File MOLINFOS = new File(OUTPUT.getAbsolutePath() + "/" + PREFIX + "_snpmolinfos.txt");
            matrix.writeIsoformMatrix(MATRIX, METRICS, MOLINFOS, null);
        }
        else
            log.info(new Object[]{"end of processing...\tnothing has been detected, check your input parameters (if Illumina set CELLBC=CB UMITAG=UB), no output files generated\t"});
        
        
        log.info(new Object[]{"STATISTICS...\t\thits=" + total_hits + ", lowRN=" + total_lowrn + ", lowQV= " + total_lowqv});
        
        return 0;
    }
    
    public int[] stringToIntArray(String[] s)
    {
        int size = s.length;
        int [] arr = new int [size];
        for(int i=0; i<size; i++) 
            arr[i] = Integer.parseInt(s[i]);
        
        Arrays.sort(arr);
        return arr;
    }
    
    protected static String complementBase(char base)
    {
        String cc = "";
        if (base == 'A') cc = "T";
        if (base == 'C') cc = "G";
        if (base == 'G') cc = "C";
        if (base == 'T') cc = "A";
	
        return cc;
    }
    
    public String[] complementBaseArray(String[] nuc)
    {
        int size = nuc.length;
        String[] comp = new String[size];
        for(int i=0; i<size; i++) 
            comp[i] = complementBase(nuc[i].charAt(0));
        
        return comp;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new SNPMatrix().instanceMain(paramArrayOfString));
    }
}
