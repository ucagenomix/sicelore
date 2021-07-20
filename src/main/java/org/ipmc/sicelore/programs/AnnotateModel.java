package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */ 
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.tribble.annotation.Strand;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.BEDParser;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Annotate model file", oneLineSummary = "Annotate model file", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class AnnotateModel extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "MODEL", doc = "Sicelore Model file (TXT file of CollapseModel pipeline")
    public File TXT;
    @Argument(shortName = "I", doc = "The short read SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CAGE", doc = "CAGE peaks file (.bed)")
    public File CAGE;
    
    // ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M24/gencode.vM24.polyAs.gtf.gz
    // more +6 gencode.vM24.polyAs.gtf  | grep polyA_signal | awk '{ print $1 "\t" $4 "\t" $5 "\t" $1 ":" $4 ".." $5 "," $7 "\t1\t" $7 "\t" $4 "\t" $4+1 }' > gencode.vM24.polyAs.bed
    // sed -i -e "s/chr//g" gencode.vM24.polyAs.bed
    
    @Argument(shortName = "POLYA", doc = "POLYA sites file (.bed)")
    public File POLYA;
    @Argument(shortName = "O", doc = "The output file name")
    public File OUTPUT;
    @Argument(shortName = "DELTA", doc = "Allowed base number difference between start/end of exons and read block position (default=0)")
    public int DELTA = 0;
    
    public HashSet<String> allcmp = new HashSet<String>();
    
    BEDParser cage;
    BEDParser polyA;
    
    public AnnotateModel() {
        log = Log.getInstance(JunctionValidator.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        int nb=0;
        int e_prev=0;
        DataOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(TXT);
        HashMap<String, Integer> already = new HashMap<String, Integer>();
        
        cage = new BEDParser(CAGE);
        polyA = new BEDParser(POLYA);
        
        // junction validation
        // CAGE-peak position
        // polyA position
        // canonical junction
        
        /* 
        geneId  TranscriptId            chrom   strand  genomic_start_coord     genomic_end_coord       exons   evidences   category                subcategory                     junctions_list                  exons_list
        Adora2b ENSMUST00000018644.2    11      +       62248983                62266453                2       1           full_splice_match       gencode                         -                               62248984-62249436,62265062-62266453
        Rrp1    Novel.16                10      -       78400385                78413036                13      48          novel_not_in_catalog    at_least_one_novel_splicesite   nss=10:78400981-78401850,...    78400385-78400981,78401850-78401966,78402425-78402444,78402610-78402707,78403505-78403584,78404925-78
        */
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        
        try {
            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(TXT));
            String header = fichier.readLine();
            header += "\tis_validated\tsupport_reads\twithin_cagepeak\tdist_cagepeak\twithin_polya\tdist_polya\n";
            os.writeBytes(header);
            
            String line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split("\t");
                nb++;
                
                boolean is_validated = true;
                String support_info="";
                
                String chrom = tok[2];
                Strand strand = Strand.toStrand(tok[3]);
                int genomic_start_coord = Integer.parseInt(tok[4]);
                int genomic_end_coord = Integer.parseInt(tok[5]);
                String category = tok[8];
                String subcategory = tok[9];
                String junctions_list = tok[10];
                
                // get distance of isoforms to cage-peak
                int dist_cage_peak = cage.getDistanceCage(chrom, strand, (strand == Strand.POSITIVE)?genomic_start_coord:genomic_end_coord);
                boolean within_cage_peak = cage.isWithin(chrom, strand, (strand == Strand.POSITIVE)?genomic_start_coord:genomic_end_coord);
                
                // get distance of isoforms to polyA
                int dist_polyA = polyA.getDistancePolyA(chrom, strand, (strand == Strand.POSITIVE)?genomic_end_coord:genomic_start_coord);
                boolean within_polyA = polyA.isWithin(chrom, strand, (strand == Strand.POSITIVE)?genomic_end_coord:genomic_start_coord);
                
                //System.out.println(strand+"\t"+dist_cage_peak);
                //int dist_polya = polyA.getDistance(chrom, strand, (strand)?genomic_end_coord:genomic_start_coord);
                
                if(!"-".equals(junctions_list)){
                    String[] lst = junctions_list.split(",");
                    for(int index=0; index<lst.length; index++)
                    {
                        String[] tmp = lst[index].split(":");
                        String[] tmp2 = tmp[1].split("-");
                        int donor = Integer.parseInt(tmp2[0]);
                        int acceptor = Integer.parseInt(tmp2[1]);
                        String isokey = chrom+":"+donor+"-"+acceptor;
                        
                        if(! already.containsKey(isokey)){
                            int nb_reads=0;
                            int supporting_reads = 0;
                            SAMRecordIterator iter = samReader.query(chrom, donor, acceptor, false);
                            while(iter.hasNext()){
                                SAMRecord r = iter.next();

                                List<int[]> junctions = new ArrayList<int[]>();
                                List<AlignmentBlock> blocks = r.getAlignmentBlocks();
                                for (int i=0; i<blocks.size(); i++) {
                                    AlignmentBlock currBlock = blocks.get(i);
                                    int s = currBlock.getReferenceStart();
                                    int e = s + currBlock.getLength();
                                    if(i>0){ junctions.add(new int[]{e_prev-1, s}); }
                                    e_prev = e;
                                }
                                //log.info(new Object[]{lst[index] + "\t" + junctions.size()});
                                if(isIn(new int[]{donor, acceptor}, junctions)){
                                    supporting_reads++;
                                    //log.info(new Object[]{r.getReadNegativeStrandFlag()});
                                }

                                nb_reads++;
                            }
                            iter.close();
                            
                            already.put(isokey, supporting_reads);
                        }
                        if(already.get(isokey)== 0)
                            is_validated = false;
                        
                        support_info += already.get(isokey) + ",";
                    }
                }
                if("".equals(support_info))
                    support_info = "-";
                else
                    support_info = support_info.substring(0,support_info.length()-1);
                
                String info = line + "\t" + is_validated + "\t" + support_info + "\t" + within_cage_peak + "\t" + dist_cage_peak + "\t" + within_polyA + "\t" + dist_polyA + "\n";
                os.writeBytes(info);
                
                if(nb%2500 == 0)
                    log.info(new Object[]{nb});
                
                line = fichier.readLine();
            }
            samReader.close();
            fichier.close();
            os.close();
            
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { samReader.close();  os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        /*
        Iterator<String> iterator = allcmp.iterator();
        while (iterator.hasNext()) {
            String s = iterator.next();
            System.out.println(s);
        }
        */
        
        return 0;
    }
    
    public boolean isIn(int[] junc, List<int[]> list) {
        boolean bool = false;
        
        for (int[] a : list) {
            allcmp.add(junc[0] + "-" + junc[1] + "/"+ a[0] +"-" + a[1]);
            //System.out.println("comparing " + a[0] + "==" + junc[0] +" and " + a[1] +"==" +junc[1]);
            //if((a[0] == junc[0]) && (a[1] == junc[1]))
            if ((Math.abs(a[0] - junc[0]) <= DELTA) && (Math.abs(a[1] - junc[1]) <= DELTA))
                bool = true;
        }
        return bool;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new AnnotateModel().instanceMain(paramArrayOfString));
    }
}
