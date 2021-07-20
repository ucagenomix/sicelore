package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.BufferedReader;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.util.*;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Isoform junctions validation tool", oneLineSummary = "Isoform junctions validation tool", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class JunctionValidator extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The short read SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "REFFLAT", doc = "The gencode refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "TXT", doc = "SQANTI2 junctions.txt file")
    public File TXT;
    @Argument(shortName = "O", doc = "The output file name")
    public File OUTPUT;
    
    HashMap<String,Integer> statistics;
            
    public JunctionValidator() {
        log = Log.getInstance(JunctionValidator.class);
        pl = new ProgressLogger(log);
        
        statistics = new HashMap<String,Integer>();
        statistics.put("novel_canonical_gencode_validated",0);
        statistics.put("novel_non_canonical_gencode_validated",0);
        statistics.put("known_canonical_gencode_validated",0);
        statistics.put("known_non_canonical_gencode_validated",0);
        statistics.put("novel_canonical_gencode_unvalidated",0);
        statistics.put("novel_non_canonical_gencode_unvalidated",0);
        statistics.put("known_canonical_gencode_unvalidated",0);
        statistics.put("known_non_canonical_gencode_unvalidated",0);
        
        statistics.put("novel_canonical_notingencode_validated",0);
        statistics.put("novel_non_canonical_notingencode_validated",0);
        statistics.put("known_canonical_notingencode_validated",0);
        statistics.put("known_non_canonical_notingencode_validated",0);
        statistics.put("novel_canonical_notingencode_unvalidated",0);
        statistics.put("novel_non_canonical_notingencode_unvalidated",0);
        statistics.put("known_canonical_notingencode_unvalidated",0);
        statistics.put("known_non_canonical_notingencode_unvalidated",0);
    }

    protected int doWork()
    {
        int nb=0;
        HashMap<String,String> already = new HashMap<String,String>();
        DataOutputStream os = null;
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(TXT);
        
        /* 
     
        isoform chrom   strand  junction_number genomic_start_coord     genomic_end_coord       transcript_coord        junction_category       start_site_category     end_site_category       diff_
        to_Ref_start_site       diff_to_Ref_end_site    bite_junction   splice_site     canonical       RTS_junction    indel_near_junct        phyloP_start    phyloP_end      sample_with_cov total_coverage
        PB.135.19       1       -       junction_1      4774517 4777524 ?????   known   known   known   0       0       TRUE    GTAG    canonical       FALSE   FALSE   NA      NA      NA      NA
        PB.135.19       1       -       junction_2      4777649 4782567 ?????   known   known   known   0       0       TRUE    GTAG    canonical       FALSE   FALSE   NA      NA      NA      NA
        PB.135.19       1       -       junction_3      4782734 4783950 ?????   known   known   known   0       0       TRUE    GTAG    canonical
        
        */
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        HashSet<String> splices_sites = new HashSet<String>();
        
        try {
            log.info(new Object[]{"refflat loading start..."});
            
            FileInputStream in = new FileInputStream(REFFLAT);
            BufferedReader reader = new BufferedReader(new InputStreamReader(in));
            String line = null;
            while ((line = reader.readLine()) != null){
                String[] fields = line.split("\t");
                fields[9] = StringUtils.stripEnd(fields[9], ",");
                String[] vals = fields[9].split(",");
                fields[2] = fields[2].replaceAll("chr","");
                
                for (int i = 1; i < vals.length; i++){
                    splices_sites.add(fields[2]+":"+(Integer.valueOf(vals[i])));
                    //System.out.println("push:"+fields[2]+":"+(Integer.valueOf(vals[i])));
                }
                
                fields[10] = StringUtils.stripEnd(fields[10], ",");
                String[] vals2 = fields[10].split(",");
                for (int i = 0; i < vals2.length -1; i++){
                    splices_sites.add(fields[2]+":"+(Integer.valueOf(vals2[i])+1));
                    //System.out.println("push:"+fields[2]+":"+(Integer.valueOf(vals2[i])+1));
                }
            }
            reader.close();
            
            log.info(new Object[]{"refflat loading end"});
            log.info(new Object[]{"Splice sites loaded\t\t[" + splices_sites.size() + "]"});

            os = new DataOutputStream(new FileOutputStream(OUTPUT));
            os.writeBytes("isoform\tjunction\tchrom\tstrand\tstart\tend\tjunction_category\tcanonical\tboth_splice_sites_in_refflat\tsupporting_reads\n");
            
            int e_prev=0;
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(TXT));
            line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split("\t");
                nb++;
                
                if("isoform".equals(tok[0])){
                    line = fichier.readLine(); 
                    tok = line.split("\t");
                }
                
                String isoform = tok[0];
                String chrom = tok[1];
                boolean strand = ("-".equals(tok[2]))?true:false;
                String junction_number = tok[3];
                int genomic_start_coord = Integer.parseInt(tok[4]);
                int genomic_end_coord = Integer.parseInt(tok[5]);
                String junction_category = tok[7];
                String canonical = tok[14];
                
                boolean both_is_gencode = false;
                if(splices_sites.contains(tok[1]+":"+tok[4])){
                    
                    if(splices_sites.contains(tok[1]+":"+tok[5]))
                        both_is_gencode = true;
                }
                
                //System.out.println("is_in:"+tok[1]+":"+tok[4]+","+tok[1]+":"+tok[5]+" --> "+both_is_gencode);
                
                
                String isokey = chrom + ":" + genomic_start_coord + "-" + genomic_end_coord;
                if(already.containsKey(isokey)){ // already process and printed
                    os.writeBytes(isoform+"\t"+already.get(isokey));
                } 
                else{
                    int supporting_reads = 0;
                    // query short reads spanning the junction
                    int nb_reads=0;
                    SAMRecordIterator iter = samReader.query(chrom, genomic_start_coord, genomic_end_coord, false);
                    while(iter.hasNext()){
                        SAMRecord r = iter.next();
                        String name = r.getReadName();
                        
                        //System.out.println(name);
                        
                        //if(strand == r.getReadNegativeStrandFlag()){
                            //pl.record(r);
                            
                            List<int[]> junctions = new ArrayList<int[]>();
                            List<AlignmentBlock> blocks = r.getAlignmentBlocks();
                            for (int i = 0; i < blocks.size(); i++) {
                                AlignmentBlock currBlock = blocks.get(i);
                                int s = currBlock.getReferenceStart();
                                int e = s + currBlock.getLength();
                                
                                if(i>0)
                                    junctions.add(new int[]{e_prev, s-1});
                                
                                e_prev = e;
                            }
                            
                            if(isIn(new int[]{genomic_start_coord, genomic_end_coord}, junctions))
                                supporting_reads++;
                        //}
                        
                        nb_reads++;
                    }
                    iter.close();
                    
                    String info = isoform+"\t"+junction_number+"\t"+chrom+"\t"+tok[2]+"\t"+genomic_start_coord+"\t"+genomic_end_coord+"\t"+junction_category+"\t"+canonical+"\t"+both_is_gencode+"\t"+supporting_reads+"\n";
                    os.writeBytes(info);
                    already.put(isokey, junction_number+"\t"+chrom+"\t"+tok[2]+"\t"+genomic_start_coord+"\t"+genomic_end_coord+"\t"+junction_category+"\t"+canonical+"\t"+both_is_gencode+"\t"+supporting_reads+"\n");
                    
                    if(both_is_gencode){
                        if(supporting_reads > 0)
                            statistics.put(new String(junction_category+"_"+canonical+"_gencode_validated"), (Integer)statistics.get(new String(junction_category+"_"+canonical+"_gencode_validated")) + 1);
                        else
                            statistics.put(new String(junction_category+"_"+canonical+"_gencode_unvalidated"), (Integer)statistics.get(new String(junction_category+"_"+canonical+"_gencode_unvalidated")) + 1);                        
                    }
                    else{
                        if(supporting_reads > 0)
                            statistics.put(new String(junction_category+"_"+canonical+"_notingencode_validated"), (Integer)statistics.get(new String(junction_category+"_"+canonical+"_notingencode_validated")) + 1);
                        else
                            statistics.put(new String(junction_category+"_"+canonical+"_notingencode_unvalidated"), (Integer)statistics.get(new String(junction_category+"_"+canonical+"_notingencode_unvalidated")) + 1);                        
                    }
                }
                
                if(nb%1000 == 0)
                    displayStatistics(nb);

                line = fichier.readLine();
            }
            samReader.close();
            fichier.close();
            os.close();
            
        } catch (Exception e) { e.printStackTrace(); } 
        finally { try { samReader.close();  os.close(); } catch (Exception e) { System.err.println("can not close stream"); } }
        
        //                      validated |not validated
        // known_canonical:     [67191    |39] 
        // known_non_canonical: [94       |6]
        // novel_canonical:     [2852     |764] 
        // novel_non_canonical: [7        |1227]
        displayStatistics(nb);

        return 0;
    }
    
    public void displayStatistics(int nb)
    {
        String str = nb + " junctions processed\t";
        str += "known_canonical[i:"+statistics.get("known_canonical_gencode_validated")+"/"+statistics.get("known_canonical_gencode_unvalidated")+ ",o:"+statistics.get("known_canonical_notingencode_validated")+"/"+statistics.get("known_canonical_notingencode_unvalidated")+ "]";
        str += ", known_non_canonical[i:"+statistics.get("known_non_canonical_gencode_validated")+"|"+statistics.get("known_non_canonical_gencode_unvalidated")+ ",o:"+statistics.get("known_non_canonical_notingencode_validated")+"/"+statistics.get("known_non_canonical_notingencode_unvalidated")+ "]";
        str += ", novel_canonical[i:"+statistics.get("novel_canonical_gencode_validated")+"|"+statistics.get("novel_canonical_gencode_unvalidated")+ ",o:"+statistics.get("novel_canonical_notingencode_validated")+"/"+statistics.get("novel_canonical_notingencode_unvalidated")+ "]";
        str += ", novel_non_canonical[i:"+statistics.get("novel_non_canonical_gencode_validated")+"|"+statistics.get("novel_non_canonical_gencode_unvalidated")+ ",o:"+statistics.get("novel_non_canonical_notingencode_validated")+"/"+statistics.get("novel_non_canonical_notingencode_unvalidated")+ "]";
        
        log.info(new Object[]{str});
    }
    
    public List<int[]> junctionsFromExons(List<int[]> exons) {
        ArrayList lst = new ArrayList();

        for (int i = 1; i < exons.size(); i++) {
            int j = ((int[]) exons.get(i - 1))[1];
            int k = ((int[]) exons.get(i))[0];
            lst.add(new int[]{j, k});
        }

        return lst;
    }
    
    public boolean isIn(int[] junc, List<int[]> list) {
        boolean bool = false;
        for (int[] a : list) {
            
            //System.out.println("comparing " + a[0] + "==" + junc[0] +" and " + a[1] +"==" +junc[1]);
            
            if((a[0] == junc[0]) && (a[1] == junc[1]))
                bool = true;
        }
        return bool;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new JunctionValidator().instanceMain(paramArrayOfString));
    }
}
