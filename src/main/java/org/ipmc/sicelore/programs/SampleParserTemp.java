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
import java.util.HashSet;
import java.util.List;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SampleParserTemp", oneLineSummary = "SampleParserTemp", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SampleParserTemp extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "CSV", doc = "The samples to process")
    public File CSV;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "PREFIX", doc = "The output file prefix (default=out)")
    public String PREFIX = "out";
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=CB)", optional=true)
    public String CELLTAG = "CB";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=UB)", optional=true)
    public String UMITAG = "UB";
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=GN)", optional=true)
    public String GENETAG = "GN";
    
    public SampleParserTemp() {
        log = Log.getInstance(SampleParserTemp.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(CSV);

        //hts_sample_ref,mix_ref,path
        //170619_10x_vassaux_3_02,vassaux3,/data/public/hts/5148/5776_vassaux3_infected_possorted_genome.bam
        
        try {
            int nb=0;
            BufferedReader fichier = new BufferedReader(new java.io.FileReader(CSV));
            String line = fichier.readLine();
            line = fichier.readLine();
            while(line != null && !"".equals(line)){
                String[] tok = line.split(",");
                
                //log.info(new Object[]{"processing...\t\t" + line});
                
                IOUtil.assertFileIsReadable(new File(tok[2]));
                File INPUT = new File(tok[2]);
                SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
                htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
                htsjdk.samtools.SAMSequenceDictionary dictionnary = samFileHeader.getSequenceDictionary() ;
                List<SAMProgramRecord> progs = samFileHeader.getProgramRecords();
                
                String build = "undef";
                for (SAMProgramRecord p : progs){
                    String cmd = p.getCommandLine();
                    if(cmd != null && !"".equals(cmd)){
                        if(java.util.regex.Pattern.matches(".*GRCh38.*", cmd))
                            build = "hg38";
                        if(java.util.regex.Pattern.matches(".*hg38.*", cmd))
                            build = "hg38";
                        else if(java.util.regex.Pattern.matches(".*hg19.*", cmd))
                            build = "hg19";
                        else if(java.util.regex.Pattern.matches(".*mm10.*", cmd))
                            build = "mm10";
                        else if(java.util.regex.Pattern.matches(".*Sscrofa11-2.*", cmd))
                            build = "Sscrofa11";
                        else if(java.util.regex.Pattern.matches(".*canFam3.*", cmd))
                            build = "canFam3";
                        
                    }
                }
                
                if(!"undef".equals(build)){
                    String gene = "ACE2";
                    String chromosome = "1";
                    int pos1 = 1;
                    int pos2 = 20;
                    
                    HashSet<String> mol = new HashSet();
                    int reads=0;
                    
                    if("hg19".equals(build)){ chromosome="X"; pos1=15579155; pos2=15580136; gene="ACE2"; }
                    if("hg38".equals(build)){ chromosome="X"; pos1=15561032; pos2=15562013; gene="ACE2"; }
                    if("mm10".equals(build)){ chromosome="X"; pos1=164187491 ; pos2=164188420 ; gene="Ace2"; }
                    if("Sscrofa11".equals(build)){ chromosome="X"; pos1=12099852 ; pos2=12100175 ; gene="ACE2"; }
                    if("canFam3".equals(build)){ chromosome="X"; pos1=11838383; pos2=11838668 ; gene="ACE2"; }
                    
                    SAMRecordIterator iter = samReader.query(chromosome, pos1, pos2, false);
                    while(iter.hasNext()){
                        SAMRecord r = iter.next();
                        
                        String CB = (String)r.getAttribute(CELLTAG);
                        String UB = (String)r.getAttribute(UMITAG);
                        String GN = (String)r.getAttribute(GENETAG);
                        int txStart = r.getAlignmentStart();
                        int txEnd = r.getAlignmentEnd();
                        
                        if(gene.equals(GN)){
                            reads++;
                            mol.add(CB+"-"+UB);
                        }
                    }
                    iter.close();
                    
                    log.info(new Object[]{tok[0] + "\t" + tok[1] + "\t" + tok[2] + "\t" + build + "\t" + reads + "\t" + mol.size()});
                }
                else
                    log.info(new Object[]{line + " not found"});
                    
                samReader.close();
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); }

        return 0;
    }
    
    public static void main(String[] paramArrayOfString) {
        System.exit(new SampleParserTemp().instanceMain(paramArrayOfString));
    }
}
