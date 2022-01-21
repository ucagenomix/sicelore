package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*; 
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import java.util.HashMap;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Isobam production.", oneLineSummary = "Isobam production.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class Isobam extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "MOLINFOS", doc = "The _molinfos.txt file from IsoformMatrix")
    public File MOLINFOS;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";
    @Argument(shortName = "GENETAG", doc = "The gene tag (default=IG)")
    public String GENETAG = "IG";
    @Argument(shortName = "ISOTAG", doc = "The isoform tag (default=IT)")
    public String ISOTAG = "IT";
    @Argument(shortName = "UNDEF", doc = "Wether or not to keep molecules not define at the isoform level (ISOTAG=\"undef\") (Default=true)")
    public boolean UNDEF = true;
    
    private ProgressLogger pl;
    private final Log log;

    public Isobam() {
        log = Log.getInstance(Isobam.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(MOLINFOS);
        process();
        return 0;
    }

    protected void process()
    {
        HashMap<String, String> geneIds = new HashMap<String, String>();
        HashMap<String, String> transcriptIds = new HashMap<String, String>();
        //cellBC  UMI     nbReads nbSupportingReads       mappingPctId    snpPhredScore   geneId  transcriptId
        //ATAGGCTGTGTAGCAG        GCTCGAAGCCCC    2       0       0.9594          Wdr38   undef

        try {
            BufferedReader fichier = new BufferedReader(new FileReader(MOLINFOS));
            String line = fichier.readLine();
            while(line != null) {
                String[] tmp = line.split("\t");
                if(!UNDEF && !"undef".equals(tmp[7])){
                    geneIds.put(tmp[0]+":"+tmp[1],tmp[6]);
                    transcriptIds.put(tmp[0]+":"+tmp[1],tmp[7]);
                }
                else if(UNDEF){
                    geneIds.put(tmp[0]+":"+tmp[1],tmp[6]);
                    transcriptIds.put(tmp[0]+":"+tmp[1],tmp[7]);
                }
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        SamReader samReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = samReader.getFileHeader();
        samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, OUTPUT);
        try{
            for(SAMRecord r : samReader){
                pl.record(r);
                String isokey = (String)r.getAttribute(CELLTAG)+":"+(String)r.getAttribute(UMITAG);
                
                if(geneIds.containsKey(isokey)){
                    r.setAttribute(GENETAG, (String)geneIds.get(isokey));
                    r.setAttribute(ISOTAG, (String)transcriptIds.get(isokey));
                    samFileWriter.addAlignment(r);
                }
            }
            samReader.close();
            samFileWriter.close();
        }catch (Exception e) { e.printStackTrace(); }
    }

    public static void main(String[] args) {
        System.exit(new Isobam().instanceMain(args));
    }
}
