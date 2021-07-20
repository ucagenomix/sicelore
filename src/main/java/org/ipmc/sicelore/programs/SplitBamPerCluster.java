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
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Bam file cluster-by-cluster splitter.", oneLineSummary = "Bam file cluster-by-cluster splitter.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SplitBamPerCluster extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "CSV", doc = "The \"cellBC,clusterName\" cluster file (.csv)")
    public File CSV;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";

    public SplitBamPerCluster() {
        log = Log.getInstance(SplitBamPerCluster.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        String str1 = null;

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(CSV);

        HashMap localHashMap1 = new HashMap();
        HashMap localHashMap2 = new HashMap();

        SamReader localSamReader = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        try {
            BufferedReader br = new BufferedReader(new java.io.FileReader(CSV));
            str1 = br.readLine();
            str1 = br.readLine();
            while(str1 != null) {
                
                if("".equals(str1)){}
                else{
                    str1 = str1.replaceAll("\"","");
                    str1 = str1.replaceAll(" ","_");
                    String[] line = str1.split(",");
                    line[0] = line[0].replaceAll("-1","");
                    
                    localHashMap1.put(line[0], line[1]);
                    localHashMap1.put(line[0] + "-1", line[1]);
                    
                    if (!localHashMap2.containsKey(line[1]))
                        localHashMap2.put(line[1], new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, new File(OUTPUT.getAbsolutePath() + "/" + line[1] + ".bam")));
                }
                str1 = br.readLine();
            }
            br.close();

            for (Iterator localIterator = localSamReader.iterator(); localIterator.hasNext();) {
                SAMRecord localObject = (SAMRecord) localIterator.next();
                pl.record((SAMRecord) localObject);
                String str2 = (String) ((SAMRecord) localObject).getAttribute(CELLTAG);
                String str3 = (String) localHashMap1.get(str2);

                if (str3 != null) {
                    ((SAMFileWriter) localHashMap2.get((String) localHashMap1.get(str2))).addAlignment((SAMRecord) localObject);
                }
                //else{
                //    System.out.println(str2 + " not found");
                //}
            }

            localSamReader.close();
            Object localObject = localHashMap2.keySet();
            Iterator localIterator = ((Set) localObject).iterator();
            while (localIterator.hasNext()) {
                ((SAMFileWriter) localHashMap2.get((String) localIterator.next())).close();
            }
        } catch (Exception e) { e.printStackTrace(); }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new SplitBamPerCluster().instanceMain(paramArrayOfString));
    }
}
