package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import gnu.trove.THashSet;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.BufferedReader;
import java.io.File;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "export records from reads ID list received.", oneLineSummary = "export records from reads ID list received.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SplitBam extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file to analyze")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output directory")
    public File OUTPUT;
    @Argument(shortName = "ID", doc = "The reads ids (one per line")
    public File ID;

    public SplitBam() {
        log = Log.getInstance(SplitBam.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        String line = null;

        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(ID);

        THashSet<String> ids = new THashSet<String>();

        SamReader input = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader header = input.getFileHeader();
        header.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);
        SAMFileWriter yes = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, new File(OUTPUT.getAbsolutePath() + "/yes.bam"));
        SAMFileWriter no = new SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, new File(OUTPUT.getAbsolutePath() + "/no.bam"));
        
        try {
            BufferedReader br = new BufferedReader(new java.io.FileReader(ID));
            line = br.readLine();
            while(line != null) {
                
                if("".equals(line)){}
                else{
                    line = line.replaceAll("@","");
                    ids.add(line);
                }
                line = br.readLine();
            }
            br.close();

            for (Iterator localIterator = input.iterator(); localIterator.hasNext();)
            {
                SAMRecord rec = (SAMRecord) localIterator.next();
                pl.record((SAMRecord) rec);
                String name = (String)rec.getReadName();
                String[] tab = name.split("_");
                
                if(ids.contains(tab[0]))
                    yes.addAlignment(rec);
                else
                    no.addAlignment(rec);
            }

            input.close();
            yes.close();
            no.close();
            
        } catch (Exception e) { e.printStackTrace(); }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new SplitBam().instanceMain(paramArrayOfString));
    }
}
