package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import htsjdk.samtools.SAMFileHeader.SortOrder; 
import java.util.HashSet;
import java.util.Set;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.ipmc.sicelore.utils.LongreadParser;
import org.ipmc.sicelore.utils.MoleculeDataset;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Count unique cellBC/UMI in bam file.", oneLineSummary = "Count unique cellBC/UMI in bam file.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class MoleculeCounter extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The reference input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";

    public MoleculeCounter() {
        log = Log.getInstance(AddReadsToMolecules.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);

        HashSet<String> uniqUMIs = new HashSet<String>();
        HashSet<String> uniqReads = new HashSet<String>();

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String name = r.getReadName();
                String BC = (String)r.getAttribute(CELLTAG);
                String U8 = (String)r.getAttribute(UMITAG);
                
                uniqReads.add(name);
                
                if(BC != null && U8 != null){
                    uniqUMIs.add(BC+":"+U8);
                    
                }
            }
            localSamReader.close();
            
                        
        } catch (Exception localException) {
            localException.printStackTrace();
        }
        
        log.info(new Object[]{"\treads\t\t" + uniqReads.size()});
        log.info(new Object[]{"\tmolecules\t" + uniqUMIs.size()});

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new AddReadsToMolecules().instanceMain(paramArrayOfString));
    }
}