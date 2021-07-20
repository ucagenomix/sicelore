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

@CommandLineProgramProperties(summary = "Add reads (targeted exp.) for molecules in input bam (standard exp.).", oneLineSummary = "Add reads (targeted exp.) for molecules in input bam (standard exp.).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class AddReadsToMolecules extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The reference input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "T", doc = "The targeted input SAM or BAM file")
    public File TARGETED;
    @Argument(shortName = "O", doc = "The output SAM or BAM file with tags")
    public File OUTPUT;
    @Argument(shortName = "CELLTAG", doc = "Cell tag (default=BC)", optional=true)
    public String CELLTAG = "BC";
    @Argument(shortName = "UMITAG", doc = "UMI tag (default=U8)", optional=true)
    public String UMITAG = "U8";

    public AddReadsToMolecules() {
        log = Log.getInstance(AddReadsToMolecules.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(TARGETED);
        IOUtil.assertFileIsWritable(OUTPUT);

        HashSet<String> mol2keep = new HashSet<String>();
        int rec=0;
        int addrec=0;

        LongreadParser bam = new LongreadParser(INPUT, false, false, true, true);
        MoleculeDataset dataset = new MoleculeDataset(bam);
        Set cles = dataset.getMapMolecules().keySet();
        
        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        localSAMFileHeader.setSortOrder(SortOrder.unsorted);
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String BC = (String)r.getAttribute(CELLTAG);
                String U8 = (String)r.getAttribute(UMITAG);
                if(cles.contains(BC+":"+U8)){
                    rec++;
                    localSAMFileWriter.addAlignment(r);
                    
                }
            }
            localSamReader.close();
            
            log.info(new Object[]{"\tinput bam file\trecords\t[" + rec + "]"});
            log.info(new Object[]{"\tinput bam file\tmolecules\t[" + cles.size() + "]"});
            
            SamReader localSamReader2 = SamReaderFactory.makeDefault().open(TARGETED);
            
            for (SAMRecord r : localSamReader2) {
                pl.record(r);
                String BC = (String)r.getAttribute(CELLTAG);
                String U8 = (String)r.getAttribute(UMITAG);
                if(cles.contains(BC+":"+U8)){
                    addrec++;
                    localSAMFileWriter.addAlignment(r);
                    mol2keep.add(BC+":"+U8);
                }
            }
            
            log.info(new Object[]{"\toutput bam file\tadded records\t[" + addrec + "]"});
            log.info(new Object[]{"\toutput bam file\tfor molecules\t[" + mol2keep.size() + "]"});
            
            localSamReader2.close();
            localSAMFileWriter.close();
        } catch (Exception localException) {
            localException.printStackTrace();
        }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new AddReadsToMolecules().instanceMain(paramArrayOfString));
    }
}