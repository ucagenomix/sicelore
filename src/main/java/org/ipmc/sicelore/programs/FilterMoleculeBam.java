package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Molecule Bam filtering on UMI depth and isoform level annotation.", oneLineSummary = "Molecule Bam filtering on UMI depth and isoform level annotation.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class FilterMoleculeBam extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file")
    public File OUTPUT;
    @Argument(shortName = "CHECKRN", doc = "Wether or not to check RN value (Default=true)")
    public boolean CHECKRN = true;
    @Argument(shortName = "CHECKUNDEF", doc = "Wether or not to check and filter molecules isoform level undefined (ISOTAG=\"undef\") (Default=true)")
    public boolean CHECKUNDEF = true;
    @Argument(shortName = "CHECKGEIG", doc = "Wether or not to check GE/IG compatibility, if true should filter misasigned chimeric SAM records (Default=true)")
    public boolean CHECKGEIG = true;
    @Argument(shortName = "RNTAG", doc = "The read number tag (default=RN)")
    public String RNTAG = "RN";
    @Argument(shortName = "VALUE", doc = "The value cutoff (default=1)")
    public int VALUE=1; 
    @Argument(shortName = "ISOTAG", doc = "The isoform tag (default=IT)")
    public String ISOTAG = "IT";
    @Argument(shortName = "ISOGENETAG", doc = "The isoform gene name tag (default=IG)")
    public String ISOGENETAG = "IG";

    public FilterMoleculeBam() {
        log = Log.getInstance(FilterMoleculeBam.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            for (SAMRecord r : localSamReader) {
                pl.record(r);
                String IG = (String)r.getAttribute(ISOGENETAG);
                String IT = (String)r.getAttribute(ISOTAG);
                
                boolean keepit = true;
                
                if(CHECKGEIG){
                    if(r.getAttribute("GE") != null){
                        String[] GEs = ((String)r.getAttribute("GE")).split(",");
                        Set<String> set = new HashSet<>(Arrays.asList(GEs));
                        if(!set.contains(IG))
                            keepit = false;
                    }
                    else
                        keepit = false; 
                }
                if(CHECKRN){
                    int rn = (Integer) r.getAttribute(RNTAG);
                    if(rn < VALUE)
                        keepit = false;
                }
                if(CHECKUNDEF && "undef".equals(IT))
                    keepit = false;
                    
                if(keepit)
                    localSAMFileWriter.addAlignment(r);
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception localException) { localException.printStackTrace(); }

        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new FilterMoleculeBam().instanceMain(paramArrayOfString));
    }
}
