package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "SAMrecord read sequence and QV tagging.", oneLineSummary = "SAMrecord read sequence and QV tagging.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class AddBamReadSequenceTag extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output BAM file")
    public File OUTPUT;
    @Argument(shortName = "FASTQDIR", doc = "The .fastq files directory")
    public File FASTQDIR;
    @Argument(shortName = "SEQTAG", doc = "The sequence tag <default=US, use CS if added using ParseFastq input>")
    public String SEQTAG = "US";
    @Argument(shortName = "QVTAG", doc = "The QV tag <default=UQ>")
    public String QVTAG = "UQ";
    @Argument(shortName = "addQV", doc = "Wheter or not add the QV tag <default=true>")
    public boolean addQV = true;
    
    public AddBamReadSequenceTag() {
        log = Log.getInstance(AddBamReadSequenceTag.class);
        pl = new ProgressLogger(log);
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        //IOUtil.assertFileIsReadable(FASTQDIR);
        IOUtil.assertFileIsWritable(OUTPUT);

        log.info(new Object[]{"loadFastq\tSTART..."});
        FastqLoader localFastqLoader = new FastqLoader(FASTQDIR, addQV);
        
        log.info(new Object[]{"loadFastq\t" + localFastqLoader.getMap().size() + " reads loaded"});

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            log.info(new Object[]{"addTag\t\tSTART..."});
            for (SAMRecord localSAMRecord : localSamReader) {
                pl.record(localSAMRecord);
                String name = localSAMRecord.getReadName();
                
                //log.info(new Object[]{name});
                
                if(localFastqLoader.getMap().containsKey(name)){
                    String seq = new String((byte[]) localFastqLoader.getMap().get(name));
                    localSAMRecord.setAttribute(SEQTAG, seq);
                    if(addQV){
                        String qv = new String((byte[]) localFastqLoader.getMapQV().get(name));
                        localSAMRecord.setAttribute(QVTAG, qv);
                    }
                }
                else{
                    log.info(new Object[]{"Error: read " + name + " not found !!!"});
                }
                
                localSAMFileWriter.addAlignment(localSAMRecord);
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception localException1) {
            localException1.printStackTrace();
            try {
                localSamReader.close();
                localSAMFileWriter.close();
            } catch (Exception localException2) {
                System.err.println("can not close stream");
            }
        } finally {
            try {
                localSamReader.close();
                localSAMFileWriter.close();
            } catch (Exception localException3) {
                System.err.println("can not close stream");
            }
        }
        
        return 0;
    }

    public static void main(String[] paramArrayOfString) {
        System.exit(new AddBamReadSequenceTag().instanceMain(paramArrayOfString));
    }
}
