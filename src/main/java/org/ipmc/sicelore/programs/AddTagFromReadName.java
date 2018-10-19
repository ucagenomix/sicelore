/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import org.ipmc.sicelore.utils.*;
import gnu.trove.THashMap;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Add a sam flag from read name", oneLineSummary = "Add a sam flag from read name", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class AddTagFromReadName extends CommandLineProgram {

    private final Log log;
    private ProgressLogger pl;
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output BAM file")
    public File OUTPUT;
    @Argument(shortName = "TAG", doc = "The tag <default=XC>")
    public String TAG;

    public AddTagFromReadName()
    {
        log = Log.getInstance(AddTagFromReadName.class);
        pl = new ProgressLogger(log);
        TAG = "XC";
    }

    protected int doWork()
    {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        SamReader localSamReader = SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader localSAMFileHeader = localSamReader.getFileHeader();
        SAMFileWriter localSAMFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, OUTPUT);
        try {
            for (SAMRecord localSAMRecord : localSamReader) {
                pl.record(localSAMRecord);
                String name = localSAMRecord.getReadName();
                String[] tags = name.split(":");
                
                localSAMRecord.setAttribute(TAG, tags[7]);
                localSAMFileWriter.addAlignment(localSAMRecord);
            }
            localSamReader.close();
            localSAMFileWriter.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        return 0;
    }

    public static void main(String[] paramArrayOfString)
    {
        System.exit(new AddTagFromReadName().instanceMain(paramArrayOfString));
    }
}

