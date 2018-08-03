package org.ipmc.sicelore.programs;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.BufferedReader;
import java.io.File;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary="Split a bam according to the cell types provided in .csv file", oneLineSummary="Split a bam according to the cell types provided in .csv file", programGroup=org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SplitBamPerCluster extends CommandLineProgram
{
  private final Log log;
  private ProgressLogger pl;
  @Argument(shortName="I", doc="The input SAM or BAM file to analyze")
  public File INPUT;
  @Argument(shortName="O", doc="The output directory")
  public File OUTPUT;
  @Argument(shortName="CSV", doc="The .csv cluster file")
  public File CSV;
  @Argument(shortName="CELL_FLAG", doc="The cell barcode flag (default CB)")
  public String CELL_FLAG = "CB";
  
  public SplitBamPerCluster()
  {
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
    try
    {
      BufferedReader localBufferedReader = new BufferedReader(new java.io.FileReader(CSV));
      str1 = localBufferedReader.readLine();
      str1.split(",");
      str1 = localBufferedReader.readLine();
      while (str1 != null) {
        String[] localObject = str1.split(",");

        localHashMap1.put(localObject[0], localObject[1]);
        localHashMap1.put(localObject[0] + "-1", localObject[1]);
        
        if (!localHashMap2.containsKey(localObject[1])) {
          localHashMap2.put(localObject[1], new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader, true, new File(OUTPUT.getAbsolutePath() + "/" + localObject[1] + ".bam")));
        }
        str1 = localBufferedReader.readLine();
      }
      localBufferedReader.close();
      
      for (Iterator localIterator = localSamReader.iterator(); localIterator.hasNext();) { 
    	SAMRecord localObject = (SAMRecord)localIterator.next();
        pl.record((SAMRecord)localObject);
        String str2 = (String)((SAMRecord)localObject).getAttribute(CELL_FLAG);
        String str3 = (String)localHashMap1.get(str2);
        
        if (str3 != null) {
          ((SAMFileWriter)localHashMap2.get((String)localHashMap1.get(str2))).addAlignment((SAMRecord)localObject);
        }
      }
      
      localSamReader.close();
      Object localObject = localHashMap2.keySet();
      Iterator localIterator = ((Set)localObject).iterator();
      while (localIterator.hasNext())
        ((SAMFileWriter)localHashMap2.get((String)localIterator.next())).close();
    } catch (Exception localException) {
      localException.printStackTrace();
    }
    return 0;
  }
  
  public static void main(String[] paramArrayOfString)
  {
    System.exit(new SplitBamPerCluster().instanceMain(paramArrayOfString));
  }
}
