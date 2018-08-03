package org.ipmc.sicelore.programs;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary="Split a bam cells-by-cells provided in .csv file", oneLineSummary="Split a bam cells-by-cells provided in .csv file", programGroup=org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class SplitBamPerCell extends CommandLineProgram
{
  private final Log log;
  private ProgressLogger pl;
  @Argument(shortName="I", doc="The input SAM or BAM file to analyze")
  public File INPUT;
  @Argument(shortName="O", doc="The output directory")
  public File OUTPUT;
  @Argument(shortName="CSV", doc="The .csv cell barcode file (barcodes.tsv)")
  public File CSV;
  @Argument(shortName="CELL_FLAG", doc="The cell barcode flag (default CB)")
  public String CELL_FLAG;
  
  public SplitBamPerCell()
  {
    log = Log.getInstance(SplitBamPerCell.class);
    pl = new ProgressLogger(log);
  }

  protected int doWork()
  {
    String str1 = null;
    
    IOUtil.assertFileIsReadable(INPUT);
    IOUtil.assertFileIsReadable(CSV);
    
    HashMap localHashMap = new HashMap();
    
    SamReader localSamReader = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
    SAMFileHeader localSAMFileHeader1 = localSamReader.getFileHeader();
    SAMFileHeader localSAMFileHeader2 = localSAMFileHeader1.clone();
    localSAMFileHeader2.setSortOrder(htsjdk.samtools.SAMFileHeader.SortOrder.coordinate);

    try{
      BufferedReader localBufferedReader = new BufferedReader(new java.io.FileReader(CSV));
      str1 = localBufferedReader.readLine();
      str1.split(",");
      str1 = localBufferedReader.readLine();
      while (str1 != null) {
        String[] localObject = str1.split(",");
        
        if (!java.util.regex.Pattern.matches(".*-1", localObject[0])) {
          localObject[0] = (localObject[0] + "-1");
        }
        if (!localHashMap.containsKey(localObject[0])) {
          localHashMap.put(localObject[0], new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(localSAMFileHeader2, true, new File(OUTPUT.getAbsolutePath() + "/" + localObject[0] + ".bam")));
        }
        str1 = localBufferedReader.readLine();
      }
      localBufferedReader.close();
      
      for (Iterator localIterator = localSamReader.iterator(); localIterator.hasNext();) { 
    	  SAMRecord localObject = (SAMRecord)localIterator.next();
        pl.record((SAMRecord)localObject);
        String str2 = (String)((SAMRecord)localObject).getAttribute(CELL_FLAG);
        if ((SAMFileWriter)localHashMap.get(str2) != null) {
          ((SAMFileWriter)localHashMap.get(str2)).addAlignment((SAMRecord)localObject);
        }
      }
      localSamReader.close();
      
      Object localObject = localHashMap.keySet();
      Iterator localIterator = ((Set)localObject).iterator();
      while (localIterator.hasNext())
        ((SAMFileWriter)localHashMap.get((String)localIterator.next())).close();
    } catch (Exception localException) {
      localException.printStackTrace(); }
    return 0;
  }
  
  public static void main(String[] paramArrayOfString)
  {
    System.exit(new SplitBamPerCell().instanceMain(paramArrayOfString));
  }
}
