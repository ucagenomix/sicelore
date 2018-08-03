package org.ipmc.sicelore.programs;

import java.io.*;
import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.util.*;
import java.util.regex.Pattern;
import java.math.BigInteger;
import gnu.trove.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary="Filter for reads associated with an Illumina cell barcode and UMI", oneLineSummary="Filter for reads associated with an Illumina cell barcode and UMI", programGroup=org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class FilterGetBcUmiReads extends CommandLineProgram
{
	  private final Log log;
	  private htsjdk.samtools.util.ProgressLogger pl;
	  
	  @Argument(shortName="I", doc="The input SAM or BAM file")
	  public File INPUT;
	  @Argument(shortName="O", doc="The output BAM file")
	  public File OUTPUT;
    
	  public FilterGetBcUmiReads()
	  {
		  log = Log.getInstance(FilterGetBcUmiReads.class);
		  pl = new htsjdk.samtools.util.ProgressLogger(log);
	  }
	  
	  protected int doWork()
	  {
		  int ww=0;

		  IOUtil.assertFileIsReadable(INPUT);
		  IOUtil.assertFileIsWritable(OUTPUT);
			
		  htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
		  htsjdk.samtools.SAMFileHeader header = inputSam.getFileHeader();
		  SAMFileWriter outputSam = new htsjdk.samtools.SAMFileWriterFactory().makeSAMOrBAMWriter(header, true, OUTPUT);
			
		  try{
			  for(SAMRecord r : inputSam){
				  pl.record(r);
				  String read_id = r.getReadName();
				  String barcode = (String)r.getAttribute("BC");
				  String umi = (String)r.getAttribute("U8");
				  
				  if(barcode != null && umi != null)
					  outputSam.addAlignment(r);
			  }
			  inputSam.close();
			  outputSam.close();
		  }catch(Exception e){ e.printStackTrace(); }
		  finally{ try{inputSam.close(); outputSam.close();} catch(Exception e){System.err.println("can not close stream");} }
		  
		return 0;
	}
	
	public static void main(String[] args)
	{
		System.exit(new FilterGetBcUmiReads().instanceMain(args));
	}

}
