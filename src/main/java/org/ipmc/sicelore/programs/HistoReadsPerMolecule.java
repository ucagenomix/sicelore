package org.ipmc.sicelore.programs;

import htsjdk.samtools.*;
import htsjdk.samtools.util.*;
import java.io.*;
import java.util.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary="Histogram of reads per molecule for cells in .csv", oneLineSummary="Histogram of reads per molecule for cells in .csv", programGroup=org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class HistoReadsPerMolecule extends CommandLineProgram
{
  private final Log log;
  private ProgressLogger pl;
  @Argument(shortName="I", doc="The input SAM or BAM file to analyze")
  public File INPUT;
  @Argument(shortName="O", doc="The output histogram.txt")
  public File OUTPUT;
  @Argument(shortName="CSV", doc="The cell barcodes .csv file")
  public File CSV;
  @Argument(shortName="GENE_TAG", doc="The gene name TAG (default IG, GN for 10x illumina bam)")
  public String GENE_TAG = "IG";
  @Argument(shortName="CELL_TAG", doc="The cell barcode TAG (default CB)")
  public String CELL_TAG = "CB";
  @Argument(shortName="UMI_TAG", doc="The UMI TAG (default U8, UB for 10x illumina bam)")
  public String UMI_TAG = "U8";
  
  HashMap<String, List<String>> map;
  public HashSet<String> DTEcells;
  
  public HistoReadsPerMolecule()
  {
    log = Log.getInstance(HistoReadsPerMolecule.class);
    pl = new ProgressLogger(log);
    DTEcells = new HashSet<String>();
    map = new HashMap<String, List<String>>();
  }

  protected int doWork()
  {
	  List<String> lst;
	  DataOutputStream os = null;
	  IOUtil.assertFileIsReadable(INPUT);
	  IOUtil.assertFileIsWritable(OUTPUT);
		
	  loadDTEcells();
	  log.info(new Object[] { "Cells number\t\t[" + DTEcells.size() + "]" });
		
	  htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
	  try{
		  
		  log.info(new Object[] { "Parsing bam file\t\tstart..."});
		  
		  for(SAMRecord r : inputSam){
			  pl.record(r);
			  String read_name = r.getReadName();
			  String gene = (String)r.getAttribute(GENE_TAG);
			  String barcode = (String)r.getAttribute(CELL_TAG);
			  String umi = (String)r.getAttribute(UMI_TAG);
			  
			  if(gene != null && barcode != null && umi != null && DTEcells.contains(barcode)){
				  if((lst = (List<String>)this.map.get(barcode+":"+umi)) == null)
					  map.put(barcode+":"+umi, new ArrayList());
		        	   
				  ((List<String>)map.get(barcode+":"+umi)).add(read_name);
			  }
		  }
		  inputSam.close();
		  
		  log.info(new Object[] { "Total molecules\t\t["+map.size()+"]"});
		  
		  int[] counter = new int[101];
		  Set cles = map.keySet();
		  Iterator<String> it = cles.iterator();
		  while(it.hasNext()){
				String key = (String)it.next();
				int nb = ((List<String>)map.get(key)).size();
				if(nb > 100){ nb=100;}
				counter[nb]++;
		  }
		  
		  os = new DataOutputStream(new FileOutputStream(OUTPUT));
		  os.writeBytes("xtimes\tnumber\n");
		  for(int i=1; i<101; i++){
			  os.writeBytes(i+"\t"+counter[i]+"\n");
			  log.info(new Object[] { "x"+i+"\t["+counter[i]+"]"});
		  }
		  os.close();
		  
		  
	  }catch(Exception e){ e.printStackTrace(); }
	  finally{ try{ inputSam.close(); os.close(); } catch(Exception e){System.err.println("can not close stream");} }
	  
	  return 0;
  
  }
	public void loadDTEcells()
	{
		try{
			BufferedReader fichier = new BufferedReader(new FileReader(CSV));
			String line = fichier.readLine();
			while(line != null){
				DTEcells.add(line + "-1");
				line = fichier.readLine();
			}
			fichier.close();
		}catch(Exception e){ e.printStackTrace(); }
	}
  
  public static void main(String[] paramArrayOfString)
  {
    System.exit(new HistoReadsPerMolecule().instanceMain(paramArrayOfString));
  }
}