package org.ipmc.sicelore.programs;

import java.io.*;
import java.util.*;
import java.util.regex.Pattern;
import htsjdk.samtools.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.CloseableIterator;
import org.ipmc.sicelore.utils.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary="Produce Isoforms Expression Matrix", oneLineSummary="Produce Isoforms Expression Matrix", programGroup=org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class IsoformExpressionMatrix extends CommandLineProgram
{
	@Argument(shortName="I", doc="The input SAM or BAM file")
	public File INPUT;
	@Argument(shortName="REFFLAT", doc="The refFlat gene model file")
	public File REFFLAT;
	@Argument(shortName="CSV", doc="The cell barcodes .csv file")
	public File CSV;
	@Argument(shortName="MATRIX", doc="The transcripts count matrix file")
	public File MATRIX;
	@Argument(shortName="DELTA", doc="Allowed base number difference between start/end of exons and read block position (default=10)")
	public int DELTA = 10;
	@Argument(shortName="SOFT", doc="Transcripts exons can be smaller than LongReadRecord exons (detection of specific alternative exons like flip/flop gria2 of Pkm1/Pkm2)")
	public boolean SOFT = false;
	@Argument(shortName="METRICS", doc="The output metrics file")
	public File METRICS;
	
	public HashSet<String> DTEcells;
	private final Log log;
	
	public IsoformExpressionMatrix()
	{
		log = Log.getInstance(IsoformExpressionMatrix.class);
		this.DTEcells = new HashSet<String>();
	}
	
	protected int doWork()
	{
		IOUtil.assertFileIsReadable(REFFLAT);
		IOUtil.assertFileIsReadable(INPUT);
		IOUtil.assertFileIsReadable(CSV);
		process();
		
		return 0;
	}
	
	protected void process()
	{
            loadDTEcells();
            log.info(new Object[] { "Cells number\t\t[" + DTEcells.size() + "]" });
		
            // 4mn and 9.6Gb for 1.450.000 SAMrecords [747.000 molecules]
            UCSCRefFlatParser model = new UCSCRefFlatParser(REFFLAT);
            LongreadParser bam = new LongreadParser(INPUT);
            MoleculeDataset dataset = new MoleculeDataset(bam, model, DELTA, SOFT);
            Matrix matrix = dataset.DTEMatrix(model);
            matrix.write(MATRIX, DTEcells);
            dataset.displayMetrics(METRICS);
	}
	
	public void loadDTEcells()
	{
		try{
			BufferedReader fichier = new BufferedReader(new FileReader(CSV));
			String line = fichier.readLine();
			while(line != null){
				DTEcells.add(line);
				line = fichier.readLine();
			}
			fichier.close();
		}catch(Exception e){ e.printStackTrace(); }
	}
	
	public static void main(String[] args)
	{
		System.exit(new IsoformExpressionMatrix().instanceMain(args));
	}
}