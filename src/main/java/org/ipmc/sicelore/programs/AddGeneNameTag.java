package org.ipmc.sicelore.programs;

/**
 *
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import java.io.*;
import java.util.*;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.Log;
import htsjdk.samtools.util.ProgressLogger;
import htsjdk.samtools.util.Interval;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;
import picard.annotation.GeneAnnotationReader;
import htsjdk.samtools.util.OverlapDetector;
import picard.annotation.Gene;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecordIterator;
import picard.annotation.LocusFunction;

@CommandLineProgramProperties(summary = "Add a gene name tag to SAM records case it overlaps an exon", oneLineSummary = "Add a gene name tag to SAM records case it overlaps an exon.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoRe.class)
@DocumentedFeature
public class AddGeneNameTag extends CommandLineProgram
{ 
    @Argument(shortName = "I", doc = "The input SAM or BAM file")
    public File INPUT;
    @Argument(shortName = "O", doc = "The output SAM or BAM file with tags")
    public File OUTPUT;
    @Argument(shortName = "REFFLAT", doc = "The refFlat gene model file")
    public File REFFLAT;
    @Argument(shortName = "GENETAG", doc = "Gene name tag (default=GE)", optional=true)
    public String GENETAG = "GE";
    @Argument(shortName = "STRANDTAG", doc = "Gene strand tag (default=GS)", optional=true)
    public String STRANDTAG = "GS";
    @Argument(shortName = "FUNCTIONTAG", doc = "Gene function tag (default=XF)", optional=true)
    public String FUNCTIONTAG = "XF";
    @Argument(shortName = "USE_STRAND_INFO", doc = "Whether or ot nuse the strand information", optional=true)
    public boolean USE_STRAND_INFO = true;
    @Argument(shortName = "ALLOW_MULTI_GENE_READS", doc = "Whether or ot allow for multi-gene reads", optional=true)
    public boolean ALLOW_MULTI_GENE_READS = true;
    
    private ProgressLogger pl;
    private final Log log;
    private ReadTaggingMetric metrics;

    public AddGeneNameTag() {
        log = Log.getInstance(AddGeneNameTag.class);
        pl = new htsjdk.samtools.util.ProgressLogger(log, 1000000, "\tProcessed\t", "Records");
        metrics = new ReadTaggingMetric();  
    }

    protected int doWork() {
        IOUtil.assertFileIsReadable(REFFLAT);
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);
        
        //System.setProperty("samjdk.use_async_io_read_samtools", "true");
        //System.setProperty("samjdk.samjdk.use_async_io_write_samtools", "true");
        
        System.out.println(htsjdk.samtools.Defaults.allDefaults());
        
        process();

        return 0;
    }

    protected void process()
    {
        
        
        htsjdk.samtools.SamReader inputSam = htsjdk.samtools.SamReaderFactory.makeDefault().open(INPUT);
        htsjdk.samtools.SAMFileHeader samFileHeader = inputSam.getFileHeader();
        samFileHeader.setSortOrder(SAMFileHeader.SortOrder.unsorted);
        SAMFileWriter samFileWriter = new SAMFileWriterFactory().makeSAMOrBAMWriter(samFileHeader, true, OUTPUT);
        
        OverlapDetector<Gene> geneOverlapDetector = GeneAnnotationReader.loadRefFlat(REFFLAT, samFileHeader.getSequenceDictionary());
        log.info(new Object[]{"Loaded " + geneOverlapDetector.getAll().size() + " transcripts."});
        
        try{
            for(SAMRecord r : inputSam){
                pl.record(r);
                //log.info(new Object[]{"processing :" + r.getReadName()});
                
                if (!r.getReadUnmappedFlag())
                    r = setGeneExons(r, geneOverlapDetector);
                
                samFileWriter.addAlignment(r);
            }
            inputSam.close();
            samFileWriter.close();
        } catch (Exception e) { e.printStackTrace(); }
        
        log.info(new Object[] { metrics.toString() });
    }
  
    public SAMRecord setGeneExons(SAMRecord r, OverlapDetector<Gene> geneOverlapDetector)
    {
        Map<Gene, LocusFunction> map = getLocusFunctionForReadByGene(r, geneOverlapDetector);
        
        //log.info(new Object[]{"map.size :" + map.size() + "," + map.keySet()});
        
        java.util.Set<Gene> exonsForRead = getConsistentExons(r, map.keySet(), ALLOW_MULTI_GENE_READS);
        
        //log.info(new Object[]{"exonsForRead.size :" + exonsForRead.size()});
        
        List<Gene> genes = new ArrayList();
        for (Gene g : exonsForRead){
            LocusFunction f = (LocusFunction)map.get(g);
            if ((f == LocusFunction.CODING) || (f == LocusFunction.UTR)){
                genes.add(g);
                //log.info(new Object[]{"add gene: " + g});
            }
        }
        LocusFunction f = getLocusFunction(map.values());

        if(USE_STRAND_INFO) {
            genes = getGenesConsistentWithReadStrand(genes, r);
        }
        //if ((genes.size() > 1) && (!ALLOW_MULTI_GENE_READS)) {
        //    log.error(new Object[] { "There should only be 1 gene assigned to a read for DGE purposes." });
        //}
        
        String finalGeneStrand = getCompoundStrand(genes);
        String finalGeneName = getCompoundGeneName(genes);
        
        if (f != null)
            r.setAttribute(FUNCTIONTAG, f.toString());
        
        if ((finalGeneName != null) && (finalGeneStrand != null)) {
            
            r.setAttribute(GENETAG, finalGeneName);
            //log.info(new Object[]{"finalGeneName :" + finalGeneName});
            r.setAttribute(STRANDTAG, finalGeneStrand);
        }
        else {
            r.setAttribute(GENETAG, null);
            r.setAttribute(STRANDTAG, null);
        }
        return r;
    }
    
    private List<Gene> getGenesConsistentWithReadStrand(List<Gene> genes, SAMRecord r)
    {
        metrics.TOTAL_READS += 1;
        
        List<Gene> sameStrand = new ArrayList();
        List<Gene> oppositeStrand = new ArrayList();
        boolean negativeStrandRead = r.getReadNegativeStrandFlag();
        for (Gene g : genes) {
          boolean geneNegativeStrand = g.isNegativeStrand();

          if (((negativeStrandRead) && (geneNegativeStrand)) || ((!negativeStrandRead) && (!geneNegativeStrand))) {
            sameStrand.add(g);
          } else {
            oppositeStrand.add(g);
          }
        }
        if(sameStrand.isEmpty() && (oppositeStrand.size() > 0)) {
            metrics.READS_WRONG_STRAND += 1;
            return new ArrayList();
        }

        if (sameStrand.size() > 1) {
            // update KL 21/04/2020 --> authorize same strand gene ambiguous, we be resolve latter with isoform attribution
            //metrics.AMBIGUOUS_READS_REJECTED += 1;
            //return new ArrayList();
        }
        if (oppositeStrand.size() > 0) {
          metrics.READ_AMBIGUOUS_GENE_FIXED += 1;
        }
        metrics.READS_RIGHT_STRAND += 1;

        return sameStrand;
    }
 
    public Set<Gene> getConsistentExons(SAMRecord rec, Set<Gene> genes, boolean allowMultiGeneReads)
    {
        Set<Gene> result = new HashSet();
        String refName = rec.getReferenceName();
        List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
        for (AlignmentBlock b : alignmentBlocks) {
            Set<Gene> blockGenes = getAlignmentBlockonGeneExon(refName, b, genes);
            
            //if ((result.size() > 0) && (blockGenes.size() > 0)) {
            
            if (blockGenes.size() > 0) { // KL 19/09/2019
                if (allowMultiGeneReads) result.addAll(blockGenes);
                if (!allowMultiGeneReads) result.retainAll(blockGenes);
            }
            //else {
            //    result = blockGenes;
            //}
        }
        
        //log.error(new Object[] { "result.size" + result.size() });
        return result;
    }

    private Set<Gene> getAlignmentBlockonGeneExon(String refName, AlignmentBlock b, Set<Gene> genes)
    {
        Set<Gene> result = new HashSet();
        for (Gene g : genes) {
            if (getAlignmentBlockOverlapsExon(refName, b, g)) {
                result.add(g);
            }
        }
        //log.error(new Object[] { "result.size" + result.size() });
        return result;
    }
  
    private boolean getAlignmentBlockOverlapsExon(String refName, AlignmentBlock b, Gene g)
    {
        Interval ib = getInterval(refName, b);
        //log.error(new Object[] { "block interval" + ib.toString() });
        
        for (Gene.Transcript t : g) {
            for (Gene.Transcript.Exon e : t.exons) {
                Interval ei = getInterval(refName, e);
                
                if (ib.intersects(ei)) {
                    //log.error(new Object[] { "intersect gene interval" + ei.toString() });
                    return true;
                }
            }
        }
        return false;
    }
  
    private Interval getInterval(String refName, Gene.Transcript.Exon e)
    {
        Interval i = new Interval(refName, e.start, e.end);
        return i;
    }
  
    private Interval getInterval(String refName, AlignmentBlock b)
    {
        int s = b.getReferenceStart();
        int e = s + b.getLength() - 1;
        Interval i = new Interval(refName, s, e);
        return i;
    }
  
    public Map<Gene, LocusFunction> getLocusFunctionForReadByGene(SAMRecord rec, OverlapDetector<Gene> geneOverlapDetector)
    {
        Map<Gene, LocusFunction> result = new HashMap();
        Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
        Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(readInterval);

        for (Gene g : overlappingGenes) {
            LocusFunction f = getLocusFunctionForRead(rec, g);
            result.put(g, f);
        }
        return result;
    }

    private LocusFunction getLocusFunctionForRead(SAMRecord rec, Gene g)
    {
        List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();

        LocusFunction[] blockSummaryFunction = new LocusFunction[alignmentBlocks.size()];
        Set<Gene> temp = new HashSet();
        temp.add(g);

        for (int i = 0; i < alignmentBlocks.size(); i++) {
            AlignmentBlock alignmentBlock = (AlignmentBlock)alignmentBlocks.get(i);

            LocusFunction[] blockFunctions = getLocusFunctionsByBlock(alignmentBlock, temp);
            LocusFunction blockFunction = getLocusFunction(blockFunctions);
            blockSummaryFunction[i] = blockFunction;
        }
        LocusFunction readFunction = getLocusFunction(blockSummaryFunction);
        return readFunction;
    }
    
    public LocusFunction getLocusFunctionForRead(SAMRecord rec, OverlapDetector<Gene> geneOverlapDetector)
    {
        Interval readInterval = new Interval(rec.getReferenceName(), rec.getAlignmentStart(), rec.getAlignmentEnd());
        Collection<Gene> overlappingGenes = geneOverlapDetector.getOverlaps(readInterval);
        List<AlignmentBlock> alignmentBlocks = rec.getAlignmentBlocks();
        LocusFunction[] blockSummaryFunction = new LocusFunction[alignmentBlocks.size()];

        for (int i = 0; i < alignmentBlocks.size(); i++) {
            AlignmentBlock alignmentBlock = (AlignmentBlock)alignmentBlocks.get(i);
            LocusFunction[] blockFunctions = getLocusFunctionsByBlock(alignmentBlock, overlappingGenes);
            LocusFunction blockFunction = getLocusFunction(blockFunctions);
            blockSummaryFunction[i] = blockFunction;
        }

        LocusFunction readFunction = getLocusFunction(blockSummaryFunction);
        return readFunction;
    }
    
    private String getCompoundGeneName(Collection<Gene> genes)
    {
        if (genes.isEmpty()) { return null; }
    
        StringBuilder result = new StringBuilder();
        Iterator<Gene> iter = genes.iterator();
        result.append(((Gene)iter.next()).getName());
    
        while (iter.hasNext()) {
            result.append(",");
            result.append(((Gene)iter.next()).getName());
        }
        return result.toString();
    }
  
    private String getCompoundStrand(Collection<Gene> genes)
    {
        if (genes.isEmpty()) { return null; }
        StringBuilder result = new StringBuilder();
        Iterator<Gene> iter = genes.iterator();
        result.append((((Gene)iter.next()).isPositiveStrand())?"+":"-");
    
        while (iter.hasNext()) {
            result.append(",");
            result.append((((Gene)iter.next()).isPositiveStrand())?"+":"-");
        }
        return result.toString();
    }
    
    public LocusFunction getLocusFunction(Collection<LocusFunction> locusFunctions)
    {
        if (locusFunctions.size() == 0) return LocusFunction.INTERGENIC;
        LocusFunction[] array = (LocusFunction[])locusFunctions.toArray(new LocusFunction[locusFunctions.size()]);
        return getLocusFunction(array);
    }
  
    public LocusFunction getLocusFunction(LocusFunction[] locusFunctions)
    {
        Map<LocusFunction, Integer> functionScores;
        functionScores = new HashMap();
        functionScores.put(LocusFunction.CODING, new Integer(4));
        functionScores.put(LocusFunction.UTR, new Integer(3));
        functionScores.put(LocusFunction.INTRONIC, new Integer(2));
        functionScores.put(LocusFunction.INTERGENIC, new Integer(1));
      
        if ((locusFunctions == null) || (locusFunctions.length == 0)) { return null; }
        int bestScore = Integer.MIN_VALUE;
        LocusFunction bestFunction = LocusFunction.INTERGENIC;
        for (LocusFunction f : locusFunctions) {
            int score = ((Integer)functionScores.get(f)).intValue();
            if (score > bestScore) {
                bestScore = score;
                bestFunction = f;
            }
        }
        return bestFunction;
    }
  
    public LocusFunction[] getLocusFunctionsByBlock(AlignmentBlock b, Collection<Gene> overlappingGenes)
    {
        LocusFunction[] locusFunctions = new LocusFunction[b.getLength()];
        Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);
    
        for (Gene gene : overlappingGenes) {
            for (Gene.Transcript transcript : gene) {
                transcript.assignLocusFunctionForRange(b.getReferenceStart(), locusFunctions);
            }
        }
        return locusFunctions;
    }
  
    public LocusFunction[] getLocusFunctionsByInterval(Interval i, Collection<Gene> genes)
    {
        LocusFunction[] locusFunctions = new LocusFunction[i.length()];
        Arrays.fill(locusFunctions, 0, locusFunctions.length, LocusFunction.INTERGENIC);
    
        for (Gene gene : genes) {
            for (Gene.Transcript transcript : gene) {
                transcript.assignLocusFunctionForRange(i.getStart(), locusFunctions);
            }
        }
        return locusFunctions;
    }
    
    public class ReadTaggingMetric extends htsjdk.samtools.metrics.MetricBase
    {
      public int TOTAL_READS = 0;
      public int READS_WRONG_STRAND = 0;
      public int READS_RIGHT_STRAND = 0;
      public int READ_AMBIGUOUS_GENE_FIXED = 0;
      public int AMBIGUOUS_READS_REJECTED = 0;

      public ReadTaggingMetric() {}

      public String toString() { return "TOTAL READS [" + TOTAL_READS + "] CORRECT_STRAND [" + READS_RIGHT_STRAND + "]  WRONG_STRAND [" + READS_WRONG_STRAND + "] AMBIGUOUS_STRAND_FIXED [" + READ_AMBIGUOUS_GENE_FIXED + "] AMBIGUOUS REJECTED READS [" + AMBIGUOUS_READS_REJECTED + "]"; }
    }

    public static void main(String[] args) {
        System.exit(new AddGeneNameTag().instanceMain(args));
    }
}
