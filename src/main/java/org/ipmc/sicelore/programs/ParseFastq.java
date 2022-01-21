package org.ipmc.sicelore.programs;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.util.*;
import java.io.*;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Parse passed .fastq file and export cDNA.", oneLineSummary = "Parse passed .fastq file and export cDNA.", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtils.class)
@DocumentedFeature
public class ParseFastq extends CommandLineProgram {

    private final Log log;
    @Argument(shortName = "I", doc = "The input .fastq directory")
    public File FASTQDIR;
    @Argument(shortName = "O", doc = "The output .fastq directory")
    public File OUTPUT;
    @Argument(shortName = "offset", doc = "Barcodes length after adaptor end flag (default=28) removed from output")
    public int offset = 28;
    @Argument(shortName = "min_cdna", doc = "Minimum cDNA length (default=20) otherwise keep original read sequence")
    public int min_cdna = 20;
            
    public ParseFastq() {
        log = Log.getInstance(ParseFastq.class);
    }

    protected int doWork()
    {
        //IOUtil.assertFileIsReadable(INPUT);
        //IOUtil.assertFileIsWritable(OUTPUT);
        BufferedOutputStream os = null;
        
        // 16.914 fastq files files parsed in 30 minutes.
        // need then to add the US tag to the bam file using AddBamReadSequenceTag and addQV=false
        // then ComputeConsensus on a per chromosome basis
        
        try {
            for(File file : FASTQDIR.listFiles()) {
                
                log.info(new Object[]{"processing..."+file});
                
                File out   = new File(OUTPUT.getAbsolutePath() + "/" + file.getName());
                os = new BufferedOutputStream(new java.io.FileOutputStream(out));
                
                BufferedReader br = new BufferedReader(new java.io.FileReader(file));
                String readName = br.readLine();
                while(readName != null){
                    int PAst=0;
                    int Aend=0;
                    String sequence = br.readLine();
                    br.readLine();
                    String qv = br.readLine();
                    readName = readName.replace("@", "");
                    String[] arr = readName.split("\\s+");
                    String[] arr2 = arr[0].split("_");
                    for(int i=0;i<arr2.length;i++){
                        String[] keyval = arr2[i].split("=");
                        if(keyval.length>1){
                            if(keyval[0].equals("PAst")){
                                PAst = new Integer(keyval[1]).intValue();
                            }
                            if(keyval[0].equals("AEnd")){
                                Aend = new Integer(keyval[1]).intValue();
                            }
                        }
                    }
                    
                    String cDNA = "";
                    if(PAst>0 && Aend > 0 && (PAst-1-(Aend+offset)>min_cdna)){
                        cDNA = sequence.substring(Aend+offset,PAst-1);
                    }
                    else{
                        cDNA = sequence;
                        //log.info(new Object[]{"error..."+readName});
                    }
                    
                    os.write(new String("@"+arr[0]+"\n"+cDNA+"\n+\n\n").getBytes());
                    
                    //map.put(arr[0], sequence.getBytes());
                    //mapQV.put(arr[0], qv.getBytes());
                    
                    readName = br.readLine();
                }
                br.close();
                os.close();
            }
        }
        catch (Exception localException1) { localException1.printStackTrace(); } 
        finally {}
        
        return 0;
    }

    
    public static void main(String[] paramArrayOfString) {
        System.exit(new ParseFastq().instanceMain(paramArrayOfString));
    }
}
