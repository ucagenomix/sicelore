package org.ipmc.sicelore.programs;

/**
 *
 * @author rainer waldmann
 *
 */

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.LogManager;
import org.apache.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Illumina to Nanopore cellBC/UMI merger (c.f. NanoporeBC_UMI_finder-1.0.jar).", oneLineSummary = "Illumina to Nanopore cellBC/UMI merger (c.f. NanoporeBC_UMI_finder-1.0.jar).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtilsExternal.class)
@DocumentedFeature
public class IlluminaOxfordBCUmiMerger extends CommandLineProgram {

//    public static final Level LogLevel = Level.FINER;
//    public static final Level ConsoleLogLevel = Level.FINER;
    /**
     *
     * @param args the command line arguments
     *
     */
    private static final Logger LOGGER = LogManager.getLogger(IlluminaOxfordBCUmiMerger.class.getName());

    //static final Logger LOGGER = Logger.getLogger(IlluminaOxfordBCUmiMerger.class.getName());
    public IlluminaOxfordBCUmiMerger(String[] args) {

        Options options = cli_otions(/*xparametersx*/);

        CommandLineParser parser = new DefaultParser();

        CommandLine cmd = null;

        try {

            cmd = parser.parse(options, args);

        } catch (ParseException ex) {

            // Logger.getLogger(IlluminaOxfordBCUmiMerger.class.getName()).log(Level.SEVERE, null, ex);
            HelpFormatter formatter = new HelpFormatter();

            formatter.printHelp("Usage ", options);

            LOGGER.fatal("Command line parsing error:\n  ", ex);

            System.exit(1);

        }

        /*x      

         parameters.files.inFileNanopore = new File(cmd.getOptionValue("i"));

        if (parameters.files.inFileNanopore.exists() == false) {

            LOGGER.fatal("nanopore input File/Directory does not exist  " + parameters.files.inFileNanopore.getName());

            System.exit(1);

        }

        parameters.files.inFileIllumina = new File(cmd.getOptionValue("k"));

        if (parameters.files.inFileIllumina.exists() == false) {

            LOGGER.fatal(parameters.files.inFileIllumina.getName() + " input File does not exist");

            System.exit(1);

        }

               if (cmd.hasOption("o")) {

            parameters.files.outFile = new File(cmd.getOptionValue("o"));

        }

        if (parameters.files.outFile == null) {

            String path = System.getProperty("user.home") + File.separator + parameters.general.output_directory;

            File outDir = new File(path);

            if (outDir.exists() == false) {

                outDir.mkdirs();

            }

            parameters.files.outFile = new File(outDir, parameters.files.inFileNanopore.getName().substring(0, parameters.files.inFileNanopore.getName().lastIndexOf('.')) + parameters.general.output_filesuffix + ".bam");

            for (int i = 1; parameters.files.outFile.exists(); i++) {

                parameters.files.outFile = new File(outDir, parameters.files.inFileNanopore.getName().substring(0, parameters.files.inFileNanopore.getName().lastIndexOf('.')) + parameters.general.output_filesuffix + i + ".bam");

            }

            }

            String s = parameters.files.outFile.getPath();

            if(s.endsWith(".bam"))

                s = s.substring(0, s.lastIndexOf(".bam"));

            parameters.files.outFileUmiFoundOnly = new File(s+ "_umifound_" + ".bam");

       

        

        if (cmd.hasOption("l")) {

            if (cmd.getOptionValue("l") != null && cmd.getOptionValue("l").length() > 2) {

                parameters.general.logOutput = new File(cmd.getOptionValue("l"));

                if (parameters.general.logOutput.getParentFile() == null) {

                    File outDir = new File(System.getProperty("user.home") + File.separator + parameters.general.output_directory);

                    if (outDir.exists() == false) {

                        outDir.mkdirs();

                        parameters.general.logOutput = new File(outDir, cmd.getOptionValue("l"));              

                    }

                }

            } else { // no file name specified, write to input file + .log

                parameters.general.logOutput = new File(parameters.files.outFile.getPath() + ".log");

            }

            try {

                            parameters.general.logStream = new PrintStream(new FileOutputStream(parameters.general.logOutput));

                        } catch (FileNotFoundException ex) {

                            java.util.logging.Logger.getLogger(IlluminaOxfordBCUmiMerger.class.getName()).log(Level.SEVERE, null, ex);

                        }

        } else {

           parameters.general.logStream = System.out;

        }

 

 

        if (cmd.hasOption("h")) {

            HelpFormatter formatter = new HelpFormatter();

            formatter.printHelp("Usage ", options);

            System.exit(0);

        }

        if (cmd.hasOption("d")) {

            parameters.general.debugMode = true;

        }

 

        if (cmd.hasOption("q")) {

            parameters.adapter.maxAdapterNeedlemanMismatches = Integer.parseInt(cmd.getOptionValue('q'));

        }

        if (cmd.hasOption("p")) {

            parameters.adapter.setAdaptersequence(cmd.getOptionValue("p"));

        }

        if (parameters.adapter.adapterBytesequence == null) {

            HelpFormatter formatter = new HelpFormatter();

            formatter.printHelp("adapter sequence needs to be provided either in config.xml or command line", options);

            System.exit(1);

        }

        if (cmd.hasOption("s")) {

            parameters.tso.setTSOsequence(cmd.getOptionValue("s"));

        }

        if (cmd.hasOption("r")) {

           if (cmd.hasOption("s") == false) {

                 parameters.general.logStream.println("*** INFO -r option ignored since no TSO sequence provided");

            } else {

                //  parameters.tso.maxTSO_NeedlemanMismatches = Integer.parseInt(cmd.getOptionValue('r'));

                parameters.tso.maxTSO_NeedlemanMismatches = new Integer(cmd.getOptionValue('r'));

            }

        }

        if (parameters.tso.tsoBytesequence != null && parameters.tso.maxTSO_NeedlemanMismatches == null) {

            HelpFormatter formatter = new HelpFormatter();

            formatter.printHelp("must provide max TSO mismatches either in config.xml or in command line", options);

            System.exit(1);

        }

            if (cmd.hasOption("t")) {

                parameters.general.nThreads = Integer.parseInt(cmd.getOptionValue('t'));

            }

        if (cmd.hasOption("a")) {

            parameters.polyAT.windowSearchForPolyA = Integer.parseInt(cmd.getOptionValue('a'));

        }

        if (cmd.hasOption("b")) {

            parameters.barcodes.cellBC_editdistance = Integer.parseInt(cmd.getOptionValue('b'));

        }

        if (cmd.hasOption("e")) {

            parameters.barcodes.simulateRandomBCs = true;

            parameters.general.logStream.println("******************************************************************");

             parameters.general.logStream.println("***********  SIMULATING RANDOM BARCODES **************************");

             parameters.general.logStream.println("******************************************************************");

        }

        if (cmd.hasOption("f")) {

            parameters.umis.simulateRandomUmis = true;

             parameters.general.logStream.println("******************************************************************");

             parameters.general.logStream.println("***********  SIMULATING RANDOM UMIs ******************************");

             parameters.general.logStream.println("******************************************************************");

        }

        if (cmd.hasOption("u")) {

            parameters.umis.umi_editdistance = Integer.parseInt(cmd.getOptionValue('u'));

        }

        if (cmd.hasOption("j")) {

            parameters.umis.maxUMIfalseAssignmentPcnt = Integer.parseInt(cmd.getOptionValue('j'));

            parameters.umis.umiEditDistances = readDynamicEditDistances(Parameters.UMIs.edit_distance_xml, "umi edit distances",

                    parameters.umis.maxUMIfalseAssignmentPcnt);

                   //System.out.println("O.K. Found and using dynamic UMI edit distances: " + Parameters.UMIs.edit_distance_xml);

            if (parameters.umis.umiEditDistances.entries.get(parameters.umis.umi_length).containsKey(parameters.umis.maxUMIfalseAssignmentPcnt) == false) {

                 parameters.general.logStream.println(" max percent error in --maxUMIfalseMatchPercent should included in "

                        + Parameters.UMIs.edit_distance_xml);

                System.exit(1);

            }

        }

        if (cmd.hasOption("m")) {

            parameters.umis.incrementUMI_ED_ifFewUmis = true;

        }

        if (cmd.hasOption("n")) {

            parameters.barcodes.incrementBC_ED_ifFewBCs = true;

        }

        if (cmd.hasOption("c")) {

            parameters.barcodes.cell_bc_length = Integer.parseInt(cmd.getOptionValue('c'));

        }

        if (cmd.hasOption("v")) {

            parameters.umis.umi_length = Integer.parseInt(cmd.getOptionValue('v'));

        }

        if (cmd.hasOption("w")) {

            parameters.polyAT.polyATlength = Integer.parseInt(cmd.getOptionValue('w'));

        }

        if (cmd.hasOption("x")) {

            parameters.polyAT.fractionATInPolyAT = Float.parseFloat(cmd.getOptionValue('x'));

        }

        if (cmd.hasOption("y")) {

            parameters.barcodes.maxbcfalseAssignmentPcnt = Integer.parseInt(cmd.getOptionValue('y'));

 

            parameters.barcodes.bcEditDistances = readDynamicEditDistances(Parameters.Barcodes.edit_distance_xml, "barcode edit distances", parameters.barcodes.maxbcfalseAssignmentPcnt);

             //System.out.println("O.K. Found and using dynamic barcode edit distances: " + Parameters.Barcodes.edit_distance_xml.toString());

            if (parameters.barcodes.bcEditDistances.entries.get(parameters.barcodes.cell_bc_length).containsKey(parameters.barcodes.maxbcfalseAssignmentPcnt) == false) {

                 parameters.general.logStream.println(" max percent error in --maxBCfalseMatchPercent should included in "

                        + Parameters.Barcodes.edit_distance_xml);

                System.exit(1);

            }

        }

        if (cmd.hasOption("z")) {

            parameters.barcodes.cell_BC_bailout_after_ED = Integer.parseInt(cmd.getOptionValue('z'));

        }

        if (cmd.hasOption("g")) {

            String val = cmd.getOptionValue("g");

            if (val.length() != 2 || val.chars().allMatch(Character::isLetter) == false) {

                HelpFormatter formatter = new HelpFormatter();

                formatter.printHelp("-g option should have two letters", options);

                System.exit(1);

            } else {

                parameters.files.gene_name_attribute = cmd.getOptionValue('g');

            }

        }

        if (parameters.validate( parameters.general.logStream) == false) {

            System.exit(1);

        }

        //       validateArguments(parameters, options);

        NanoporeSeqAnalyzer worker = new NanoporeSeqAnalyzer(parameters, readIlluminaData(parameters));

        worker.go();

        // finish up

        parameters.general.logStream.close();

        x*/
    }

    /**
     *
     *
     *
     * @return
     *
     */
    private static Options cli_otions(/*Parameters params*/) {

        Options options = new Options();

        options.addOption(Option.builder("a")
                .required(false)
                .longOpt("polyawin")
                .desc("n bases from start to search polyA, defaults to :" /*x+ params.polyAT.windowSearchForPolyAx*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("b")
                .required(false)
                .longOpt("bcedit")
                .desc("max errors for cell barcode (edit distance, default 3")
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("c")
                .required(false)
                .longOpt("bclength")
                .desc("cell barcode length, default " /*x+ params.barcodes.cell_bc_lengthx*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("d")
                .required(false)
                .longOpt("debug")
                .desc("debug mode , allows stepwise execution")
                .numberOfArgs(0)
                .build());

        options.addOption(Option.builder("e").
                longOpt("randomBarcode").
                //argName("UMI seq").

                required(false).
                desc("Random Barcode Simulation").
                numberOfArgs(0)
                .build());

        options.addOption(Option.builder("f").
                longOpt("randomUMI").
                //argName("UMI seq").

                required(false).
                desc("Random UMI simulation").
                numberOfArgs(0)
                .build());

        options.addOption(Option.builder("g")
                // .argName("n cpu cores")

                .required(false)
                .longOpt("ONTgene")
                .desc("2 char SAM attribute for gene name in Nanopore SAM. Default: " /*x+ params.files.gene_name_attributex*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("h")
                .required(false)
                .longOpt("help")
                .desc("help")
                .numberOfArgs(0)
                .build());

        options.addOption(Option.builder("i")
                .required(true)
                .longOpt("inFileNanopore")
                .desc("Nanopore input file")
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("j")
                .required(false)
                .longOpt("maxUMIfalseMatchPercent")
                .desc("Maximal percentage of false UMI association. Will dynamically adjust UMI edit distance depending \n"
                        + "on the numbers of umis to compare with. Uses data from similation in file " /*x+ params.umis.getEdit_distance_xml()x*/ + " in current working directory and if not found there in application root\n"
                        + "overides fixed umi edit distance parameters")
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("k").
                longOpt("inFile10x").
                //argName("19X input file").

                required(true).
                desc("10x Illumina File parsed by 10xBamParser").
                numberOfArgs(1)
                .build());

        options.addOption(Option.builder("l").
                longOpt("logFile").
                required(false).
                desc("log output file, defaults to stdout. If just one character is given as name will create log with same name as outfile with .log appended").
                //numberOfArgs(1).

                build());

        options.addOption(Option.builder("m")
                .required(false)
                .longOpt("incrumiedit")
                .desc("increase umi edit distance by one for genes wit less than " + /*xparams.umis.maxUMIsforincrementUMI_ED + x*/ " umis for gene in given cell")
                .numberOfArgs(0)
                .build());

        options.addOption(Option.builder("n")
                .required(false)
                .longOpt("incrBCedit")
                .desc("increase BC edit distance by one for genes expressed in less than " + /*x params.barcodes.maxBCsforincrementBC_ED + x*/ " cells. Increase won't affect testing for all szlected 10x barcodes or empty drop barcodes")
                .numberOfArgs(0)
                .build());

        options.addOption(Option.builder("o")
                //.argName("output file")

                .required(false)
                .longOpt("outfile")
                .desc("output bam file, if ommited will write file (inputfilename + \"10xAttributes.bam\" into <user home>/" /*x+ params.general.output_directoryx*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("p").
                longOpt("adapterSeq").
                //argName("UMI seq").

                required(false).
                desc("adapter sequence preceeeding UMI or barcode").
                numberOfArgs(1)
                .build());

        options.addOption(Option.builder("q").
                longOpt("adapterMaxErr").
                //argName("UMI seq").

                required(false).
                desc("max mismatches in adapter (needleman alignment), default " /*x+ params.adapter.maxAdapterNeedlemanMismatchesx*/).
                numberOfArgs(1)
                .build());

        options.addOption(Option.builder("r")
                .required(false)
                .longOpt("tsMaxErr")
                .desc("max number of TSO mismatches, default: " /*x+ params.tso.maxTSO_NeedlemanMismatchesx*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("s")
                .required(false)
                .longOpt("tsoseq")
                .desc("TSO sequence, optional will be searched for and flagged if provided, not used for read selection")
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("t")
                // .argName("n cpu cores")

                .required(false)
                .longOpt("ncpu")
                .desc("n threads generated - if not provided will use all available cpu cores")
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("u")
                .required(false)
                .longOpt("umiedit")
                .desc("max errors for umi")
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("v")
                .required(false)
                .longOpt("umilength")
                .desc("umi length, default " /*x+ params.umis.umi_lengthx*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("w")
                .required(false)
                .longOpt("polyAlength")
                .desc("min polyA/T length, default " /*x+ params.polyAT.polyATlengthx*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("x")
                .required(false)
                .longOpt("polyATfraction")
                .desc("fraction  T  in poly T tail, default: " /*x+ params.polyAT.fractionATInPolyATx*/)
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("y")
                // .argName("n cpu cores")

                .required(false)
                .longOpt("maxBCfalseMatchPercent")
                .desc("Maximal percentage of false BC association. Will dynamically adjust BC edit distance depending \n"
                        + "on the numbers of umis to compare with. Uses data from similation in file " /*x+ Parameters.Barcodes.edit_distance_xmlx*/ + " in current working directory and if not found there in application root\n"
                        + "overides fixed BC edit distance parameters")
                .numberOfArgs(1)
                .build());

        options.addOption(Option.builder("z")
                .required(false)
                .longOpt("edBCbailout")
                .desc("mut cycles (ED) after which testing is aborted for Barcode when match was found.\n avoids unnecessary testing but results in less info on secondary matches.\n"
                        + "Defaults to always testing full ED cycles")
                .numberOfArgs(1)
                .build());

        return options;

    }

    @Override
    protected int doWork() {
        return 0;
    }

    /**
     *
     * read config file try first in dir where application is if not found try
     *
     * in current working directory
     *
     *
     *
     * @return
     *
     */

    /*x

    private static Parameters readConfigFile() {

        File configXml = checkConfigFilePath(Parameters.CONFIG_XML, "config");

        if (configXml == null) {

            LOGGER.fatal(Parameters.CONFIG_XML + " not found" + " ABORTING");

            System.exit(1);

        }

        Parameters parameters = JAXB.unmarshal(configXml, Parameters.class);

        try {

            JaxbValidator.validateRequired(parameters, Parameters.class);

        } catch (JaxbValidator.ValidationException ex) {

            LOGGER.fatal("TROUBLE WITH CONFIG.XML FORMAT !!!!!!!  EXITING");

        }

        // System.out.println("O.K. Found and used " + configXml.toString());

        parameters.files.configXMLFile = configXml.toString();

        return parameters;

    }

x*/
    /**
     *
     * searches for file in current working directory and if not found in
     *
     * application root
     *
     *
     *
     * @param fileName
     *
     * @param message Breaf description of file e.g. config ....
     *
     * @return
     *
     */

    /*x

    static private File checkConfigFilePath(String fileName, String message) {

PrintStream stream = System.out;

File retval;

if(fileName.contains(File.separator)){

    retval = new File(fileName);

    if(retval.exists() == false){

      stream.println(message + " " + fileName + " not found ABORTING !!!!");

      System.exit(1);

    }

    stream.println("USING " + message + " FILE: " + retval.getName());

} else {   

        retval = new File(new File("").getAbsolutePath(), fileName);

         stream.println("READING " + message + " FILE: " + retval.getName());

        stream.println("Trying current working directory: " + retval.toString());

        if (retval.exists() == false) {

            retval = new File(new File(IlluminaOxfordBCUmiMerger.class.getProtectionDomain().getCodeSource().getLocation().getPath()).getParent(), fileName);

            stream.println("Not found trying application root:" + retval.toString());

        }

        if (retval.exists() == false) {

            stream.println(message + " file " + fileName + "should be either in");

            stream.println("current working directory: " + new File("").getAbsolutePath());

            stream.println("or in application directory: " + retval.getPath());

 

            retval = null;

        }

        stream.println("O.K. Found and using " + retval.toString());

}

        return retval;

    }

x*/
    /**
     *
     * reads umiMaxEditDistances.xml if required and checks consistency of umi
     *
     * parameters
     *
     *
     *
     * @return
     *
     */

    /*x

    private static DynamicEditDistances readDynamicEditDistances(String fileName, String message, int maxErrorPercent) {

        File umiconfigXml = checkConfigFilePath(fileName, message);

        if (umiconfigXml == null) {

            LOGGER.fatal("Trouble reading " + message + " Umi edit distance config file" + umiconfigXml.toString() + "\nABORTING!!!!!");

            System.exit(1);

        }

        DynamicEditDistances editDistances = JAXB.unmarshal(umiconfigXml, DynamicEditDistances.class);

        try {

            JaxbValidator.validateRequired(editDistances, DynamicEditDistances.class);

        } catch (JaxbValidator.ValidationException ex) {

            LOGGER.fatal("Trouble reading " + message + " Umi edit distance config file" + umiconfigXml.toString());

            System.exit(1);

        }

        editDistances.fileName = umiconfigXml.toString();

        return editDistances;

    }

x*/
    /**
     *
     * Read Illumina data
     *
     *
     *
     * @param params
     *
     */

    /*x

    @SuppressWarnings({"unchecked"}) // otherwise warning: [unchecked] unchecked cast retval.allUMIsforEachCell = (Long2ObjectOpenHashMap<LongOpenHashSet>) ois.readObject();

    private static Illumina10xData readIlluminaData(Parameters params) {

        Illumina10xData retval = new Illumina10xData();

        ObjectInputStream ois = null;

        try {

            InputStream streamIn;

            streamIn = new BufferedInputStream(new FileInputStream(params.files.inFileIllumina), 1000000);

            ois = new ObjectInputStream(streamIn);

            retval.illuminaDataGene_CellBC_UMI = (IlluminaData) ois.readObject();

            retval.unusedCellBCs = (UnusedCellBCs) ois.readObject();

            long emptDropcount = retval.unusedCellBCs.values().stream().flatMap((e)->e.stream()).distinct().count();

             params.general.logStream.println("Empty Drop Cell barcodes in Illumina file:" + emptDropcount);

            retval.allpassed10xCellBCs = (LongOpenHashSet) ois.readObject();

            retval.allUMIsforEachCell = (Long2ObjectOpenHashMap<LongOpenHashSet>) ois.readObject();

//            Object allUMIsforEachCell = ois.readObject();

//            if(allUMIsforEachCell instanceof Long2ObjectOpenHashMap)

//                retval.allUMIsforEachCell =  (Long2ObjectOpenHashMap<LongOpenHashSet>) allUMIsforEachCell;

//            else

//               throw (new ClassNotFoundException());

 

        } catch (FileNotFoundException ex) {

            LOGGER.fatal(params.files.inFileIllumina.getName() + " input File does not exist");

        } catch (IOException | ClassNotFoundException ex) {

            LOGGER.fatal(params.files.inFileIllumina.getName() + " I/O exceptiont");

        } finally {

            if (ois != null) {

                try {

                    ois.close();

                } catch (IOException ex) {

                    LOGGER.fatal(params.files.inFileIllumina.getName() + "error closing file");

                }

            }

        }

        return retval;

    }

x*/
}
