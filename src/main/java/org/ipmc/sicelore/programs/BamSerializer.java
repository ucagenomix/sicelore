package org.ipmc.sicelore.programs;

/**
 * 
 * @author rainer waldmann
 * 
 */
import java.io.File;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import picard.cmdline.CommandLineProgram;

@CommandLineProgramProperties(summary = "Possorted_genome 10x bam file preprocessor (c.f.IlluminaParser-1.0.jar).", oneLineSummary = "Possorted_genome 10x bam file preprocessor (c.f.IlluminaParser-1.0.jar).", programGroup = org.ipmc.sicelore.cmdline.SiCeLoReUtilsExternal.class)
@DocumentedFeature
public class BamSerializer extends CommandLineProgram {

    /**
     *
     * @param args the command line arguments
     *
     */
    public BamSerializer(String[] args) {

        File inFile;

        File outFile;

        Options options = cli_otions();

        CommandLineParser parser = new DefaultParser();

        CommandLine cmd = null;

        try {

            cmd = parser.parse(options, args);

        } catch (ParseException ex) {

            HelpFormatter formatter = new HelpFormatter();

            formatter.printHelp("10x bam preprocessor\n parses 10x bam files and creates Hashmaps that are used to identify corresponding Nanopore reads ", options);

            System.exit(1);

        }

        inFile = new File(cmd.getOptionValue("i"));

        if (inFile.exists() == false) {

            System.out.println(inFile.getName() + "input File/Directory does not exist");

            System.exit(1);

        }

        outFile = new File(cmd.getOptionValue("o"));

        if (outFile.exists()) {

            System.out.println(outFile.getName() + " !!!!!!!!!!!!!!!!!!!!!!!!!!!!  output File exists !!!!!!!!!!!!!!!!!!");

            System.exit(1);

        }

        Integer nCells = null;

        if (cmd.hasOption("n")) {

            nCells = Integer.parseInt(cmd.getOptionValue("n"));

        }

        String tsv = cmd.getOptionValue("t");

        File tsvFile = null;

        if (tsv != null) {

            tsvFile = new File(tsv);

        }

        if (nCells == null && tsvFile == null || nCells != null && tsvFile != null) {

            HelpFormatter formatter = new HelpFormatter();

            formatter.printHelp("Must provide either n cells to use (-n) or cellranger tsv file(-t)", options);

            System.exit(1);

        }

        /*x

        new Parser(inFile, outFile, tsvFile, nCells, cmd.hasOption("x")).parse();

        x*/
    }

    /**
     *
     *
     *
     * @return
     *
     */
    private static Options cli_otions() {

        Options options = new Options();

        options.addOption(Option.builder("i").
                longOpt("inFileIllumina").
                required(true).
                desc("full path of 10x bam file").
                numberOfArgs(1)
                .build());

        options.addOption(Option.builder("o").
                longOpt("outFile").
                required(true).
                desc("full path of object output file").
                numberOfArgs(1)
                .build());

        options.addOption(Option.builder("n").
                longOpt("nCells").
                required(false).
                desc("use n cells with the most umis\nsupply this or the 10x tsv file with the list of cell BCs"
                        + "bam file contains all barcodes also from drops without cells").
                numberOfArgs(1)
                .build());

        options.addOption(Option.builder("t").
                longOpt("tsv").
                required(false).
                desc("use this 10x tsv to define the cell barcodes to use \n"
                        + "supply this or the 10x tsv file with the list of cell BCs\n"
                        + "bam file contains all barcodes also from drops without cells").
                numberOfArgs(1)
                .build());

        options.addOption(Option.builder("x")
                .longOpt("xml")
                .required(false)
                .desc("output xml file, gets pretty big and xml generation slows down program")
                .hasArg(false)
                .build());

        return options;

    }

    @Override
    protected int doWork() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
