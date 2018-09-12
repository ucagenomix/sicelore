package org.ipmc.sicelore.cmdline;

import java.util.ArrayList;
import java.util.List;
import org.ipmc.sicelore.programs.IlluminaOxfordBCUmiMerger;
import org.ipmc.sicelore.programs.BamSerializer;
import picard.cmdline.PicardCommandLine;

public class SiCeLoReMain extends PicardCommandLine {

    private static final String COMMAND_LINE_NAME = SiCeLoReMain.class.getSimpleName();

    public SiCeLoReMain() {
    }

    protected static List<String> getPackageList() {
        List<String> packageList = new ArrayList();
        packageList.add("org.ipmc.sicelore");
        return packageList;
    }

    public static void main(String[] args) {
        if (args.length > 0 && args[0].equals("IlluminaOxfordBCUmiMerger")) {
            String[] argument = new String[args.length - 1];
            int j = 0;
            for (int i = 1; i < args.length; i++) {
                argument[j] = args[i];
                j++;
            }
            new IlluminaOxfordBCUmiMerger(argument);
        } else if (args.length > 0 && args[0].equals("BamReader")) {
            String[] argument = new String[args.length - 1];
            int j = 0;
            for (int i = 1; i < args.length; i++) {
                argument[j] = args[i];
                j++;
            }
            new BamSerializer(argument);
        } else {
            System.exit(new SiCeLoReMain().instanceMain(args, getPackageList(), COMMAND_LINE_NAME));
        }
    }
}
