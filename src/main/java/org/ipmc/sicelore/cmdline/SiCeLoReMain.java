package org.ipmc.sicelore.cmdline;

import java.util.ArrayList;
import java.util.List;
import org.ipmc.sicelore.programs.IlluminaOxfordBCUmiMerger;
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
        if(args.length > 0 && args[0].equals("IlluminaOxfordBCUmiMerger")) {
            new IlluminaOxfordBCUmiMerger(args);
        }
        else {
            System.exit(new SiCeLoReMain().instanceMain(args, getPackageList(), COMMAND_LINE_NAME));
        }
    }
}
