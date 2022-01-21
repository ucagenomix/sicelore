package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import gnu.trove.THashMap;
import htsjdk.samtools.util.Log;
import java.io.BufferedReader;
import java.io.File;

public class FastqLoader {

    private final Log log;
    THashMap<String, byte[]> map = new THashMap();
    THashMap<String, byte[]> mapQV = new THashMap();

    public FastqLoader(File directory, boolean addQV)
    {
        map = new THashMap();
        mapQV = new THashMap();
        String str1 = null;
        String str2 = null;
        String qv = null;
        
        log = Log.getInstance(FastqLoader.class);
        
        /*
        LineIterator it = FileUtils.lineIterator(theFile, "UTF-8");
        try {
            while (it.hasNext()) {
                String line = it.nextLine();
                // do something with line
            }
        } finally {
            LineIterator.closeQuietly(it);
        }
        */
        
        /*
        Stream<String> lines = Files.lines(Paths.get("c:\myfile.txt"));
        lines.forEach(l -> {
               // Do anything line by line   
        });
        */
        
        try {
            int i=0;
            for(File file : directory.listFiles()) {
                i++;
                if(i%250 == 0)
                    log.info(new Object[]{"read "+i+" files"});
                
                BufferedReader br = new BufferedReader(new java.io.FileReader(file));
                str2 = br.readLine();
                while (str2 != null) {
                    str1 = br.readLine();
                    br.readLine();
                    qv = br.readLine();

                    str2 = str2.replace("@", "");
                    String[] arrayOfString = str2.split("\\s+");
                    //System.out.println(arrayOfString[0]);
                    map.put(arrayOfString[0], str1.getBytes());
                    if(addQV)
                        mapQV.put(arrayOfString[0], qv.getBytes());
                    str2 = br.readLine();
                }
                br.close();
            }
            
        } catch (Exception localException) { localException.printStackTrace(); }
    }

    public THashMap<String, byte[]> getMap() {
        return map;
    }
    public THashMap<String, byte[]> getMapQV() {
        return mapQV;
    }
}
