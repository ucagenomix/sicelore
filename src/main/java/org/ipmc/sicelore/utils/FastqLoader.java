package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import gnu.trove.THashMap;
import java.io.BufferedReader;
import java.io.File;

public class FastqLoader {

    THashMap<String, byte[]> map = new THashMap();
    THashMap<String, byte[]> mapQV = new THashMap();

    public FastqLoader(File f) {
        map = new THashMap();
        mapQV = new THashMap();
        String str1 = null;
        String str2 = null;
        String qv = null;
        
        try {
            BufferedReader br = new BufferedReader(new java.io.FileReader(f));
            str2 = br.readLine();
            while (str2 != null) {
                str1 = br.readLine();
                br.readLine();
                qv = br.readLine();

                str2 = str2.replace("@", "");
                String[] arrayOfString = str2.split("\\s+");
                //System.out.println(arrayOfString[0]);
                map.put(arrayOfString[0], str1.getBytes());
                mapQV.put(arrayOfString[0], qv.getBytes());
                str2 = br.readLine();
            }
            br.close();
        } catch (Exception localException) { localException.printStackTrace(); }
    }

    public THashMap<String, byte[]> getMap() {
        return map;
    }
    public THashMap<String, byte[]> getMapQV() {
        return mapQV;
    }
}
