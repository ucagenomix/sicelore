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

    public FastqLoader(File paramFile) {
        map = new THashMap();
        String str1 = null;
        String str2 = null;

        try {
            BufferedReader localBufferedReader = new BufferedReader(new java.io.FileReader(paramFile));
            str2 = localBufferedReader.readLine();
            while (str2 != null) {
                str1 = localBufferedReader.readLine();
                localBufferedReader.readLine();
                localBufferedReader.readLine();

                str2 = str2.replace("@", "");
                String[] arrayOfString = str2.split(" ");
                map.put(arrayOfString[0], str1.getBytes());
                str2 = localBufferedReader.readLine();
            }
            localBufferedReader.close();
        } catch (Exception localException) {
            localException.printStackTrace();
        }
    }

    public THashMap<String, byte[]> getMap() {
        return map;
    }
}
