package org.ipmc.sicelore.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashSet;
import java.io.File;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public class CellList extends HashSet<String>
{
    public CellList(File CSV)
    {
        try {
            BufferedReader fichier = new BufferedReader(new FileReader(CSV));
            String line = fichier.readLine();
            while(line != null) {
                line=line.replace("-1","");
                this.add(line);
                line = fichier.readLine();
            }
            fichier.close();
        } catch (Exception e) { e.printStackTrace(); }
   }
}
