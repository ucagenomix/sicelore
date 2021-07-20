package org.ipmc.common.utils;

import java.io.BufferedReader;
import java.io.Closeable;

public class ExecuteCmd {

    private Runtime r = null;

    private Process p = null;

    private String[] commande;

    private String[] envp;

    private String pathExecution;

    public ExecuteCmd(String[] paramArrayOfString1, String[] paramArrayOfString2, String paramString)
    {
        commande = paramArrayOfString1;
        envp = paramArrayOfString2;
        pathExecution = paramString;
    }

    public boolean run()
    {
        boolean bool = true;
        try {
            //System.out.println(commande[2]);
            commande[2] = commande[2].replaceAll("\\)", "\\\\)");
            commande[2] = commande[2].replaceAll("\\(", "\\\\(");

            r = Runtime.getRuntime();
            p = r.exec(commande, envp, new java.io.File(pathExecution));
            p.waitFor();
        } catch (Exception localException) {
            System.out.println("Erreur Execution : ");
            localException.printStackTrace();
            bool = false;
        }
        return bool;
    }

    public void stop()
    {
        if (p != null) {
            close(p.getOutputStream());
            close(p.getInputStream());
            close(p.getErrorStream());
            p.destroy();
        }
    }

    private static void close(Closeable paramCloseable)
    {
        if (paramCloseable != null) {
            try {
                paramCloseable.close();
            } catch (java.io.IOException localIOException) {
            }
        }
    }

    public String getInput()
    {
        StringBuffer localStringBuffer = new StringBuffer();
        String str = "";
        try {
            BufferedReader localBufferedReader = new BufferedReader(new java.io.InputStreamReader(p.getInputStream()));
            while ((str = localBufferedReader.readLine()) != null) {
                localStringBuffer.append(str + "\n");
            }
            localBufferedReader.close();
        } catch (Exception localException) {
            localException.printStackTrace();
        }
        return localStringBuffer.toString();
    }
    

    public void getError()
    {
        try {
            BufferedReader localBufferedReader = new BufferedReader(new java.io.InputStreamReader(p.getErrorStream()));
            while (localBufferedReader.readLine() != null) {
                 System.out.println(localBufferedReader.readLine());
            }

            localBufferedReader.close();
        } catch (Exception localException) {
            localException.printStackTrace();
        }
    }
}
