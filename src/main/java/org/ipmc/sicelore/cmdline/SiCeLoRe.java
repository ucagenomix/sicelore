package org.ipmc.sicelore.cmdline;

import org.broadinstitute.barclay.argparser.CommandLineProgramGroup;

public class SiCeLoRe implements CommandLineProgramGroup {
  public SiCeLoRe() {}
   
  public String getName() { return "SiCeLoRe Pipeline"; }
  
  public String getDescription()
  {
    return "SiCeLoRe Pipeline - [Si]ngle [Ce]ll [Lo]ng [Re]ads tools";
  }
}

