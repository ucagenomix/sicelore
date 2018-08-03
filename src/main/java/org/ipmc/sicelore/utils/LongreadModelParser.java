package org.ipmc.sicelore.utils;

import htsjdk.samtools.SAMRecord;

public interface LongreadModelParser {
    public LongreadRecord parseSAMRecord(SAMRecord r) throws LongreadParseException;
}