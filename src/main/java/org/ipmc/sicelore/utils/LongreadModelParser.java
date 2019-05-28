package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
import htsjdk.samtools.SAMRecord;

public interface LongreadModelParser {

    public LongreadRecord parseSAMRecord(SAMRecord r, boolean primary_only) throws LongreadParseException;
}
