package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public interface GeneModelParser {

    public TranscriptRecord parseLine(String line) throws GTFParseException;
}
