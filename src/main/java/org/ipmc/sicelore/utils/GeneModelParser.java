package org.ipmc.sicelore.utils;

public interface GeneModelParser {
    public TranscriptRecord parseLine(String line) throws GTFParseException;
}
