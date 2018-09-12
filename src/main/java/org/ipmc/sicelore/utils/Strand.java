package org.ipmc.sicelore.utils;

/**
 * 
 * @author kevin lebrigand
 * 
 */
public enum Strand {
    FORWARD {
        public String toString() {
            return "+";
        }
    },
    REVERSE {
        public String toString() {
            return "-";
        }
    };

    public static Strand fromString(String str) throws GTFParseException {
        if (str.equalsIgnoreCase(Strand.FORWARD.toString())) {
            return Strand.FORWARD;
        } else if (str.equalsIgnoreCase(Strand.REVERSE.toString())) {
            return Strand.REVERSE;
        } else {
            throw new GTFParseException("Invalid strand '" + str + "'");
        }
    }
}
