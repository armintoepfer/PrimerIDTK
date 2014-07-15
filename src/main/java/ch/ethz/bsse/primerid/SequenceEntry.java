/**
 * Copyright (c) 2014 Armin Töpfer
 *
 * This file is part of PrimerIDTK.
 *
 * PrimerIDTK is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * PrimerIDTK is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * PrimerIDTK. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.primerid;

import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class SequenceEntry {

    final int HAMMING_DISTANCE = 2;
    final String NSH = "GCCTTGCACG";
//    final String NSH = "GCCTTGCCAGCACGCTCAGGCCTTGCACG";
    public String sequence;
    public String header;
    public String quality;
    public int hamming_distance;
    public String primer_id;
    public SequenceEntry crick;

    public SequenceEntry(String sequence, String header, String quality) {
        this.sequence = sequence;
        this.header = header;
        this.quality = quality;
    }

    private boolean isGAP(char c) {
        return c == '-' || c == 'N';
    }

    public void split(String prefix) {
        Alignment alignNSH = SmithWatermanGotoh.align(
                new Sequence(this.sequence.substring(0, 100), "", "", Sequence.NUCLEIC),
                new Sequence(NSH, "", "", Sequence.NUCLEIC),
                Globals.MATRIX, 30, 3);
        if (alignNSH == null) {
            System.err.println("NO MATCH AT ALL");
            System.exit(0);
            return;
        }

        char[] nsh_cigar = computeCigar(alignNSH.getMarkupLine(), alignNSH.getSequence2(), alignNSH.getSequence1());
        int nsh_end = alignNSH.getStart1();
        int nsh_mismatches = 0;
        int nsh_hit = 0;
        int nsh_del = 0;
        int nsh_ins = 0;
        for (char h : nsh_cigar) {
            if (h != 0) {
                switch (h) {
                    case 'X':
                        nsh_hit++;
                        nsh_mismatches++;
                        nsh_end++;
                        break;
                    case 'D':
                        nsh_del++;
                        nsh_end++;
                        break;
                    case 'M':
                        nsh_hit++;
                        nsh_end++;
                        break;
                    case 'I':
                        nsh_ins++;
                        nsh_hit++;
                }
            }
        }
        if (nsh_hit < 10) {
            Globals.NSH_TOO_SMALL++;
            return;
        }
        int maxL = Math.min(nsh_end + prefix.length() + 20, sequence.length());
        Alignment alignPrefix = SmithWatermanGotoh.align(
                new Sequence(this.sequence.substring(nsh_end, maxL), "", "", Sequence.NUCLEIC),
                new Sequence(prefix, "", "", Sequence.NUCLEIC),
                Globals.MATRIX, 30, 3);
        char[] prefix_cigar = computeCigar(alignPrefix.getMarkupLine(), alignPrefix.getSequence2(), alignPrefix.getSequence1());
        int prefix_length = 0;
        int prefix_mismatches = 0;
        int prefix_del = 0;
        int prefix_ins = 0;
        for (char h : prefix_cigar) {
            if (h != 0) {
                switch (h) {
                    case 'X':
                        prefix_mismatches++;
                        prefix_length++;
                        break;
                    case 'D':
                        prefix_del++;
                        prefix_length++;
                        break;
                    case 'M':
                        prefix_length++;
                        break;
                    case 'I':
                        prefix_ins++;
                }
            }
        }
        int prefixStart = alignPrefix.getStart1();
        if (alignPrefix.getStart2() > 0) {
            Globals.RTP_NOT_FIRST++;
            return;
        }
        boolean failed = false;
        if (nsh_mismatches > 2 && nsh_del > 1 && nsh_ins > 1) {
            Globals.NSH_MISMATCH++;
            failed = true;
        }
        if (prefixStart >= 20) {
            Globals.PID_TOO_LARGE++;
            failed = true;
        }
        if (prefix_mismatches > 2 && prefix_del > 1 && prefix_ins > 1) {
            Globals.RTP_MISMATCH++;
            failed = true;
        }
        if (prefix_length < prefix.length() / 2d) {
            Globals.RTP_TOO_SMALL++;
            failed = true;
        }
        if (!failed) {
            this.primer_id = Utils.reverseComplement(sequence.substring(nsh_end, nsh_end + prefixStart));
        }
    }

    private char[] computeCigar(char[] m, char[] c, char[] g) {
        char[] cigar = new char[m.length];
        for (int j = 0; j < m.length; j++) {
            char currentConsensus = '*';
            try {
                if (m[j] == '|') {
                    currentConsensus = c[j];
                    cigar[j] = 'M';
                } else if (m[j] == ' ') {
                    if (isGAP(c[j]) && isGAP(g[j])) {
                        cigar[j] = 'D';
                        currentConsensus = '-';
                    } else if (isGAP(c[j])) {
                        currentConsensus = '-';
                        cigar[j] = 'D';
                    } else if (isGAP(g[j])) {
                        cigar[j] = 'I';
                        currentConsensus = c[j];
                    }
                } else if (m[j] == '.') {
                    if (isGAP(g[j])) {
                        currentConsensus = c[j];
                        if (currentConsensus == '-') {
                            cigar[j] = 'D';
                        } else {
                            cigar[j] = 'X';
                        }
                    } else if (isGAP(c[j])) {
                        cigar[j] = 'D';
                        currentConsensus = c[j];
                    } else {
                        cigar[j] = 'X';
                        currentConsensus = c[j];
                    }
                } else if (m[j] == ':') {
                    if (isGAP(c[j])) {
                        cigar[j] = 'D';
                        currentConsensus = '-';
                    } else {
                        cigar[j] = 'X';
                        currentConsensus = c[j];
                    }
                }

                if (c[j] == 'N') {
                    currentConsensus = c[j];
                    if (isGAP(g[j])) {
                        cigar[j] = 'I';
                    } else {
                        cigar[j] = 'M';
                    }
                }
            } catch (Exception e) {
                System.err.println("FUTURE: " + e);
            }
        }
        return cigar;
    }
}
