/**
 * Copyright (c) 2013 Armin Töpfer
 *
 * This file is part of PrimerIDTK.
 *
 * InDelFixer is free software: you can redistribute it and/or modify it under
 * the terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or any later version.
 *
 * InDelFixer is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * InDelFixer. If not, see <http://www.gnu.org/licenses/>.
 */
package ch.ethz.bsse.primerid;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Preprocessing {

    public void compute() throws FileNotFoundException, IOException {
        BufferedReader brWatson = new BufferedReader(new FileReader(new File(Globals.INPUT)));
        List<SequenceEntry> sequences = new LinkedList<>();
        int j = 0;
        for (;;) {
            j++;
            try {
                SequenceEntry watsonQ = parseFastq(brWatson);
                if (watsonQ != null) {
                    sequences.add(watsonQ);
                    watsonQ.split(Globals.RT_PRIMER);
                    StatusUpdate.print("Reads:\t ", j);
                }
            } catch (IllegalAccessError e) {
                break;
            }
        }
        StatusUpdate.println("Reads:\t " + j);
        Map<String, List<SequenceEntry>> primer_sequences = new HashMap<>();
        for (SequenceEntry sequence : sequences) {
            if (!Globals.PID_OCCURENCE.containsKey(sequence.primer_id)) {
                Globals.PID_OCCURENCE.put(sequence.primer_id, 1);
            } else {
                Globals.PID_OCCURENCE.put(sequence.primer_id, Globals.PID_OCCURENCE.get(sequence.primer_id) + 1);
            }
            if (!primer_sequences.containsKey(sequence.primer_id)) {
                primer_sequences.put(sequence.primer_id, new LinkedList<SequenceEntry>());
            }
            primer_sequences.get(sequence.primer_id).add(sequence);
        }
        StringBuilder out = new StringBuilder(Globals.PID_OCCURENCE.size() * 10);
        int readCount = 0;
        int uniques = 0;
        for (Map.Entry<String, Integer> e : Globals.PID_OCCURENCE.entrySet()) {
            if (e.getKey() != null) {
                readCount += e.getValue();
                uniques++;
                out.append(e.getKey()).append("\t").append(e.getValue()).append("\n");
                if (!Globals.PID_SIZE.containsKey(e.getKey().length())) {
                    Globals.PID_SIZE.put(e.getKey().length(), 0);
                }
                Globals.PID_SIZE.put(e.getKey().length(), Globals.PID_SIZE.get(e.getKey().length()) + 1);
            }
        }
        Utils.saveFile(Globals.INPUT + ".dist", out.toString());
        StringBuilder stats = new StringBuilder();
        stats.append("Reads with primerID:\t ").append(readCount).append("\n");
        stats.append("Unique primerIDs:\t ").append(uniques).append("\n");
        stats.append("NSH too short:\t\t ").append(Globals.NSH_TOO_SMALL).append("\n");
        stats.append("NSH did not match:\t ").append(Globals.NSH_MISMATCH).append("\n");
        stats.append("RT Primer too short:\t ").append(Globals.RTP_TOO_SMALL).append("\n");
        stats.append("RT Primer did not match: ").append(Globals.RTP_MISMATCH).append("\n");
        stats.append("RT Primer soft clipped:\t ").append(Globals.RTP_NOT_FIRST).append("\n");
        stats.append("Primer ID too long:\t ").append(Globals.PID_TOO_LARGE).append("\n");
        Utils.saveFile(Globals.INPUT + ".stats", stats.toString());

        System.out.println(" === STATS ===");
        System.out.println(stats.toString());
//            new File(files[i] + "_sequences/").mkdirs();
//            for (Map.Entry<String, List<SequenceEntry>> e : primer_sequences.entrySet()) {
//                StringBuilder sb = new StringBuilder();
//                if (e.getValue().size() >= 3) {
//                    int x = 0;
//                    for (SequenceEntry se : e.getValue()) {
//                        sb.append("@").append(x++).append("\n");
//                        sb.append(se.sequence).append("\n");
//                        sb.append("+").append("\n");
//                        sb.append(se.quality).append("\n");
//                    }
//                    saveFile(files[i] + "_sequences/" + e.getValue().size() + "_" + e.getKey() + ".fastq", sb.toString());
//                }
//            }
        //            saveFile(files[i] + "_" + k + ".dist", out.toString());
//        }
//        }
        StringBuilder lengthDistribution = new StringBuilder();
        for (Map.Entry<Integer, Integer> e : Globals.PID_SIZE.entrySet()) {
            lengthDistribution.append(e.getKey()).append("\t").append(e.getValue()).append("\n");
        }
        Utils.saveFile(Globals.INPUT + ".lengthDistribution", lengthDistribution.toString());
    }

    /**
     * Parses a single fastq block of given BufferedReader.
     *
     * @param br BufferedReader
     * @return Single fastq block of type SequenceEntry
     * @throws IllegalAccessError Thrown if BufferedReader has reached EOF.
     */
    private static SequenceEntry parseFastq(BufferedReader br) throws IOException, IllegalAccessError {
        //head
        String header = br.readLine();
        if (header == null) {
            throw new IllegalAccessError();
        }
        //sequence
        String seq = br.readLine();
        char[] c = seq.toCharArray();
        //description
        br.readLine();
        //quality
        String qualityString = br.readLine();
        char[] quality = qualityString.toCharArray();
        int[] p = new int[quality.length];
        for (int i = 0; i < quality.length; i++) {
            p[i] = ((int) quality[i]) - 33;
        }

        return new SequenceEntry(seq,
                header,
                qualityString);
    }

}
