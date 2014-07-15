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

import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 *
 * @author toepfera
 */
public class KmerStats {

    private Map<Integer, Map<String, Integer>> kmerMap = new HashMap<>();

    public void splitPIDs() {
        for (Map.Entry<String, Integer> e : Globals.PID_OCCURENCE.entrySet()) {
            if (e.getKey() == null || e.getKey().isEmpty()) {
                continue;
            }
            int c = 0;
            if (e.getValue() < 100) {
                c = 0;
//            } else if (e.getValue() >= 10 && e.getValue() < 20) {
//                c = 10;
//            } else if (e.getValue() >= 20 && e.getValue() < 30) {
//                c = 20;
//            } else if (e.getValue() >= 30 && e.getValue() < 40) {
//                c = 30;
//            } else if (e.getValue() >= 50 && e.getValue() < 100) {
//                c = 50;
            } else if (e.getValue() >= 100) {
                c = 100;
            }
            splitSingle(e, c);
        }
        Map<String, Double> comparisonMap = new HashMap<>();
        boolean first = true;
        StringBuilder zSB = new StringBuilder();
        for (Map.Entry<Integer, Map<String, Integer>> e : kmerMap.entrySet()) {
            DescriptiveStatistics ds = new DescriptiveStatistics();
            for (Map.Entry<String, Integer> e2 : e.getValue().entrySet()) {
                ds.addValue(e2.getValue());
            }
            System.out.println(e.getKey() + ": " + ds.getMean() + " " + ds.getStandardDeviation());
            double mean = ds.getMean();
            double sd = ds.getStandardDeviation();
            ValueComparator bvc = new ValueComparator(e.getValue());
            TreeMap<String, Integer> sorted_map = new TreeMap<>(bvc);
            sorted_map.putAll(e.getValue());
            for (String k : sorted_map.descendingKeySet()) {
                double z = (e.getValue().get(k) - mean) / sd;
                if (first) {
                    comparisonMap.put(k, z);
                } else {
                    if (z >= 2 || z <= -2) {
                        zSB.append(k).append("\t").append(Math.round(z*100)/100d).append("\t").append(Math.round(comparisonMap.get(k)*100)/100d).append("\n");
                    }
//                    System.out.println(k + " " + Math.round(z) + " " + Math.round(comparisonMap.get(k)));
//                    if (z > 2 || z < -2) {
//                        if (z - comparisonMap.get(k) >= 2)
//                        System.out.println(k + " " + Math.round(z) + " " + Math.round(comparisonMap.get(k)));
//                    }
                }
//                if (z > 1) {
//                    System.out.println(k + " " + Math.round(z));
//                }
            }
            first = false;
        }
        System.out.println(zSB.toString());
        Utils.saveFile(Globals.INPUT + ".z", zSB.toString());
    }

    private void splitSingle(Map.Entry<String, Integer> e, int c) {
        if (!kmerMap.containsKey(c)) {
            kmerMap.put(c, new HashMap<String, Integer>());
        }
        Map<String, Integer> map = kmerMap.get(c);
        for (int i = Globals.KMER_LENGTH; i <= e.getKey().length(); i++) {
            String kmer = e.getKey().substring(i - Globals.KMER_LENGTH, i);
            if (!map.containsKey(kmer)) {
                map.put(kmer, 0);
            }
            map.put(kmer, map.get(kmer) + 1);
        }
    }
}

class ValueComparator implements Comparator<String> {

    Map<String, Integer> base;

    public ValueComparator(Map<String, Integer> base) {
        this.base = base;
    }

    // Note: this comparator imposes orderings that are inconsistent with equals.    
    public int compare(String a, String b) {
        if (base.get(a) >= base.get(b)) {
            return 1;
        } else {
            return -1;
        } // returning 0 would merge keys
    }
}
