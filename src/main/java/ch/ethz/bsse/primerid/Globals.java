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

import jaligner.matrix.Matrix;
import jaligner.matrix.MatrixLoader;
import jaligner.matrix.MatrixLoaderException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class Globals {

    public static String INPUT;
    public static String INPUT2;
    public static String RT_PRIMER;

    public static int NSH_TOO_SMALL = 0;
    public static int NSH_MISMATCH = 0;
    public static int RTP_TOO_SMALL = 0;
    public static int RTP_MISMATCH = 0;
    public static int RTP_NOT_FIRST = 0;
    public static int PID_TOO_LARGE = 0;
    public static int KMER_LENGTH = 6;

    public static Map<String, Integer> PID_OCCURENCE = new HashMap<>();
    public static Map<Integer, Integer> PID_SIZE = new HashMap<>();
    
    public static Matrix MATRIX = loadMatrix();

    private static Matrix loadMatrix() {
        Matrix m = null;
        try {
            m = MatrixLoader.load("EDNAFULL");
        } catch (MatrixLoaderException ex) {
            System.err.println("Matrix error loading");
        }
        return m;
    }
}
