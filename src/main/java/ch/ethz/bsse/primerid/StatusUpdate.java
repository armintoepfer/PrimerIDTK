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

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.TimeZone;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class StatusUpdate {

    private static final long start = System.currentTimeMillis();
    private static final DateFormat df;
    public static int readCount = 0;
    private static String oldOut = "";

    static {
        df = new SimpleDateFormat("HH:mm:ss");
        df.setTimeZone(TimeZone.getTimeZone("GMT"));
    }

    public static void print(String text, int count) {
        if (!oldOut.equals(time())) {
            oldOut = time();
            System.out.print("\r" + time() + "  " + text + count + " ("+ (count-readCount) +" rpm)");
            readCount = count;
        }
    }

    public static void print(String text) {
        System.out.print("\r" + time() + "  " + text);
    }

    public synchronized static void processReads() {
        readCount++;
        if (!oldOut.equals(time())) {
            oldOut = time();
            System.out.print("\r                                                                                                                      ");
            System.out.print("\r" + time() + " Mapped: " + readCount);
        }
    }

    public static void println(String text) {
        System.out.println("\r" + time() + "  " + text);
    }

    private static String time() {
        return df.format(new Date(System.currentTimeMillis() - start));
    }
}
