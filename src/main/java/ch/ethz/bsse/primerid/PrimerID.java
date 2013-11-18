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

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.logging.Handler;
import java.util.logging.Logger;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class PrimerID {
    
    @Option(name = "-i")
    private String input;
    @Option(name = "-rt")
    private String rt_primer;
    @Option(name = "-s")
    private String saved;
    @Option(name = "-xlr")
    private int xlr = -1;
    
    {
        Logger rootLogger = Logger.getLogger("");
        Handler[] handlers = rootLogger.getHandlers();
        if (handlers.length > 0) {
            rootLogger.removeHandler(handlers[0]);
        }
    }
    
    public void doMain(String[] args) throws IOException {
        CmdLineParser parser = new CmdLineParser(this);
        
        parser.setUsageWidth(80);
        try {
            parser.parseArgument(args);
            this.setGlobal();
            if (saved == null) {
                new Preprocessing().compute();
                try {
                    FileOutputStream fos = new FileOutputStream(Globals.INPUT + ".saved");
                    try (ObjectOutputStream out = new ObjectOutputStream(fos)) {
                        out.writeObject(new InformationHolder(Globals.PID_SIZE));
                    }
                } catch (IOException ex) {
                }
            } else {
                InformationHolder ih = null;
                try {
                    FileInputStream fis = new FileInputStream(this.saved);
                    try (ObjectInputStream in = new ObjectInputStream(fis)) {
                        ih = (InformationHolder) in.readObject();
                    }
                } catch (IOException | ClassNotFoundException ex) {
                    System.err.println(ex);
                }
                Globals.PID_SIZE = ih.pid_size;
                ih = null;
                System.out.println(Globals.PID_SIZE.keySet());
            }
        } catch (CmdLineException cmderror) {
            System.err.println(cmderror.getMessage());
            System.err.println("");
            System.err.println("QuasiRecomb version: " + PrimerID.class.getPackage().getImplementationVersion());
            System.err.println("Get latest version from FILLME");
            System.err.println("");
            System.err.println("USAGE: java -jar PrimerIDTK.jar options...\n");
            System.err.println(" -------------------------");
            System.err.println(" === GENERAL options ===");
            System.err.println("  -i INPUT\t\t: Fastq input file.");
            System.err.println(" -------------------------");
            System.err.println(" === Technical options ===");
            System.err.println("  -XX:NewRatio=9\t: Reduces the memory consumption (RECOMMENDED to use).");
            System.err.println("  -Xms2G -Xmx10G\t: Increase heap space.");
            System.err.println("  -XX:+UseParallelGC\t: Enhances performance on multicore systems.");
            System.err.println("  -XX:+UseNUMA\t\t: Enhances performance on multi-CPU systems.");
            System.err.println(" -------------------------");
            System.err.println(" === EXAMPLES ===");
            System.err.println("   java -XX:+UseParallelGC -Xms2g -Xmx10g -XX:+UseNUMA -XX:NewRatio=9 -jar PrimerIDTK.jar -i R2.fastq");
            System.err.println(" -------------------------");
        }
    }
    
    private void setGlobal() {
        if (xlr != -1) {
            String[] files = new String[]{
                "/Users/XLR/Dropbox/Projects/PrimerID/raw_data/a.fastq",
                "3223a_S1_L001_R2_001.fastq",
                "3223b_S2_L001_R2_001.fastq",
                "3223c_S3_L001_R2_001.fastq",
                "3236a_S4_L001_R2_001.fastq",
                "3236b_S5_L001_R2_001.fastq",
                "3236c_S6_L001_R2_001.fastq"};
            String[] primer = new String[]{"GGTTCTTTCTGATG", "GGTTCTTTCTGATG", "GGTTCTTTCTGATG", "GGTTCTTTCTGATG", "CCAAAGGAATGGAGGTTCTTTCTGATG", "CCAAAGGAATGGAGGTTCTTTCTGATG", "CCAAAGGAATGGAGGTTCTTTCTGATG"};
            Globals.INPUT = files[xlr];
            Globals.RT_PRIMER = primer[xlr];
        } else {
            Globals.INPUT = this.input;
            Globals.RT_PRIMER = this.rt_primer;
        }
    }

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {
        new PrimerID().doMain(args);
        System.exit(0);
    }
}
