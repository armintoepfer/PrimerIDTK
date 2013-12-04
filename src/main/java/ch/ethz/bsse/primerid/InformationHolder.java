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

import java.io.Serializable;
import java.util.Map;

/**
 * @author Armin Töpfer (armin.toepfer [at] gmail.com)
 */
public class InformationHolder implements Serializable {

    private static final long serialVersionUID = 12L;
    public Map<Integer, Integer> pid_size;
    public Map<String, Integer> pid_occurence;

    public InformationHolder(Map<Integer, Integer> pid_size, Map<String, Integer> pid_occurence) {
        this.pid_size = pid_size;
        this.pid_occurence = pid_occurence;
    }

}
