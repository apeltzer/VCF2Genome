package FastAIO;

/*
 * Copyright (c) 2016. VCF2Genome Alexander Herbig
 * This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;
import java.util.Map;

/**
 * Reads a multiple FASTA file containing DNA sequences in FASTA format
 * and stores the entries as FASTAEntry objects.
 *
 * @author Alexander Herbig
 */
public class FASTAParser {
    /**
     * Reads a multiple FASTA file containing DNA sequences in FASTA format.
     *
     * @param br the BufferedReader which reads from the multiple FASTA file
     * @throws Exception
     * @return a list containing the resulting FASTAEntry objects
     */
    public static Map<String, String> parseDNA(String filename) throws Exception {
        BufferedReader br = new BufferedReader(new FileReader(filename));

        Map<String, String> fastaEntries = new HashMap<String, String>();

        String tmpID;
        StringBuffer tmpSeqString = new StringBuffer();
        String tmpLine = br.readLine();

        //jump to first header line
        while (tmpLine != null && (tmpLine.charAt(0) != '>' || tmpLine.length() == 0)) {
            tmpLine = br.readLine();
        }
        tmpID = tmpLine;

        tmpLine = br.readLine();
        while (tmpLine != null) {
            if (tmpLine.length() != 0) {
                if (tmpLine.charAt(0) == '>') {
                    fastaEntries.put(toID(tmpID), tmpSeqString.toString());

                    tmpSeqString = new StringBuffer();

                    tmpID = tmpLine;
                } else {
                    tmpSeqString.append(tmpLine);
                }
            }
            tmpLine = br.readLine();
        }
        fastaEntries.put(toID(tmpID), tmpSeqString.toString());

        br.close();

        return fastaEntries;
    }

    private static String toID(String fastaID) {
        return (fastaID.substring(1));
    }

}
