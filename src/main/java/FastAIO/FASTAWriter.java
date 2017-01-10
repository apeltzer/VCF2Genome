package FastAIO;
/*
 * Copyright (c) 2016. VCF2Genome Alexander Herbig, Alexander Peltzer
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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Writes a list of sequences to a file in FASTA format.
 *
 * @author Alexander Herbig
 */
public class FASTAWriter {
    /**
     * Writes a list of sequences to a file in FASTA format.
     *
     * @param bw           the BufferedWriter writing the FASTA file
     * @param fastaEntries list of FASTAEntry objects holding the sequences
     *                     which are written to the FASTA file
     * @throws Exception
     */

    private BufferedWriter bfw;
    private String genomeID;
    private String sequence;

    public FASTAWriter(String filePath, String genomeID, String sequence) throws Exception{
        this.bfw = new BufferedWriter(new FileWriter(new File(filePath)));
        this.genomeID = genomeID;
        this.sequence = sequence;
        writeFastAFile();
    }


    /**
     * Method that writes FastA entries into a file..
     *
     * @throws IOException
     */
    private void writeFastAFile() throws IOException {
        int charsInLine;
        String tmpSeqString;

        tmpSeqString = sequence;
        bfw.write(">" + genomeID);
        bfw.newLine();

        charsInLine = 0;
        for (int i = 0; i < tmpSeqString.length(); i++) {
            bfw.write(tmpSeqString.charAt(i));
            charsInLine++;

            if (charsInLine == 60) {
                bfw.newLine();
                charsInLine = 0;
            }
        }
        if (charsInLine != 0) {
            bfw.newLine();
        }

        bfw.flush();
        bfw.close();
    }
}
