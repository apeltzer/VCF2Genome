
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

import java.io.BufferedWriter;

/**
 * Writes a list of sequences to a file in FASTA format.
 * 
 * @author Alexander Herbig
 *
 *
 */
public class FASTAWriter 
{
	/**
	 * Writes a list of sequences to a file in FASTA format.
	 * 
	 * @param bw	the BufferedWriter writing the FASTA file
	 * @param fastaEntries	list of FASTAEntry objects holding the sequences 
	 * 						which are written to the FASTA file
	 * @throws Exception
	 */
	public static void write(BufferedWriter bw, String genomeID, String sequence) throws Exception
	{	
		int charsInLine;
		String tmpSeqString;
		
		tmpSeqString=sequence;
		bw.write(">"+genomeID);
		bw.newLine();
		
		charsInLine=0;
		for(int i=0; i<tmpSeqString.length();i++)
		{
			bw.write(tmpSeqString.charAt(i));
			charsInLine++;
			
			if(charsInLine==60)
			{
				bw.newLine();
				charsInLine=0;
			}
		}
		if(charsInLine!=0)
		{
			bw.newLine();
		}
		
		bw.flush();
	}
}
