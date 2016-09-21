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
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;



public class Main {

	/**
	 * @author Alexander Herbig
	 * 
	 */
	public static void main(String[] args) throws Exception
	{
		String inFile;
		String refFile;
		String outFileDraft;
		String outFileRefMod;
		String outFileNR1234;
		double minQual;
		int minCov;
		double minSNPallelFreq;
		String draftName = "new draft sequence";
		
		
		String programName = "VCF2Genome";
		String author = "Alexander Herbig";
		String authorsEmail = "alexander.herbig@uni-tuebingen.de";
		String version = "0.84 (2015-03-09)";
		
		System.out.println(programName+" - "+version+"\nby "+author+"\n");
				
		
		if(args.length==0 || args[0].equalsIgnoreCase("--help") || args[0].equalsIgnoreCase("-help") || args[0].equalsIgnoreCase("-?") || args[0].equalsIgnoreCase("-h"))
		{
			System.out.println("--- How to use "+programName+" ---\n");
			printUsage();
			System.out.println("\nIn case of any questions contact "+author+" ("+authorsEmail+").");
			return;
		}
		
		//paras
		inFile = args[0];
		refFile = args[1];
		outFileDraft = args[2];
		outFileRefMod = args[3];
		outFileNR1234 = args[4];
		minQual = Double.parseDouble(args[5]);
		minCov = Integer.parseInt(args[6]);
		minSNPallelFreq = Double.parseDouble(args[7]);
		
		if(args.length>=9)
			draftName=args[8];
		
	
		double minHomSNPallelFreq = minSNPallelFreq;
		double minHetSNPallelFreq = minSNPallelFreq;

		
		String refGenome = FASTAParser.parseDNA(refFile).values().iterator().next();
			
		
		//SNP array
		char[] calls = new char[refGenome.length()];
		char[] uncertainCalls = new char[refGenome.length()];
		
		Set<Integer> snpPositions = new HashSet<Integer>();
		
		boolean[] missingDataPos = new boolean[refGenome.length()];
		
		
		//Fill SNP array
		for(int i=0; i<calls.length; i++)
		{
			calls[i] = refGenome.charAt(i);
			uncertainCalls[i] = refGenome.charAt(i);
		}
			
		
		
		////////////////////////
		// BEGIN -- Parse VCF //
		////////////////////////
		BufferedReader br;
		
		System.out.println("SNP statistics\nQuality Threshold: "+minQual+"\nCoverage Threshold: "+minCov+"\nMinimum SNP allele frequency: "+minHomSNPallelFreq);
		System.out.println("sample\tSNP Calls\tcoverage(fold)\tcoverage(percent)\trefCall\tallPos\tnoCall\tdiscardedRefCall\tdiscardedVarCall\tfilteredVarCall\tunhandledGenotype");
		
		char nChar='N';
		char rChar='R';
		//counter
		int covCount = 0;
		int allPos = 0;
		int noCallPos = 0;
		int nonStandardRefChars = 0;
		int refCallPos = 0;
		int varCallPos = 0;
		int discardedRefCall = 0;
		int discardedVarCall = 0;
		int filteredVarCall = 0;
		int unknownCall = 0;
		
		int ns;
		double nperc;
		double covround;
		
		String line;
		String[] cols;
		
		double qual;
		int cov;
		double SNPallelFreq;
		String[] allelCols;
		
		int lastPos1based = 0;
		int currPos1based = 0;
		
		
			br = new BufferedReader(new FileReader(inFile));
						
			
			//counter
			covCount=0;
			allPos = 0;
			noCallPos = 0;
			nonStandardRefChars = 0;
			refCallPos = 0;
			varCallPos = 0;
			discardedRefCall = 0;
			discardedVarCall = 0;
			filteredVarCall = 0;
			unknownCall = 0;
			
			lastPos1based = 0;
			currPos1based = 0;
			
			while((line=br.readLine()) != null)
			{
				if(line.startsWith("#"))
				{
					continue;
				}
				
				//count
				allPos++;
				
				
				cols = line.split("\t");
				
				//insert Ns at not handled sites, which are left out by GATK
				lastPos1based=currPos1based;
			
				currPos1based = Integer.parseInt(cols[1]);
				
				
				if(currPos1based-lastPos1based!=1)
				{
					if(lastPos1based>=currPos1based)
						throw new Error("ERROR: Base calls in the vcf file are not sorted! (Note that multiple chromosomes are not supported.)");
					
					for(int i=lastPos1based+1; i<currPos1based; i++)
					{
						allPos++;
						nonStandardRefChars++;
						
						calls[currPos1based-1] = nChar;
						uncertainCalls[currPos1based-1] = nChar;
						missingDataPos[currPos1based-1] = true;
						
					}
				}
				
				allelCols = cols[9].split(":");
				
				//No call
				if(allelCols[0].equals("./."))
				{
					noCallPos++;
					
					calls[currPos1based-1] = nChar;
					uncertainCalls[currPos1based-1] = nChar;
					missingDataPos[currPos1based-1] = true;
					
				}
				//Reference Call
				else if(allelCols[0].equals("0/0"))
				{
					qual = Double.parseDouble(cols[5]);
					cov=Integer.parseInt(allelCols[1]);
					covCount+=cov;
					
					if(qual>=minQual && cov >= minCov)
					{
						refCallPos++;
						
						// do nothing since reference is called
					}
					else
					{
						discardedRefCall++;
						
						calls[currPos1based-1] = nChar;
						uncertainCalls[currPos1based-1] = rChar;
						missingDataPos[currPos1based-1] = true;
						
					}
				}
				//variant call
				else if(allelCols[0].equals("0/1") || allelCols[0].equals("1/1"))
				{
					qual = Double.parseDouble(cols[5]);
					cov=Integer.parseInt(allelCols[1].split(",")[1]);
					SNPallelFreq=Math.min((double)cov/(cov+Integer.parseInt(allelCols[1].split(",")[0])-1) , 1); // -1 because once doesn't count
					
					covCount+=cov+Integer.parseInt(allelCols[1].split(",")[0]);
					
					if(qual>=minQual && cov >= minCov && SNPallelFreq >= minHomSNPallelFreq)
					{

						varCallPos++;
						
						calls[currPos1based-1] = cols[4].charAt(0);
						uncertainCalls[currPos1based-1] = cols[4].charAt(0);
						
						snpPositions.add(currPos1based);
						
					}
					else if(qual>=minQual && cov >= minCov && SNPallelFreq >= minHetSNPallelFreq)
					{

						varCallPos++;
						
						calls[currPos1based-1] = getAmbiguousBase(calls[currPos1based-1] , cols[4].charAt(0));
						uncertainCalls[currPos1based-1] = getAmbiguousBase(calls[currPos1based-1] , cols[4].charAt(0));

						snpPositions.add(currPos1based);
						
					}
					else
					{
						cov=Integer.parseInt(allelCols[1].split(",")[0]);
						
						SNPallelFreq=(double)cov/(cov+Integer.parseInt(allelCols[1].split(",")[1])-1); // -1 because once doesn't count
						
						if(qual>=minQual && cov >= minCov)
						{
							if(SNPallelFreq>=minHomSNPallelFreq)
							{
								refCallPos++;
							}
							else
							{
								discardedVarCall++;
								
								calls[currPos1based-1] = nChar;
								uncertainCalls[currPos1based-1] = nChar;

								missingDataPos[currPos1based-1] = true;
							}
							// do nothing since reference is called
						}
						else
						{
							discardedVarCall++;
							
							calls[currPos1based-1] = nChar;
							
							if(allelCols[0].equals("1/1"))
								uncertainCalls[currPos1based-1] = Character.toLowerCase(cols[4].charAt(0));
							else
								uncertainCalls[currPos1based-1] = nChar;

							missingDataPos[currPos1based-1] = true;
							
						}
					}
					
				}
				//unhandled genotype
				else
				{
					unknownCall++;
					
					calls[currPos1based-1] = nChar;
					uncertainCalls[currPos1based-1] = nChar;
					missingDataPos[currPos1based-1] = true;
				}
				
			}
			br.close();
			
			//write stats
			ns = discardedRefCall+discardedVarCall+noCallPos+unknownCall+nonStandardRefChars;
			nperc = Math.round(((ns*100d)/allPos)*100)/100d;
			covround = Math.round((covCount/(double)allPos)*100)/100d;

			System.out.println(getSampleNameFromPath(inFile)+"\t"+varCallPos+"\t"+covround+"\t"+(100-nperc)+"\t"+refCallPos+"\t"+allPos+"\t"+noCallPos+"\t"+discardedRefCall+"\t"+discardedVarCall+"\t"+filteredVarCall+"\t"+unknownCall);
		
		
		
		//////////////////////
		//END -- Parse VCFs //
		//////////////////////

		
		//Write draftN
		BufferedWriter bw = new BufferedWriter(new FileWriter(outFileDraft));
		
		
		StringBuffer tmpSeq;
		
		tmpSeq=new StringBuffer();
			for(int pos=1; pos <= calls.length; pos++)
				tmpSeq.append(calls[pos-1]);
			
			FASTAWriter.write(bw, draftName+"_draftN", tmpSeq.toString());
		
		
		bw.close();
		
		//Write draftUncertain
		bw = new BufferedWriter(new FileWriter(outFileNR1234));
		
		
		tmpSeq = new StringBuffer();
		
		for(int pos=1; pos <= calls.length; pos++)
			tmpSeq.append(uncertainCalls[pos-1]);
			
		FASTAWriter.write(bw, draftName+"_draftUncertain", tmpSeq.toString());
		
		
		bw.close();
		
		//Write draftRefMod
		bw = new BufferedWriter(new FileWriter(outFileRefMod));
		
		
		tmpSeq = new StringBuffer();
		char tmpchar;
		
			
			for(int pos=1; pos <= calls.length; pos++)
			{
				tmpchar = calls[pos-1];
				if(tmpchar == nChar)
					tmpSeq.append(refGenome.charAt(pos-1));
				else
					tmpSeq.append(tmpchar);
			}
			
			FASTAWriter.write(bw, draftName+"_draftRefMod", tmpSeq.toString());
		
		
		bw.close();
	
}
	
	private static void printUsage()
	{
		System.out.println("USAGE: java -jar VCF2Genome.jar <in.vcf> <reference-genome.fasta> <draft.fasta> <refMod.fasta> <uncertain.fasta> <minimum quality score> <minimum coverage> <minimum SNP allel frequency> [draft sequence name]");
		System.out.println("EXAMPLE: java -jar VCF2Genome.jar in.vcf reference.fasta draft.fasta refMod.fasta uncertain.fasta 40 7 0.9 my_new_draft_genome");
		System.out.println("");
		System.out.println("VCF2Genome generates a draft genome sequence from a GATK vcf file. The file has to contain a call for each site (also non-variants).");
		System.out.println("VCF2Genome is only applicable to single reference sequences. Multiple chromosomes are not supported.");
		System.out.println("VCF2Genome is only recommended for bacterial genomes.");
		System.out.println("Example GATK command line:");
		System.out.println("GenomeAnalysisTK.jar -T UnifiedGenotyper -R <referenceGenome.fasta> --output_mode EMIT_ALL_SITES -o inputForVCF2Genome.vcf -I <input.bam>");
		System.out.println("GATK needs bam files with read groups. Read groups can be added using PicardTools:");
		System.out.println("AddOrReplaceReadGroups.jar I=In.bam O=Out.bam LB=NA PL=illumina PU=NA SM=NA VALIDATION_STRINGENCY=SILENT");
		System.out.println("");
		System.out.println("More details:");
		System.out.println("draft vs. refMod: 'draft' contains 'N's where no call can be made, 'refMod' contains reference bases in this case.");
		System.out.println("uncertain.fasta: More precise uncertainty encoding. N: Not covered or ambiguous. R: Low coverage but looks like Ref. a,c,t,g (lower case): Low coverage but looks like SNP.");
		System.out.println("minimum quality score: The score is given in the 6th column of the vcf file. Phred-scaled quality score for the call. High QUAL scores indicate high confidence calls. -10log_10 p(no variant)");
		System.out.println("minimum coverage: minimum number of reads confirming the call.");
		System.out.println("minimum SNP allel frequency: minimum fraction of reads containing the called nucleotide.");
	}
	
	
	public static boolean parseBoolean(String boo)
	{
		boolean res = false;
		if(boo.equalsIgnoreCase("true") || boo.equalsIgnoreCase("T") || boo.equalsIgnoreCase("1"))
			res=true;
		return res;
	}
	
	public static char getAmbiguousBase(char c1, char c2)
	{
		Set<Character> chars = new HashSet<Character>();
		chars.add(c1);
		chars.add(c2);
		
		if(c1==c2)
			return c1;
		
		if(chars.contains('G') && chars.contains('A'))
			return 'R';
		if(chars.contains('T') && chars.contains('C'))
			return 'Y';
		if(chars.contains('G') && chars.contains('T'))
			return 'K';
		if(chars.contains('A') && chars.contains('C'))
			return 'M';
		if(chars.contains('G') && chars.contains('C'))
			return 'S';
		if(chars.contains('A') && chars.contains('T'))
			return 'W';
		
		System.err.println("Illegal arguments in function getAmbiguousBase: "+c1+", "+c2);
		
		return 'N';
	}
	
	public static String getSampleNameFromPath(String sampleNameWithPath)
	{
		String[] res = sampleNameWithPath.split("/");
		
		if(res.length>=2)
			return res[res.length-2];
		else
			return sampleNameWithPath;
	}
	

}
