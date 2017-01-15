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

import FastAIO.FASTAParser;
import FastAIO.FASTAWriter;
import com.google.common.io.Files;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.kohsuke.args4j.*;
import org.kohsuke.args4j.spi.BooleanOptionHandler;
import org.kohsuke.args4j.spi.StringOptionHandler;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;


public class VCF2Genome {
    /**
     * Important program information
     */

    public static final String version = "0.91";
    public static final String date = "2017-01-13";
    public static final String programName = "VCF2Genome (v. "+ version + " " + date + ")";
    public static final String author = "Alexander Herbig (<= v0.84) and Alexander Peltzer (>v0.84)";
    public static final String authorsEmail = "herbig@shh.mpg.de, peltzer@shh.mpg.de";
    public static final String programInfo = programName + "\nby " + author + "\n" + authorsEmail + "\n";


    /**
     * Option handling with args4j, to keep things clean in the code.
     */

    @Option(name="-in", handler= StringOptionHandler.class, required = true, usage="input VCF file")
    private String inFile;

    @Option(name="-ref", handler= StringOptionHandler.class, required = true, usage="reference genome in FastA format")
    private String refFile;

    @Option(name="-draft", handler= StringOptionHandler.class, required = true, usage="draft contains Ns where no call can be made. RefMod contains reference calls instead at these positions.")
    private String outFileDraft;

    @Option(name="-refMod", handler= StringOptionHandler.class, required = true, usage="More precise uncertainty encoding. N: Not covered or ambiguous. R: Low coverage but looks like Ref. a,c,t,g (lower case): Low coverage but looks like SNP.")
    private String outFileRefMod;

    @Option(name="-uncertain", handler= StringOptionHandler.class, required = true, usage="Special 1234 encoded FastA output.")
    private  String outFileNR1234;

    @Option(name="-minq", required = true, usage="Minimum quality score. For UG: Phred scaled quality score. For HC genome quality score.", metaVar="MIN_QUAL_SCORE")
    private double minQual;

    @Option(name="-minc", required = true, usage="Minimum coverage / reads confirming the call.", metaVar="MIN_COVERAGE_FOR_SNP")
    private int minCov;

    @Option(name="-minfreq", required = true, usage="Minimum fraction of reads supporting the called nucleotide.", metaVar="MIN_SNP_FREQUENCY")
    private  double minSNPalleleFreq;

    @Option(name="-draftname", required = true, usage="Name of the draft sequence.", metaVar="DRAFT_SEQ_NAME")
    private  String draftName = "Test";

    @Option(name="-h", usage="Display this help information and exit.", required=false, handler= BooleanOptionHandler.class)
    private boolean help = false;

    @Argument
    private boolean hadArguments = false;

    /**
     * File Variables etc
     */


    /**
     * Main method, used for calling class method
     * @param args
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        VCF2Genome vcf2genome = new VCF2Genome(args);
    }


    /**
     * @author Alexander Herbig & Alexander Peltzer
     */
    public VCF2Genome(String[] args) throws Exception {
        //State basic program Information on CLI
        System.out.println(programInfo);

        //Parse our CLI part
        doMain(args);

        if(!hadArguments){
            return;
        } else {
            runAnalysis();
        }

    }

    private void runAnalysis() throws Exception {

        double minHomSNPallelFreq = minSNPalleleFreq;
        double minHetSNPallelFreq = minSNPalleleFreq;


        FASTAParser fastaparser = new FASTAParser(refFile);
        String refGenome = fastaparser.getFastaFile().values().iterator().next();

        //SNP array
        String[] calls = new String[refGenome.length()];
        String[] uncertainCalls = new String[refGenome.length()];

        Set<Integer> snpPositions = new HashSet<Integer>();

        String[] missingDataPos = new String[refGenome.length()]; //TODO this is required


        //Fill SNP array
        for (int i = 0; i < calls.length; i++) {
            calls[i] = String.valueOf(refGenome.charAt(i));
            uncertainCalls[i] = String.valueOf(refGenome.charAt(i));
        }


        ////////////////////////
        // BEGIN -- Parse VCF //
        ////////////////////////
        BufferedReader br;

        System.out.println("SNP statistics\nQuality Threshold: " + minQual + "\nCoverage Threshold: " + minCov + "\nMinimum SNP allele frequency: " + minHomSNPallelFreq);
        System.out.println("sample\tSNP Calls\tcoverage(fold)\tcoverage(percent)\trefCall\tallPos\tnoCall\tdiscardedRefCall\tdiscardedVarCall\tfilteredVarCall\tunhandledGenotype");

        char nChar = 'N';
        char rChar = 'R';
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

        double qual;
        int cov;
        double SNPallelFreq;

        int lastPos1based = 0;
        int currPos1based = 0;


        //counter
        covCount = 0;
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

        VCFFileReader vcfFileReader = new VCFFileReader(new File(inFile), false);
        Iterator iter = vcfFileReader.iterator();

        while (iter.hasNext()) {
            allPos++;
            VariantContext variantContext = (VariantContext) iter.next();

            lastPos1based = currPos1based;

            currPos1based = variantContext.getStart();

            //Introduce Ns for non-covered sites
            if (currPos1based - lastPos1based != 1) {
                if (lastPos1based >= currPos1based) {
                    throw new IOException("Error: Base calls in the vcf file are not sorted! (Note that we currently don't support multiple chromosomes, too!");
                }

                for (int i = lastPos1based + 1; i < currPos1based; i++) {
                    allPos++;
                    nonStandardRefChars++;

                    calls[currPos1based - 1] = String.valueOf(nChar);
                    uncertainCalls[currPos1based - 1] = String.valueOf(nChar);
                    missingDataPos[currPos1based - 1] = "t";
                }
            }


            //Check what we have here!!!

            boolean isNoCall = false;
            boolean isRefCall = false;
            boolean isVarCall = false;
            boolean isVarHomCall = false;
            boolean isUnhandledCall = false;

            // HET HOM_REF HOM_VAR MIXED NO_CALL UNAVAILABLE

            for (Genotype g : variantContext.getGenotypes()) {
                variantContext.getAttribute("AD");
                GenotypeType gtt = g.getType();
                if (gtt == GenotypeType.NO_CALL) { //TODO check with test cases here!!!
                    isNoCall = true;
                }
                if (gtt == GenotypeType.HOM_REF) { //TODO check with test  cases here!!!
                    isRefCall = true;
                }
                if (gtt == GenotypeType.HET) { //TODO check with test cases here, too!
                    isVarCall = true;
                }
                if (gtt == GenotypeType.HOM_VAR) {
                    isVarHomCall = true; //TODO check with 1/1 string for example!!!
                }
                if(variantContext.getAlleles().size() > 2){ //We currently handle haploid calls only!
                    isUnhandledCall = true;
                }
            }


            //No call

            if (isNoCall & !isUnhandledCall) {
                noCallPos++;

                calls[currPos1based - 1] = String.valueOf(nChar);
                uncertainCalls[currPos1based - 1] = String.valueOf(rChar);
                missingDataPos[currPos1based - 1] = "t";

            }


            //Ref Call

            else if (isRefCall & !isUnhandledCall) {
                qual = variantContext.getPhredScaledQual();

                Genotype g = variantContext.getGenotype(0);
                cov = g.getDP();
                //cov = Integer.parseInt((String) variantContext.getAttribute("DP")); //This can be wrong!! (109 vs 4!)

                covCount += cov;

                if (qual >= minQual && cov >= minCov) {
                    refCallPos++;

                    // do nothing since reference is called
                } else {
                    discardedRefCall++;

                    calls[currPos1based - 1] = String.valueOf(nChar);
                    uncertainCalls[currPos1based - 1] = String.valueOf(rChar);
                    missingDataPos[currPos1based - 1] = "t";

                }
            }


            //VariantCall
            //TODO what to do with something like 1/2 - cant handle this with this kind of code right now!
            else if ((isVarCall | isVarHomCall) & !isUnhandledCall)  {
                qual = variantContext.getPhredScaledQual();
                //ned to get the second allele here, non ref
                Allele ref_allele = variantContext.getReference();
                Allele call_allele = variantContext.getAlternateAlleles().get(0);



                //Get coverage of alternative allele
                int index_call = variantContext.getAlleleIndex(call_allele);
                int cov_call = variantContext.getGenotype(0).getAD()[index_call];

                //Get coverage of reference allele
                int index_ref = variantContext.getAlleleIndex(ref_allele);
                int cov_ref = variantContext.getGenotype(0).getAD()[index_ref];

                int full_cov = cov_call + cov_ref;

                //Calculate SNPalellefrequency properly...

                SNPallelFreq = Math.min((double) cov_call / (full_cov - 1), 1); // -1 because once doesn't count

                covCount += full_cov;

                //Check homozygous case

                if (qual >= minQual && cov_call >= minCov && SNPallelFreq >= minHomSNPallelFreq) {

                    varCallPos++;

                    calls[currPos1based - 1] = String.valueOf(call_allele.getDisplayString().charAt(0));
                    uncertainCalls[currPos1based - 1] = String.valueOf(call_allele.getDisplayString().charAt(0));

                    snpPositions.add(currPos1based);
                } else if (qual >= minQual && cov_call >= minCov && SNPallelFreq >= minHetSNPallelFreq) {

                    varCallPos++;

                    calls[currPos1based - 1] = String.valueOf(getAmbiguousBase(calls[currPos1based - 1].charAt(0), call_allele.getDisplayString().charAt(0)));
                    uncertainCalls[currPos1based - 1] = String.valueOf(getAmbiguousBase(calls[currPos1based - 1].charAt(0), call_allele.getDisplayString().charAt(0)));

                    snpPositions.add(currPos1based);

                } else {

                    SNPallelFreq = (double) cov_ref / (full_cov - 1); // -1 because once doesn't count

                    if (qual >= minQual && cov_ref >= minCov) {
                        if (SNPallelFreq >= minHomSNPallelFreq) {
                            refCallPos++;
                        } else {
                            discardedVarCall++;

                            calls[currPos1based - 1] = String.valueOf(nChar);
                            uncertainCalls[currPos1based - 1] = String.valueOf(nChar);

                            missingDataPos[currPos1based - 1] = "t";
                        }
                        // do nothing since reference is called
                    } else {
                        discardedVarCall++;
                        call_allele = variantContext.getAlternateAlleles().get(0);
                        calls[currPos1based - 1] = String.valueOf(nChar);
                        uncertainCalls[currPos1based - 1] = String.valueOf(Character.toLowerCase(call_allele.getDisplayString().charAt(0)));
                        missingDataPos[currPos1based - 1] = "t";
                    }
                }
            } else {
                unknownCall++;

                calls[currPos1based - 1] = String.valueOf(nChar);
                uncertainCalls[currPos1based - 1] = String.valueOf(nChar);
                missingDataPos[currPos1based - 1] = "t";
            }
        }


        //write stats
        ns = discardedRefCall + discardedVarCall + noCallPos + unknownCall + nonStandardRefChars;
        nperc = Math.round(((ns * 100d) / allPos) * 100) / 100d;
        covround = Math.round((covCount / (double) allPos) * 100) / 100d;

        System.out.println(getSampleNameFromPath(inFile) + "\t" + varCallPos + "\t" + covround + "\t" + (100 - nperc) + "\t" + refCallPos + "\t" + allPos + "\t" + noCallPos + "\t" + discardedRefCall + "\t" + discardedVarCall + "\t" + filteredVarCall + "\t" + unknownCall);


        //////////////////////
        //END -- Parse VCFs //
        //////////////////////


        //Write draftN

        StringBuffer tmpSeq;

        tmpSeq = new StringBuffer();

        for (int pos = 1; pos <= calls.length; pos++) {
            tmpSeq.append(calls[pos - 1]);
        }

        FASTAWriter fw = new FASTAWriter(outFileDraft, draftName + "_draftN", tmpSeq.toString());



        //Write draftUncertain
        tmpSeq = new StringBuffer();

        for (int pos = 1; pos <= calls.length; pos++) {
            tmpSeq.append(uncertainCalls[pos - 1]);
        }

        fw = new FASTAWriter(outFileNR1234, draftName + "_draftUncertain", tmpSeq.toString());


        //Write draftRefMod
        tmpSeq = new StringBuffer();

        String tmpchar;


        for (int pos = 1; pos <= calls.length; pos++){
            tmpchar = calls[pos - 1];
            if (tmpchar.equals(String.valueOf(nChar))) {
                tmpSeq.append(refGenome.charAt(pos - 1));
            } else {
                tmpSeq.append(tmpchar);
            }
        }

        fw = new FASTAWriter(outFileRefMod, draftName + "_draftRefMod", tmpSeq.toString());

    }

    private void printUsage(PrintStream stream) {
        stream.println("VCF2Genome generates a draft genome sequence from a GATK vcf file. The file has to contain a call for each site (also non-variants).");
        stream.println("VCF2Genome is only applicable to single reference sequences. Multiple chromosomes are not supported. UnifiedGenotyper or HaplotypeCaller output is supported though.");
    }

    public void doMain(String[] args){
        CmdLineParser parser = new CmdLineParser(this);
        parser.setUsageWidth(120);

        try {
            parser.parseArgument(args);
            this.hadArguments = true;

            if(parser.getArguments().isEmpty()){
                throw new CmdLineException(parser, "No argument is provided unfortunately.");
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
            System.err.println();
            System.err.println("    Example: java -jar VCF2Genome.jar"+parser.printExample(OptionHandlerFilter.REQUIRED));
            return;
        }

        if(help){
            printUsage(System.err);
            parser.printUsage(System.err);
            return;
        }
    }


    public char getAmbiguousBase(char c1, char c2) {
        Set<Character> chars = new HashSet<Character>();
        chars.add(c1);
        chars.add(c2);

        if (c1 == c2)
            return c1;

        if (chars.contains('G') && chars.contains('A'))
            return 'R';
        if (chars.contains('T') && chars.contains('C'))
            return 'Y';
        if (chars.contains('G') && chars.contains('T'))
            return 'K';
        if (chars.contains('A') && chars.contains('C'))
            return 'M';
        if (chars.contains('G') && chars.contains('C'))
            return 'S';
        if (chars.contains('A') && chars.contains('T'))
            return 'W';

        System.err.println("Illegal arguments in function getAmbiguousBase: " + c1 + ", " + c2);

        return 'N';
    }

    /**
     * This gets the sample name without the extension from the path.
     * @param sampleNameWithPath
     * @return
     */

    public static String getSampleNameFromPath(String sampleNameWithPath) {
        return Files.getNameWithoutExtension(sampleNameWithPath);
    }


}
