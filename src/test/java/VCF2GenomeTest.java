import FastAIO.FASTAParser;
import VCFTools.AmbigoutyCalculator;
import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.*;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by peltzer on 08/12/2016.
 */
public class VCF2GenomeTest {
    private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
    private final ByteArrayOutputStream errContent = new ByteArrayOutputStream();


    public VCF2GenomeTest(){}


    /**
     * This tests only the UG behaviour. We might need new testcases for the HaplotypeCaller behaviour as well then.
     * @throws Exception
     */

    @Test
    public void testVCF2GenomeForConsistency() throws Exception {

        VCF2Genome vcf2Genome = new VCF2Genome(new String[]{"-in", "build/classes/test/VCF2Genome_Test_Subset.vcf", "-ref", "build/classes/test/NC_021490.2.fasta", "-draft", "build/classes/test/draft.fasta", "-refMod", "build/classes/test/refMod.fasta", "-uncertain", "build/classes/test/uncertain.fasta", "-minq", "30", "-minc", "5", "-minfreq", "0.9", "-draftname", "test"});


        //Now compare output
        File draft_golden = new File("build/classes/test/draft_golden.fasta");
        File refmod_golden = new File("build/classes/test/refMod_golden.fasta");
        File uncertain_golden = new File("build/classes/main/uncertain_golden.fasta");

        File draft = new File("build/classes/test/draft.fasta");
        File refmod = new File("build/classes/test/refMod.fasta");
        File uncertain = new File("build/classes/main/uncertain.fasta");

        boolean draft_equal = FileUtils.contentEquals(draft,draft_golden);
        boolean refmod_equal = FileUtils.contentEquals(refmod,refmod_golden);
        boolean uncertain_equal = FileUtils.contentEquals(uncertain,uncertain_golden);

        assertTrue(draft_equal);
        assertTrue(refmod_equal);
        assertTrue(uncertain_equal);
    }

    @Test(expected = IOException.class)
    public void testVCF2GenomeFailMessage() throws Exception{
        VCF2Genome vcf2Genome = new VCF2Genome(new String[]{"-in", "build/classes/test/VCF2Genome_Test_SubsetFailSortOrder.vcf", "-ref", "build/classes/test/NC_021490.2.fasta", "-draft", "build/classes/test/draft_fail.fasta", "-refMod", "build/classes/test/refMod_fail.fasta", "-uncertain", "build/classes/test/uncertain_fail.fasta", "-minq", "30", "-minc", "5", "-minfreq", "0.9", "-draftname", "test"});
    }


    @Before
    public void setUpStreams(){
        System.setOut(new PrintStream(outContent));
        System.setErr(new PrintStream(errContent));
    }

    @After
    public void cleanUpStreams(){
        System.setOut(null);
        System.setErr(null);
    }

    @Test
    public void testHelpAndStuff() throws Exception {
        VCF2Genome vcf2Genome = new VCF2Genome(new String[]{"-h"});
        assertEquals("VCF2Genome (v. 0.91 2017-01-13)\nby Alexander Herbig (<= v0.84) and Alexander Peltzer (>v0.84)\nherbig@shh.mpg.de, peltzer@shh.mpg.de\n\n", outContent.toString());
    }

    @Test
    public void testFastaReader() throws Exception {
        FASTAParser fastaParser = new FASTAParser("build/classes/test/mFasta.fasta");
        assertEquals(2, fastaParser.getFastaFile().size());

    }

    @Test
    public void testAmbigousBases() throws Exception {
        char c1 = 'G';
        char c2 = 'A';
        char c3 = 'T';
        char c4 = 'C';

        AmbigoutyCalculator ambigoutyCalculator = new AmbigoutyCalculator();


        assertEquals(c1, ambigoutyCalculator.getAmbiguousBase(c1,c1));
        assertEquals('R', ambigoutyCalculator.getAmbiguousBase(c1,c2));
        assertEquals('Y', ambigoutyCalculator.getAmbiguousBase(c3,c4));
        assertEquals('K', ambigoutyCalculator.getAmbiguousBase(c1,c3));
        assertEquals('M', ambigoutyCalculator.getAmbiguousBase(c2,c4));
        assertEquals('S', ambigoutyCalculator.getAmbiguousBase(c1,c4));
        assertEquals('W', ambigoutyCalculator.getAmbiguousBase(c2,c3));
        assertEquals('N', ambigoutyCalculator.getAmbiguousBase('U', c1)); //One example!
    }
}
