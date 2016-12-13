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





    @Test
    public void testVCF2GenomeForConsistency() throws Exception {
        //<in.vcf> <reference-genome.fasta> <draft.fasta> <refMod.fasta> <uncertain.fasta> <minimum quality score> <minimum coverage> <minimum SNP allel frequency> [draft sequence name]");
        //Main.main();

        VCF2Genome vcf2Genome = new VCF2Genome(new String[]{"build/classes/test/VCF2Genome_Test_Subset.vcf", "build/classes/test/NC_021490.2.fasta", "build/classes/test/draft.fasta", "build/classes/test/refMod.fasta", "build/classes/test/uncertain.fasta", "30", "5", "0.9", "test"});


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
        VCF2Genome vcf2Genome = new VCF2Genome(new String[]{"build/classes/test/VCF2Genome_Test_SubsetFailSortOrder.vcf", "build/classes/test/NC_021490.2.fasta", "build/classes/test/draft_fail.fasta", "build/classes/test/refMod_fail.fasta", "build/classes/test/uncertain_fail.fasta", "30", "5", "0.9", "test"});
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
        assertEquals("VCF2Genome - 0.9 (2016-12-08)\n" +
                "by Alexander Herbig and Alexander Peltzer\n" +
                "\n" +
                "--- How to use VCF2Genome ---\n" +
                "\n" +
                "USAGE: java -jar VCF2Genome.jar <in.vcf> <reference-genome.fasta> <draft.fasta> <refMod.fasta> <uncertain.fasta> <minimum quality score> <minimum coverage> <minimum SNP allel frequency> [draft sequence name]\n" +
                "EXAMPLE: java -jar VCF2Genome.jar in.vcf reference.fasta draft.fasta refMod.fasta uncertain.fasta 40 7 0.9 my_new_draft_genome\n" +
                "\n" +
                "VCF2Genome generates a draft genome sequence from a GATK vcf file. The file has to contain a call for each site (also non-variants).\n" +
                "VCF2Genome is only applicable to single reference sequences. Multiple chromosomes are not supported.\n" +
                "VCF2Genome is only recommended for bacterial genomes.\n" +
                "Example GATK command line:\n" +
                "GenomeAnalysisTK.jar -T UnifiedGenotyper -R <referenceGenome.fasta> --output_mode EMIT_ALL_SITES -o inputForVCF2Genome.vcf -I <input.bam>\n" +
                "GATK needs bam files with read groups. Read groups can be added using PicardTools:\n" +
                "AddOrReplaceReadGroups.jar I=In.bam O=Out.bam LB=NA PL=illumina PU=NA SM=NA VALIDATION_STRINGENCY=SILENT\n" +
                "\n" +
                "More details:\n" +
                "draft vs. refMod: 'draft' contains 'N's where no call can be made, 'refMod' contains reference bases in this case.\n" +
                "uncertain.fasta: More precise uncertainty encoding. N: Not covered or ambiguous. R: Low coverage but looks like Ref. a,c,t,g (lower case): Low coverage but looks like SNP.\n" +
                "minimum quality score: The score is given in the 6th column of the vcf file. Phred-scaled quality score for the call. High QUAL scores indicate high confidence calls. -10log_10 p(no variant)\n" +
                "minimum coverage: minimum number of reads confirming the call.\n" +
                "minimum SNP allel frequency: minimum fraction of reads containing the called nucleotide.\n" +
                "\n" +
                "In case of any questions contact Alexander Herbig and Alexander Peltzer (herbig@shh.mpg.de, peltzer@shh.mpg.de).\n" , outContent.toString());
    }

    @Test
    public void testAmbigousBases(){
        char c1 = 'G';
        char c2 = 'A';
        char c3 = 'T';
        char c4 = 'C';

        assertEquals('R', VCF2Genome.getAmbiguousBase(c1,c2));
        assertEquals('Y', VCF2Genome.getAmbiguousBase(c3,c4));
        assertEquals('K', VCF2Genome.getAmbiguousBase(c1,c3));
        assertEquals('M', VCF2Genome.getAmbiguousBase(c2,c4));
        assertEquals('S', VCF2Genome.getAmbiguousBase(c1,c4));
        assertEquals('W', VCF2Genome.getAmbiguousBase(c2,c3));
        assertEquals('N', VCF2Genome.getAmbiguousBase('U', c1)); //One example!
    }
}
