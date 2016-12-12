import org.apache.commons.io.FileUtils;
import org.junit.Test;

import java.io.*;
import java.net.URL;

import static org.junit.Assert.assertTrue;

/**
 * Created by peltzer on 08/12/2016.
 */
public class VCF2GenomeTest {
    private InputStream is;
    private InputStreamReader isr;
    private BufferedReader bfr;

    public VCF2GenomeTest(){}


    private void setUp(String path) {
        URL url = getClass().getResource(path);
        try {
            is = url.openStream();
        } catch (IOException e) {
            e.printStackTrace();
        }
        isr = new InputStreamReader(is);
        bfr = new BufferedReader(isr);
    }


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
}
