# VCF2Genome
A tool to create a draft genome file out of a GATK VCF file

[![Build Status](https://travis-ci.org/apeltzer/VCF2Genome.svg?branch=master)](https://travis-ci.org/apeltzer/VCF2Genome)
[![codecov](https://codecov.io/gh/apeltzer/VCF2Genome/branch/master/graph/badge.svg)](https://codecov.io/gh/apeltzer/VCF2Genome)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/vcf2genome/README.html)

Download via [GitHub Releases](https://github.com/apeltzer/VCF2Genome/releases) or via [Bioconda](https://bioconda.github.io/recipes/vcf2genome/README.html).

Author: Alexander Herbig <herbig@shh.mpg.de> (v0.84), Alexander Peltzer (v0.90+).

Contact Alexander Peltzer<peltzer@shh.mpg.de> for questions regarding the tool or via GitHub and/or open a ticket here.


# Basic Usage description

## Help message
You can see a help when running the tool with `-h`. This generates the following help message:

```
Option "-draft" is required
 -draft VAL                 : draft contains Ns where no call can be made. RefMod contains reference calls instead at
                              these positions.
 -draftname DRAFT_SEQ_NAME  : Name of the draft sequence.
 -h                         : Display this help information and exit. (default: true)
 -in VAL                    : input VCF file
 -minc MIN_COVERAGE_FOR_SNP : Minimum coverage / reads confirming the call.
 -minfreq MIN_SNP_FREQUENCY : Minimum fraction of reads supporting the called nucleotide.
 -minq MIN_QUAL_SCORE       : Minimum quality score. For UG: Phred scaled quality score. For HC genome quality score.
 -ref VAL                   : reference genome in FastA format
 -refMod VAL                : More precise uncertainty encoding. N: Not covered or ambiguous. R: Low coverage but looks
                              like Ref. a,c,t,g (lower case): Low coverage but looks like SNP.
 -uncertain VAL             : Special 1234 encoded FastA output.

Example: java -jar VCF2Genome.jar -draft VAL -draftname DRAFT_SEQ_NAME -in VAL -minc MIN_COVERAGE_FOR_SNP -minfreq MIN_SNP_FREQUENCY -minq MIN_QUAL_SCORE -ref VAL -refMod VAL -uncertain VAL
```

## Option documentation

### Example call

```
java -jar VCF2Genome.jar -draft my_output_genome.fasta -draftname "My_Fancy_Genome_Name" -in my_input.vcf -minc 5 -minfreq 0.8 -minq 30 -ref myreference_genome.fasta -refMod output.refMod -uncertain 1234_output.fasta
```

### `-draft`

Name of the output file to which the FastA genome sequence should be written. Contains Ns where no call can be made. 

### `-draftname`

Name of the draft sequence inside the FastA file (header of the FastA entry that is created).

### `-in`

Name of the input VCF file in VCF4.0/4.1 format. 

### `-minc`

Minimum coverage / reads confirming the call required. 

### `-minq`

Minimum quality threshold used for filtering the calls.

### `-minfreq`
Minimum fraction of reads supporting the called nucleotide.

### `-ref`
Reference genome used in FastA format.

### `-refMod`

Path to refMod format output file. This contains a more detailed output encoding than just including `N` at unclear positions. Useful for further investigation of some sites for example. 
N: Not covered or ambiguous. R: Low coverage but looks like Reference call. a,c,t,g (lower case): Low coverage but looks like SNP.

### `-uncertain`

Path to uncertainty encoded output file in a special 1234 format for some downstream tools. 

## Additional comments

### VCF Compatibility

Note that this tool was written a couple of years ago for reconstructing genomes from GATK UnifiedGenotyper VCF output files. It may work with other genotypers providing the same kind of VCF4.0/VCF4.1 format, but might not work well with data originating for example from GATK HaplotypeCaller. The tool requires an [`EMIT_ALL_SITES`](https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php) compatible VCF input file.

### SNP/Indel handling

This tool is currently unable to handle indels properly due to the index handling procedure in the software itself. SNPs are fine.

### 
