GATK Variant Caller, Advanced Readme
====================================

Introduction
------------

This is an app which runs the UnifiedGenotyper module within the GenomeAnalysisToolkit
(GATK) to produce variant calls from a set of mapped reads.

UnifiedGenotyper
----------------

* <a href="http://www.broadinstitute.org/gatk/">GATK site</a>
* <a href="http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html">UnifiedGenotyper Documentation (includes command options)</a>

UnifiedGenotyper takes a mapped reads object and a
reference genome and produces a set of variant calls. UnifiedGenotyper will
call both SNPs and small indels. Its ability to call indels decreases as the
size of the indel increases.

The command line options listed below are specific to the UnifiedGenotyper
component of GATK. These commands can be combined with certain general GATK
options. The general options will be listed below the specific ones:

### Options

Most input parameters correspond to options that GATK would expect running in standalone linux, e.g.:

    -Xmx4g -jar GenomeAnalysisTK -T UnifiedGenotyper -I input.bam -R hg19.fa \
      -o output.vcs -out_mode EMIT_ALL_CONFIDENT_SITES

<table>
<tr><th>Name</th><th>Label</th><th>Corresponding command line option/Notes</th><th>Class/Type</th><th>Default</th></tr>
<tr><td>reference</td><td>Reference Genome</td><td>-R &lt;ref.fa&gt;</td><td>record of type:ContigSet</td><td>Mandatory</td></tr>
<tr><td>mappings</td><td>Mappings Object</td><td>- I &lt;input.bam&gt; </td><td>gtable of type LetterMappings.</td><td>Mandatory</td></tr>
<tr><td>output_mode</td><td>Output Mode</td><td>-out_mode</td><td>[EMIT_VARIANTS_ONLY, EMIT_ALL_CONFIDENT_SITES, EMIT_ALL_SITES]</td><td>EMIT_ALL_CONFIDENT_SITES (VARIANTS_ONLY is fastest. EMIT_ALL_SITES is much slower and generates very large output files. An alternative to EMIT_ALL_SITES is to EMIT_ALL_CONFIDENT_SITES and to change the "emit threshold" parameter to a low value). </td></tr>
<tr><td>call_confidence</td><td>Call Threshold</td><td>-stand_call_conf (What Phred confidence is required to consider a site confidently called)</td><td>double</td><td>30.0</td></tr>
<tr><td>emit_confidence</td><td>Emit Threshold</td><td>-stand_emit_conf (This value determines what GATK will write as a site, even if the confidence is below the call confidence). In this case it is marked as below filter quality.</td><td>double</td><td>30.0</td></tr>
<tr><td>pcr_error_rate</td><td>PCR Error Rate</td><td>-pcr_error (Used to calculate likelihoods)</td><td>double</td><td>1.0E-4</td></tr>
<tr><td>heterozygosity</td><td>SNP Heterozygosity</td><td>-hets (Used to calculate prior likelihoods)</td><td>double</td><td>0.001</td></tr>
<tr><td>indel_heterozygosity</td><td>Indel Heterozygosity</td><td>-indelHeterozygosity</td><td>double</td><td>1.25E-4</td></tr>
<tr><td>genotype_likelihood_model</td><td>Call SNPs or Indels</td><td>-glm</td><td>[SNP, INDELS, BOTH]</td><td>SNP (Joe used BOTH for Labcorp work. AC)</td></tr>
<tr><td>minimum_base_quality</td><td>Min Base Quality</td><td>-mbq (Minimum base mapping quality required to attempt variant calling)</td><td>int</td><td>17</td></tr>
<tr><td>max_alternate_alleles</td><td>Max Alleles</td><td>-maxAlleles (When possible alternates are discovered, how many are considered for genotyping)</td><td>int</td><td>3 (Computational cost scales exponentially with increasing alternate alleles)</td></tr>
<tr><td>max_deletion_fraction</td><td>Max Deletion Fraction</td><td>-deletions (Maximum number of reads with deletions spanning a locus for it to be callable)</td><td>double</td><td>0.05</td></tr>
<tr><td>min_indel_count</td><td>Min Indel</td><td>-minIndelCnt (The minimum number of consensus indels required to attempt genotyping a site as an indel)</td><td>int</td><td>5</td></tr>
<tr><td>non_reference_probability_model</td><td>Probability Model</td><td>-pnrm Which model is used to calculate the probability that a site is a variant, and the genotypes</td><td>[EXACT, GRID_SEARCH]</td><td>EXACT</td></tr>
<tr><td>calculate_BAQ</td><td>Calculate BAQ</td><td>-baq (Per-base quality scores) Another quality metric for mappings that can be used for filters</td><td>[OFF, CALCULATE_AS_NECESSARY, RECALCULATE]</td><td>OFF</td></tr>
<tr><td>BAQ_gap_open_penalty</td><td>BAQ Gap Open Penalty</td><td>-baqGOP (Phred scaled)</td><td>double</td><td>40.0 (30.0 is suggested for WGS)</td></tr>
<tr><td>no_output_SLOD</td><td>Don't Output Slod</td><td>-nosl (Recalibrates score using a model constructed from the high quality sites, outputs a "more accurate" probability for variant call</td><td>Boolean</td><td>False</td></tr>
<tr><td>intervals_to_process</td><td>Process Intervals</td><td>-L</td><td>List of regions to process in the form of any number of -L CHR:lo-hi (e.g. -L chr1:100-200 -L chr1:300-500 -L chr5:1000-2000).</td><td>None</td>
<tr><td>intervals_to_exclude</td><td>Exclude Intervals</td><td>-XL</td><td>As above, except these intervals are excluded</td><td>None</td></tr>
<tr><td>intervals_merging</td><td>Interval Merge</td><td>-im If multiple intervals are given, this controls the rules by which they specify regions. Possible options are: "UNION" and "INTERSECTION". Selecting union will add all regions while intersection will only take regions contained in all of the specified intervals</td><td>[UNION, INTERSECTION]<td>ALL</td></tr>
<tr><td>downsample_to_coverage</td><td>Downsample Coverage</td><td>-dcov (Randomly downsamples reads at each locus to a given number)</td><td>int</td><td>Not present</td></tr>
<tr><td>downsample_to_fraction</td><td>Downsample Fraction</td><td>-dfrac (Randomly downsamples reads to a fraction of the total)</td><td>double</td><td>1.0</td></tr>
<tr><td>nondeterministic</td><td>Nondeterministic</td><td>-ndrs (Makes GATK run non-deterministically)</td><td>bool</td><td>False</td></tr>
</table>
