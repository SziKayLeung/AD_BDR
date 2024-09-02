- **1.** E.Tseng Iso-Seq [IsoPhase](https://github.com/Magdoll/cDNA_Cupcake/wiki/IsoPhase:-Haplotyping-using-Iso-Seq-data)
    - Caveats
        - only developed for substitution SNPs (not indels) and for diploid organisms
        - work off from Iso-Seq collapse (require read stat files etc)
    1. Create mpileups for each gene (samtools)
    2. call variants from mpileups
        1. minimum 10 full-length reads coverage
        2. allow substitution error 0.005
        3. SNPs cannot be near homo-polymers (defined as sequences with 4 or more identical nucleotides)
        4. for variants that pass threshold, perform Fisher exact test with Bonferroni correction 
    3. for each read, assigned haplotype
    4. phase isoforms collapsed from reads
    5. clean haplotypes by using all possible combinations of haplotypes (assuming diploid) and the haplotype with the smallest difference between that and the combinations are assumed the best
    6. create a vcf file summarising all the outputs
- **2.** [Lorals](https://github.com/LappalainenLab/lorals)
  - Custom scripts to identify reads with reference and alternative allele

**Other tools**

- isoLaser: https://github.com/gxiaolab/isoLASER/issues/1 though issues and not ready for release
- https://github.com/vladimirsouza/lrRNAseqVariantCalling: publication to trial DeepVariant and Clair3 (tools that were developed for variant calling on genomic DNA) on cDN

[Transformation of alignment files improves performance of variant callers for long-read RNA sequencing data - Genome Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-02923-y)
