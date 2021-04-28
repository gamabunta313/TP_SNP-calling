# TP53_SNP-calling
This pipeline is designed to analyze TP53 custom amplicon NGS data,
performing variant calling and validating calls with a R-based workflow.
For SNP validation we adopted 2 validiation pipelines: Fisher's exact test on Variant allele frequency (VAF) and Error modelling procedure based on distribution transformation.  
Both procedures need a custom count-database generated from wild-type (WT) samples. We test the efficiency of this workflow with a database of 362 WT patients.
