CQDG denovo pipelines
======
Introduction
------
🚧 WIP : this pipeline recombine gvcf for family's samples in order to facilitate denovo identification.

(Include Damien's workflow schema?)

Usage
-----
### Samples
The workflow will accept sample data in two format (called V1 and V2). The path to the sample file must be specified with the "**sampleFile**" parameter.

1.  The first format is used by default and looks as follows:

**sampleV1.tsv**

_FAMILY_ID_ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; _Patient1_File_&nbsp; &nbsp; &nbsp;&nbsp; &nbsp;_Patient2_File_&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;_Patient3_File_
```tsv
CONGE-320       CONGE-320-01.hard-filtered.gvcf.gz   CONGE-320-02.hard-filtered.gvcf.gz   CONGE-320-03.hard-filtered.gvcf.gz
CONGE-345       CONGE-345-01.hard-filtered.gvcf.gz   CONGE-345-02.hard-filtered.gvcf.gz   CONGE-345-03.hard-filtered.gvcf.gz
```

2.  The second format is used in older data and includes the sequencing type (WGS or WES)

**sampleV2.tsv**

_FAMILY_ID_ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; _SEQUENCING_TYPE_ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;_Patient1_File_&nbsp; &nbsp; &nbsp;&nbsp; &nbsp;_Patient2_File_&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;_Patient3_File_
```tsv
CONGE-320       WES       CONGE-320-01.hard-filtered.gvcf.gz   CONGE-320-02.hard-filtered.gvcf.gz   CONGE-320-03.hard-filtered.gvcf.gz
CONGE-345       WES       CONGE-345-01.hard-filtered.gvcf.gz   CONGE-345-02.hard-filtered.gvcf.gz   CONGE-345-03.hard-filtered.gvcf.gz
```


The file format can be chosen with the "**sampleFileFormat**" parameter (either "V1" or "V2", default "V1").

Next, if the file format is "V1", the sequencing type can be specified with the "**sequencingType**" parameter (either "WGS" for Whole Genome Sequencing or "WES" for Whole Exome Sequencing, default "WGS")

> [!NOTE]
> The sequencing type also determines the type of variant filtering the pipeline will use.
> 
> In the case of Whole Genome Sequencing, VQSR (Variant Quality Score Recalibration) is used (preferred method).
> 
> In the case of Whole Exome Sequencing, Hard-filtering needs to be used.

### References
Reference files are necessary at multiple steps of the workflow, notably for joint-genotyping,the variant effect predictor (VEP) and VQSR. 

Specifically, we need a reference genome directory and filename specified with the **referenceGenome** and **referenceGenomeFasta** parameters respectively. 

⚠️ _(TO DO: I think these two parameters could be combined, I don't think we use any other file than the Fasta in the referenceGenome directory.)_ ⚠️

Generally, we use the Homo_sapiens_assembly38.fasta as referenceGenome (see Ressources)



Next, we also need broader references, which are contained in a path defined by the **broad** parameter.

The broad directory must contain the following files:

- The interval list which determines the genomic interval(s) over which we operate: filename of this list must be defined with the **intervalsFile** parameter
- Highly validated variance ressources currently required by VQSR. ***These are currently hard coded in the pipeline!***
  - HapMap file : hapmap_3.3.hg38.vcf.
  - 1000G omni2.5 file : 1000G_omni2.5.hg38.vcf
  - 1000G reference file : 1000G_phase1.snps.high_confidence.hg38.vcf
  - SNP database : Homo_sapiens_assembly38.dbsnp138.vcf
🚧

Parameters summary
-----
**sampleFile**: &nbsp; _Required_ &nbsp; path &nbsp; 

**sampleFileFormat**: &nbsp; _Optional_ &nbsp; "V1" or "V2" &nbsp; Default: "V1"

**sequencingType**: &nbsp; _Optional_ &nbsp; "WGS" or "WES" &nbsp; Default: "WGS"
🚧

Pipeline Output
-----
🚧

