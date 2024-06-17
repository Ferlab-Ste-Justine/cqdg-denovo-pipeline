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
CONGE-XXX       CONGE-XXX-01.hard-filtered.gvcf.gz   CONGE-XXX-02.hard-filtered.gvcf.gz   CONGE-XXX-03.hard-filtered.gvcf.gz
CONGE-YYY       CONGE-YYY-01.hard-filtered.gvcf.gz   CONGE-YYY-02.hard-filtered.gvcf.gz   CONGE-YYY-03.hard-filtered.gvcf.gz
```

2.  The second format is used in older data and includes the sequencing type (WGS or WES)

**sampleV2.tsv**

_FAMILY_ID_ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; _SEQUENCING_TYPE_ &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;_Patient1_File_&nbsp; &nbsp; &nbsp;&nbsp; &nbsp;_Patient2_File_&nbsp; &nbsp; &nbsp; &nbsp;&nbsp;_Patient3_File_
```tsv
CONGE-XXX       WES       CONGE-XXX-01.hard-filtered.gvcf.gz   CONGE-XXX-02.hard-filtered.gvcf.gz   CONGE-XXX-03.hard-filtered.gvcf.gz
CONGE-YYY       WES       CONGE-YYY-01.hard-filtered.gvcf.gz   CONGE-YYY-02.hard-filtered.gvcf.gz   CONGE-YYY-03.hard-filtered.gvcf.gz
```


The file format can be chosen with the "**sampleFileFormat**" parameter (either "V1" or "V2", default "V1"). Note that both types are tab-delimited (.tsv)

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

Generally, we use the Homo_sapiens_assembly38.fasta as referenceGenome (see Resources)



Next, we also need broader references, which are contained in a path defined by the **broad** parameter.

The broad directory must contain the following files:

- The interval list which determines the genomic interval(s) over which we operate: filename of this list must be defined with the **intervalsFile** parameter
- Highly validated variance ressources currently required by VQSR. ***These are currently hard coded in the pipeline!***
  - HapMap file : hapmap_3.3.hg38.vcf.gz
  - 1000G omni2.5 file : 1000G_omni2.5.hg38.vcf.gz
  - 1000G reference file : 1000G_phase1.snps.high_confidence.hg38.vcf.gz
  - SNP database : Homo_sapiens_assembly38.dbsnp138.vcf.gz

 
Finally, the vep cache directory must be specified with **vepCache**, which is usually created by vep itself on first installation.
Generally, we only need the human files obtainable from https://ftp.ensembl.org/pub/release-112/variation/vep/homo_sapiens_vep_112_GRCh38.tar.gz

### Stub run
The -stub-run option can be added to run the "stub" block of processes instead of the "script" block. This can be helpful for testing.

🚧

Parameters summary
-----

| Parameter name | Required? | Accepted input |
| --- | --- | --- |
| `sampleFile` | _Required_ | file |
| `sampleFileFormat` | _Optional_ | `V1` or `V2`, default `V1` |
| `sequencingType` | _Optional_ | `WGS` or `WES`, default `WGS` |
| `referenceGenome` | _Required_ | path |
| `referenceGenomeFasta` | _Required_ | file |
| `broad` | _Required_ | path |
| `intervalsFile` | _Required_ | list of genome intervals |
| `vepCache` | _Required_ | path |

Pipeline Output
-----
🚧

Resources
-----
The documentation of the various tools used in this workflow are available here:

[Nextflow](https://www.nextflow.io/docs/latest/index.html)

[bcftools](https://samtools.github.io/bcftools/bcftools.html)

**GATK**:
- [CombineGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037593911-CombineGVCFs)
- [GenotypeGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037057852-GenotypeGVCFs)
- [VariantRecalibrator](https://gatk.broadinstitute.org/hc/en-us/articles/360035531612-Variant-Quality-Score-Recalibration-VQSR)
- [VariantFiltration](https://gatk.broadinstitute.org/hc/enus/articles/360041850471-VariantFiltration))

[VEP](https://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html)

**Reference files**
🚧
