CQDG denovo pipelines
======
Introduction
------
🚧 WIP : this pipeline recombine gvcf for family's samples in order to facilitate denovo identification.

(Include Damien's workflow schema?)

Usage
-----
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

Next, if the file format is "V1", the sequencing type can be specified with the "**sequencingType**" parameter (either "WGS" or "WES", default "WGS")




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

