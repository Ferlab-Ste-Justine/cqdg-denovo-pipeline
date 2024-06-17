nextflow.enable.dsl = 2

include { VQSR } from "./subworkflows/vqsr"
include { hardFiltering } from './modules/hardFilter'
include { splitMultiAllelics        } from './modules/vep'
include { vep                       } from './modules/vep'
include { tabix                     } from './modules/vep'
include { validateParameters; paramsHelp; paramsSummaryLog} from 'plugin/nf-schema'


enum SequencingType {
    WGS,
    WES
}

enum SampleFileFormat {
    V1,
    V2
}
//Nf-core schema functions
if (params.help) {
    log.info paramsHelp("nextflow run Ferlab-Ste-Justine/cqdg-denovo-pipeline -r v1.0.2 -params-file params.json")
    exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

/**
Parse the sample file (tsv) and returns the 2 following channels as a map:
    - meta : (familyId, metadata_dictionary)
    - files : (familyId, [single.patient.file.gvcf.gz, single.patient.file.gvcf.gz.tbi])

Note: one line per patient is emitted for the meta channel, i.e. it contains duplicates.

The metadata dictionary contains the following keys:
- size: number of patients in the family
- sequencingType: type of sequencing data, can be either 'WGS' or 'WES'

## Input Format ##

You can specify the input format via configuration parameter `sampleFileFormat`. 
We support 2 different formats: `V1` and `V2`. The default is `V1`.

Columns in V1 format:
```
familyId    filePatient1    filePatient2    ...
```

Columns in V2 format:
```
familyId    SequencingType  filePatient1 filePatient2   ...
```

In the V1 format, the sequencing type is assumed to be the same for all samples and is specified through
the configuration parameter `sequencingType` (default to `WGS`).

## Potiential improvements ##

Consider the following improvements:
- Rethink propagation of metadata in the workflow: https://training.nextflow.io/advanced/metadata/
- Add validation step
*/
def sampleChannel() {    
    def rowMapper = getRowMapper() 

    return Channel.fromPath(file("$params.sampleFile"))
        .splitCsv(sep: '\t', strip: true)
        .map(rowMapper)
        .flatMap { it ->
            return it.files.collect{f -> [familyId: it.familyId, sequencingType: it.sequencingType, size: it.files.size(), file: f]};             
        }.multiMap { it ->
            meta: tuple(it.familyId, [size: it.size, sequencingType: it.sequencingType])
            files: tuple(it.familyId, file("${it.file}*"))
        }
}

def getSampleFileFormat() {
    if (!params.sampleFileFormat) {
        log.warn("Using default value `V1` for parameter `sampleFileFormat`")
        params.sampleFileFormat="V1"
    }
    return findParamInEnum("sampleFileFormat", params.sampleFileFormat.toUpperCase(), SampleFileFormat)
}

def getSequencingType() {
    if (!params.sequencingType) {
        log.warn("Using default value `WGS` for parameter `sequencingType`")
        params.sequencingType="WGS"
    }
    return findParamInEnum("sequencingType", params.sequencingType.toUpperCase(), SequencingType)
}

def findParamInEnum(paramName, paramValue, enumInstance) {
    def validValues = enumInstance.values()*.name()
    if (!validValues.contains(paramValue)) {
        def validValuesStr = validValues.collect{"`$it`"}.join(", ")
        error("Invalid value for parameter `$paramName`: `$paramValue`. Possible values are: $validValuesStr")
    }
    return enumInstance.valueOf(paramValue)
}

/**
Get row mapper that match the configured sample file format

Note: it is returned as a closure to guaranty the compatibility with nextflow channel operators
*/
def getRowMapper() {
    def format = getSampleFileFormat()

    if (format == SampleFileFormat.V1) {
        def sequencingType = getSequencingType()
        return {columns -> rowMapperV1(columns, sequencingType)}
    }
    return {columns -> rowMapperV2(columns)}
}

/** 
Transform a row from the sample file in V1 format from a list structure to a map structure.
*/
def rowMapperV1(columns, sequencingType) {
    return [
        familyId: columns[0],
        sequencingType: sequencingType,
        files: columns.tail()
    ]
}

/** 
Transform a row from the sample file in V2 format from a list structure to a map structure
*/
def rowMapperV2(columns) {
    return [
        familyId: columns[0],
        sequencingType: columns[1].toUpperCase() as SequencingType,
        files: columns[2..-1]
    ]
}

/**
Tag variants that are probable artifacts

In the case of whole genome sequencing data, we use the vqsr procedure.
For whole exome sequencing data, since the vqsr procedure is not supported, we use
a hard filtering approach.
*/
def tagArtifacts(inputChannel, metadataChannel, hardFilters) {
    def inputSequencingTypes = inputChannel.join(metadataChannel)
    
    def wgs = inputSequencingTypes.filter{it[2].sequencingType == SequencingType.WGS}.map(it -> it.dropRight(1))
    def wes = inputSequencingTypes.filter{it[2].sequencingType == SequencingType.WES}.map(it -> it.dropRight(1))

    def wgs_filtered = VQSR(wgs)
    def wes_filtered = hardFiltering(wes, hardFilters)

    return wgs_filtered.concat(wes_filtered)
}

/**
Exclude MNPs detected with bedtools
*/
process excludeMNPs {
    label 'medium'
    
    input:
    tuple val(familyId), path(gvcfFile)

    output:
    tuple val(familyId), path("*filtered.vcf.gz*")

    script:
    def exactGvcfFile = gvcfFile.find { it.name.endsWith("vcf.gz") }
    def uuid = UUID.randomUUID().toString()
    // --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
    // bcftools view --exclude-type mnps  ${exactGvcfFile} -O z -o ${familyId}.${uuid}.filtered.gvcf.gz
    """
    set -e
    echo $familyId > file
    bcftools filter -e 'strlen(REF)>1 & strlen(REF)==strlen(ALT) & TYPE="snp"' ${exactGvcfFile} | bcftools norm -d any -O z -o ${familyId}.${uuid}.filtered.vcf.gz
    bcftools index -t ${familyId}.${uuid}.filtered.vcf.gz
    """
    stub:

    def exactGvcfFile = gvcfFile.find { it.name.endsWith("vcf.gz") }
    def uuid = UUID.randomUUID().toString()
    """
    touch ${familyId}.${uuid}.filtered.vcf.gz
    touch ${familyId}.${uuid}.filtered.vcf.gz.tbi
    """

}

/**
Combine per-sample gVCF files into a multi-sample gVCF file 
*/
process importGVCF {

    label 'medium'

    input:
    tuple val(familyId), path(gvcfFiles)
    path referenceGenome
    path broadResource

    output:
    tuple val(familyId), path("*combined.gvcf.gz*")

    script:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("vcf.gz") }.collect { "-V $it" }.join(' ')

    """
    echo $familyId > file
    gatk -version
    gatk --java-options "-Xmx8g"  CombineGVCFs -R $referenceGenome/${params.referenceGenomeFasta} $exactGvcfFiles -O ${familyId}.combined.gvcf.gz -L $broadResource/${params.intervalsFile}
    """       

    stub:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("vcf.gz") }.collect { "-V $it" }.join(' ')

    """
    touch ${familyId}.combined.gvcf.gz
    """       

}    

/**
Keep only SNP and Indel 
*/
process genotypeGVCF {
    label 'geno'

    input:
    tuple val(familyId), path(gvcfFile)
    path referenceGenome

    output:
    tuple val(familyId), path("*genotyped.vcf.gz*")

    script:
    def exactGvcfFile = gvcfFile.find { it.name.endsWith("vcf.gz") }
    """
    echo $familyId > file
    gatk -version
    gatk --java-options "-Xmx8g" GenotypeGVCFs -R $referenceGenome/${params.referenceGenomeFasta} -V $exactGvcfFile -O ${familyId}.genotyped.vcf.gz
    """

    stub:
    def exactGvcfFile = gvcfFile.find { it.name.endsWith("vcf.gz") }
    """
    touch ${familyId}.genotyped.vcf.gz
    """
}

workflow {
    referenceGenome = file(params.referenceGenome)
    broad = file(params.broad)
    vepCache = file(params.vepCache)
    file(params.outputDir).mkdirs()

    sampleChannel().set{ sampleFile }

    //Exclude MNPs and recombine files per family id
    //(familyId, [filtered.patientUUID1.gvcf.gz, filtered.patientUUID1.gcvf.gz.tbi, filtered.patientUUID2.gcvf.gz, filtered.patientUUID2.gvcf.gz.tbi... ])
    filtered = excludeMNPs(sampleFile.files)
                    .join(sampleFile.meta)
                    .map{familyId, files, meta -> tuple( groupKey(familyId, meta.size), files)}
                    .groupTuple()
                    .map{ familyId, files -> tuple(familyId, files.flatten())}
    //Using 2 as threshold because we have 2 files per patient (gcvf.gz, gvcf.gz.tbi)
    filtered_one = filtered.filter{it[1].size() == 2}
    filtered_mult = filtered.filter{it[1].size() > 2}

    //Combine per-sample gVCF files into a multi-sample gVCF file
    DB = importGVCF(filtered_mult, referenceGenome,broad)
                    .concat(filtered_one)


    //Perform joint genotyping on one or more samples
    vcf = genotypeGVCF(DB, referenceGenome)

    //tag variants that are probable artifacts
    vcfWithTags = tagArtifacts(vcf, sampleFile.meta, params.hardFilters)

    //tag frequent mutations in the population
    s = splitMultiAllelics(vcfWithTags, referenceGenome) 

    //Annotating mutations
    vep(s, referenceGenome, vepCache)
    tabix(vep.out) 
}