

nextflow.enable.dsl = 2

/**
Exclude MNPs detected with bedtools
*/

process excludeMNPs{
    label 'medium'

    container 'staphb/bcftools'
    
    input:
    tuple val(familyId), path(gvcfFile)

    output:
    tuple val(familyId), path("*filtered.vcf.gz*")

    script:
    def exactGvcfFile = gvcfFile.find { it.name.endsWith("gvcf.gz") }
    def uuid = UUID.randomUUID().toString()
    // --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
    // bcftools view --exclude-type mnps  ${exactGvcfFile} -O z -o ${familyId}.${uuid}.filtered.gvcf.gz
    """
    set -e
    echo $familyId > file
    bcftools filter -e 'strlen(REF)>1 & strlen(REF)==strlen(ALT) & TYPE="snp"' ${exactGvcfFile} | bcftools norm -d any -O z -o ${familyId}.${uuid}.filtered.vcf.gz
    bcftools index -t ${familyId}.${uuid}.filtered.vcf.gz
    """
}


process importGVCF {

    label 'medium'

    container 'broadinstitute/gatk'
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
    gatk --java-options "-Xmx8g"  CombineGVCFs -R $referenceGenome/${params.referenceGenomeFasta} $exactGvcfFiles -O ${familyId}.combined.gvcf.gz -L ${broadResource}/${params.intervalsFile}
    """       

}    

/**
Keep only SNP and Indel 
*/

process genotypeGVCF {
    label 'geno'

    container 'broadinstitute/gatk'

    input:
    tuple val(familyId), path(gvcfFile)
    path referenceGenome

    output:
    tuple val(familyId), path("*genotyped.vcf.gz*")

    script:
    def exactGvcfFile = gvcfFile.find { it.name.endsWith("gvcf.gz") }
    """
    echo $familyId > file
    gatk --java-options "-Xmx24g" GenotypeGVCFs -R $referenceGenome/${params.referenceGenomeFasta} -V $exactGvcfFile -O ${familyId}.genotyped.vcf.gz
    """
}




/**
Parse a tsv and return multiple chanels with 
    - sizes : size by family id
    - files : files by family id
*/
def sampleChannel() {
   Channel.fromPath(file("$params.sampleFile"))
               .splitCsv(sep: '\t')
               .flatMap { it ->
                    files = it.tail().findAll{ c -> c?.trim() };
                    return files.collect{ f -> [familyId: it.first(), size: files.size(), file: f ]}; 
                }.multiMap { it ->
                    sizes: tuple(it.familyId, it.size)
                    files: tuple(it.familyId, file("${it.file}*"))
                }                
}

include { variantRecalibratorSNP    } from './modules/vqsr'
include { variantRecalibratorIndel  } from './modules/vqsr'
include { applyVQSRSNP              } from './modules/vqsr'
include { applyVQSRIndel            } from './modules/vqsr'
include { splitMultiAllelics        } from './modules/vep'
include { vep                       } from './modules/vep'
include { tabix                     } from './modules/vep'


workflow {
    referenceGenome = file(params.referenceGenome)
    broad = file(params.broad)
    vepCache = file(params.vepCache)
    file(params.outputDir).mkdirs()

    sampleChannel().set{ familiesSizeFile }

    filtered = excludeMNPs(familiesSizeFile.files)
                    .join(familiesSizeFile.sizes)
                    .map{familyId, files, size -> tuple( groupKey(familyId, size), files)}
                    .groupTuple()
                    .map{ familyId, files -> tuple(familyId, files.flatten())}
    
    
    importGVCF(filtered, referenceGenome,broad)

    vcf = genotypeGVCF(importGVCF.out, referenceGenome)

    v = variantRecalibratorSNP(vcf, referenceGenome, broad).join(vcf)
    asnp = applyVQSRSNP(v) 

    indel = variantRecalibratorIndel(vcf, referenceGenome, broad).join(asnp)
    aindel = applyVQSRIndel(indel)

    s = splitMultiAllelics(aindel, referenceGenome) 

    vep(s, referenceGenome, vepCache) 
    tabix(vep.out) 

}