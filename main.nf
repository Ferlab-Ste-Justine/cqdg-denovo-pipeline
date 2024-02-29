/**
Exclude MNPs detected with bedtools
*/

process excludeMNPs{
    label 'medium'

    container 'staphb/bcftools'
    
    input:
    tuple val(familyId), path(gvcfFile)

    output:
    tuple val(familyId), path("*filtered.gvcf.gz*")

    script:
    def exactGvcfFile = gvcfFile.find { it.name.endsWith("gvcf.gz") }
    def uuid = UUID.randomUUID().toString()
    // --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
    """
    set -e
    echo $familyId > file
    bcftools view --exclude-type mnps  ${exactGvcfFile} -O z -o ${familyId}.${uuid}.filtered.gvcf.gz
    bcftools index -t ${familyId}.${uuid}.filtered.gvcf.gz
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
    tuple val(familyId), path("genomicsdb")

    script:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("gvcf.gz") }.collect { "-V $it" }.join(' ')

    """
    echo $familyId > file
 
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport $exactGvcfFiles --genomicsdb-workspace-path genomicsdb -L ${broadResource}/${params.intervalsFile}  """    

}    

/**
Keep only SNP and Indel 
*/

process genotypeGVCF {
    label 'geno'

    container 'broadinstitute/gatk'

    input:
    tuple val(familyId), path(genomicsdb)
    path referenceGenome
    path broadResource

    output:
    tuple val(familyId), path("*genotyped.vcf.gz*")

    script:
    def workspace = genomicsdb.getBaseName()
    """
    echo $familyId > file
    gatk --java-options "-Xmx24g" GenotypeGVCFs -R $referenceGenome/${params.referenceGenomeFasta} -V gendb://$genomicsdb -O ${familyId}.genotyped.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation
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
    
    // filtered | view

    importGVCF(filtered, referenceGenome,broad)

    vcf = genotypeGVCF(importGVCF.out, referenceGenome,broad)
    // vcf | view
    

    v = variantRecalibratorSNP(vcf, referenceGenome, broad).join(vcf)
    // v | view
    asnp = applyVQSRSNP(v) 
    // asnp | view

    indel = variantRecalibratorIndel(vcf, referenceGenome, broad).join(asnp)
    // indel | view
    aindel = applyVQSRIndel(indel)
    // aindel | view

    s = splitMultiAllelics(aindel, referenceGenome) 
    // s | view 

    vep(s, referenceGenome, vepCache) 
    tabix(vep.out) 

}