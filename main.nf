

nextflow.enable.dsl = 2

/**
Exclude MNPs detected with bedtools
*/

process excludeMNPs{
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
}


process importGVCF {

    label 'medium'

    input:
    tuple val(familyId), path(gvcfFiles),val(interval)
    path referenceGenome
    path broadResource

    output:
    tuple val(familyId), path("genomicsdb*"),val(interval)

    script:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("vcf.gz") }.collect { "-V $it" }.join(' ')

    """
    echo $familyId > file
    gatk -version
    gatk --java-options "-Xmx5g -Xms5g" GenomicsDBImport $exactGvcfFiles --genomicsdb-workspace-path genomicsdb_${interval}  --genomicsdb-shared-posixfs-optimizations true --bypass-feature-reader -L ${interval} 
    
    """       

}    

/**
Keep only SNP and Indel 
*/

process genotypeGVCF_multi {
    label 'geno'

    input:
    tuple val(familyId), path(genomicsdb),val(interval)
    path referenceGenome

    output:
    tuple val(familyId), path("*${interval}*genotyped.vcf.gz*")

    script:
    """
    echo $familyId > file
    gatk --java-options "-Xmx5g" GenotypeGVCFs -R $referenceGenome/${params.referenceGenomeFasta} --genomicsdb-shared-posixfs-optimizations true -V gendb://$genomicsdb -O ${familyId}_${interval}.genotyped.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation    
    """
}

process genotypeGVCF_solo {
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
    gatk --java-options "-Xmx8g" GenotypeGVCFs -R $referenceGenome/${params.referenceGenomeFasta} -V $exactGvcfFile -O ${familyId}.genotyped.vcf.gz
    """
}

process gatherVCF {
    label 'medium'

    input:
    tuple val(familyId), path(gvcfFiles)

    output:
    tuple val(familyId), path("*gather.genotyped.vcf.gz*")

    script:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("vcf.gz") }.collect { "$it" }.join("\n")
    // gatk --java-options "-Xmx8g" GatherVcfs $exactGvcfFiles -O ${familyId}.gather.genotyped.vcf.gz
    """
    echo $familyId > file
    echo "${exactGvcfFiles}" > list_file.txt
    bcftools concat -f list_file.txt | bcftools sort - -Oz -o ${familyId}.gather.genotyped.vcf.gz
    bcftools index -t  ${familyId}.gather.genotyped.vcf.gz
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

    
    filtered_one = filtered.filter{it[1].size() == 2}
    filtered_mult = filtered.filter{it[1].size() > 2}

    interval = Channel.fromPath("${params.intervalsFile}")
            .splitText().map{it -> it.trim()}
    
    // Run Combine & GenotypeGVCF on family with > 1 files
    importGVCF(filtered_mult.combine(interval), referenceGenome,broad)
                            
    vcf_fam = genotypeGVCF_multi(importGVCF.out, referenceGenome)
                        .combine(interval.count())
                        .map{familyId, files, size -> tuple( groupKey(familyId, size), files)}
                        .groupTuple()
                        .map{ familyId, files -> tuple(familyId, files.flatten())}
                        //.map{ familyId, files -> tuple( familyId, files.toSortedList { it -> (it.name =~ /*_chr(.*?)/)[0][1].toInteger() })}

    vcf_gather = gatherVCF(vcf_fam)
    
     // Run GenotypeGVCF on family with == 1 files
    vcf_one = genotypeGVCF_solo(filtered_one, referenceGenome)
    
    vcf = vcf_gather.concat(vcf_one)
    
    v = variantRecalibratorSNP(vcf, referenceGenome, broad).join(vcf)
    asnp = applyVQSRSNP(v) 

    indel = variantRecalibratorIndel(vcf, referenceGenome, broad).join(asnp)
    aindel = applyVQSRIndel(indel)

    s = splitMultiAllelics(aindel, referenceGenome) 

    vep(s, referenceGenome, vepCache) 
    tabix(vep.out) 

}