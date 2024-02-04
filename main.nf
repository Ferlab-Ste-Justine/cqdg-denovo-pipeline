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
    """
    set -e
    echo $familyId > file
    bcftools view --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY ${exactGvcfFile} -O z -o ${familyId}.${uuid}.filtered.gvcf.gz
    bcftools index -t ${familyId}.${uuid}.filtered.gvcf.gz
    """
}


process combineGVCF {
    label 'medium'

    container 'broadinstitute/gatk'
    input:
    tuple val(familyId), path(gvcfFiles)
    path referenceGenome

    output:
    tuple val(familyId), path("*combined.gvcf.gz*")

    script:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("gvcf.gz") }.collect { "--variant $it" }.join(' ')

    """
    echo $familyId > file
    gatk CombineGVCFs -R $referenceGenome/${params.referenceGenomeFasta} $exactGvcfFiles -O ${familyId}.combined.gvcf.gz
    """    

}    

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

process splitMultiAllelics{
    label 'medium'

    container 'staphb/bcftools'
    
    input:
    tuple val(familyId), path(vcfFile)
    path referenceGenome

    output:
    tuple val(familyId), path("*splitted.vcf*")

    script:
    def exactVcfFile = vcfFile.find { it.name.endsWith("vcf.gz") }
    """
    set -e
    echo $familyId > file
    bcftools annotate -x FORMAT/PRI ${exactVcfFile} | bcftools norm --write-index -c w -m -any -f $referenceGenome/${params.referenceGenomeFasta} --old-rec-tag OLD_RECORD --output-type z --output ${familyId}.normed.vcf.gz  
    bcftools view --min-ac 1 --output-type z --output ${familyId}.splitted.vcf.gz ${familyId}.normed.vcf.gz
    bcftools index -t ${familyId}.splitted.vcf.gz
    """
}

process vep {
    label 'vep'
    container 'ensemblorg/ensembl-vep'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    tuple val(familyId), path(vcfFile)
    path referenceGenome
    path vepCache

    output:
    path "*vep.vcf.gz"

    script:
    def exactVcfFile = vcfFile.find { it.name.endsWith("vcf.gz") }
    """
    vep \
    --fork ${params.vepCpu} \
    --dir ${vepCache} \
    --offline \
    --cache \
    --fasta $referenceGenome/${params.referenceGenomeFasta} \
    --input_file $exactVcfFile \
    --format vcf \
    --vcf \
    --output_file variants.${familyId}.vep.vcf.gz \
    --compress_output bgzip \
    --xref_refseq \
    --variant_class \
    --numbers \
    --hgvs \
    --hgvsg \
    --canonical \
    --symbol \
    --flag_pick \
    --fields "Allele,Consequence,IMPACT,SYMBOL,Feature_type,Gene,PICK,Feature,EXON,BIOTYPE,INTRON,HGVSc,HGVSp,STRAND,CDS_position,cDNA_position,Protein_position,Amino_acids,Codons,VARIANT_CLASS,HGVSg,CANONICAL,RefSeq" \
    --no_stats
    """
}

process tabix {
    label 'tiny'
    container 'staphb/htslib'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    path vcfFile

    output:
    path "*.tbi"

    script:
    """
    tabix $vcfFile
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

workflow {
    referenceGenome = file(params.referenceGenome)
    vepCache = file(params.vepCache)
    file(params.outputDir).mkdirs()

    sampleChannel().set{ familiesSizeFile }
    
    filtered = excludeMNPs(familiesSizeFile.files)
                    .join(familiesSizeFile.sizes)
                    .map{familyId, files, size -> tuple( groupKey(familyId, size), files)}
                    .groupTuple()
                    .map{ familyId, files -> tuple(familyId, files.flatten())}

    combineGVCF(filtered, referenceGenome) 
    genotypeGVCF(combineGVCF.out, referenceGenome) 
    splitMultiAllelics(genotypeGVCF.out, referenceGenome) 
    vep(splitMultiAllelics.out, referenceGenome, vepCache) 
    tabix(vep.out) 

}