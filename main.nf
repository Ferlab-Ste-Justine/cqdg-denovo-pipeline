referenceGenome = file(params.referenceGenome)
vepCache = file(params.vepCache)
file(params.finalDestination).mkdirs()

process combineGVCF {
    label 'medium'

    container 'broadinstitute/gatk'
    input:
    tuple val(familyId), path(gvcfFiles)

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

    container 'broadinstitute/gatk:4.1.4.1'
    
    input:
    tuple val(familyId), path(vcfFile)

    output:
    tuple val(familyId), path("*splitted.vcf*")

    script:
    def exactVcfFile = vcfFile.find { it.name.endsWith("vcf.gz") }
    """
    echo $familyId > file
    gatk LeftAlignAndTrimVariants -R $referenceGenome/${params.referenceGenomeFasta} -V $exactVcfFile -O ${familyId}.splitted.vcf.gz --split-multi-allelics
    """
}

process vep {
    label 'vep'

    container 'ensemblorg/ensembl-vep'
    
    publishDir: "${params.finalDestination}", mode: 'copy'

    input:
    tuple val(familyId), path(vcfFile)

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

    publishDir: "${params.finalDestination}", mode: 'copy'

    input:
    path vcfFile

    output:
    path "*.tbi"

    script:
    """
    tabix $vcfFile
    """

} 

def sampleChannel() {
   return Channel.fromPath(file("$params.sampleFile"))
               .splitCsv(sep: '\t')
               .map { tuple(it.first(), it.tail().collect{ f -> file("${f}*")}.flatten()) };
}

workflow {
    
    familiesAndFiles = sampleChannel()

    familiesAndFiles | view
    combineGVCF(familiesAndFiles) | view
    genotypeGVCF(combineGVCF.out) | view
    splitMultiAllelics(genotypeGVCF.out) | view
    vep(splitMultiAllelics.out) | view
    tabix(vep.out) | view



}