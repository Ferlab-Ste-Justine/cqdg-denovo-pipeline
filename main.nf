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
    // chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
    """
    set -e
    echo $familyId > file
    bcftools view --regions chr1 ${exactGvcfFile} -O z -o ${familyId}.${uuid}.filtered.gvcf.gz
    bcftools index -t ${familyId}.${uuid}.filtered.gvcf.gz
    """
}


process importGVCF {
    label 'medium'

    container 'broadinstitute/gatk'
    input:
    tuple val(familyId), path(gvcfFiles)
    path referenceGenome

    output:
    tuple val(familyId), path("genomicsdb")

    script:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("gvcf.gz") }.collect { "-V $it" }.join(' ')

    """
    echo $familyId > file
    gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport $exactGvcfFiles --genomicsdb-workspace-path genomicsdb -L chr1
    """    

}    

process genotypeGVCF {
    label 'geno'

    container 'broadinstitute/gatk'

    input:
    tuple val(familyId), path(genomicsdb)
    path referenceGenome

    output:
    tuple val(familyId), path("*genotyped.vcf.gz*")

    script:
    def workspace = genomicsdb.getBaseName()
    """
    echo $familyId > file
    gatk --java-options "-Xmx24g" GenotypeGVCFs -R $referenceGenome/${params.referenceGenomeFasta} -V gendb://$genomicsdb -O ${familyId}.genotyped.vcf.gz -G StandardAnnotation -G AS_StandardAnnotation -L chr1
    """

}

process variantRecalibratorSNP {
    label 'geno'

    container 'broadinstitute/gatk'

    input:
    tuple val(prefix), path(vcf)
    path referenceGenome
    path broadResource

    output:
    tuple val(prefix), path("*.recal*"), path("*.tranches")
    
    script:
    def tranches = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0"].collect{"-tranche $it"}.join(' ')
    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
    $tranches \
    -R $referenceGenome/${params.referenceGenomeFasta} \
    -V ${vcf.first()} \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
    ${broadResource}/hapmap_3.3.hg38.vcf.gz  \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 \
    ${broadResource}/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 \
    ${broadResource}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR  -an DP \
    -mode SNP -O ${prefix}.recal --tranches-file ${prefix}.tranches
    """

    

}

process applyVQSRSNP {
    label 'geno'

    container 'broadinstitute/gatk'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    tuple val(prefix), path(recal), path(tranches), path(vcf)

    output:
    path("*.snp.vqsr.vcf.gz*")

    script:
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def exactRecal = recal.find { it.name.endsWith("recal") }
    """
    gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
    -V ${exactVcfFile} \
    --recal-file ${exactRecal} \
    -mode SNP \
    --tranches-file ${tranches} \
    --truth-sensitivity-filter-level 99.9 \
    --create-output-variant-index true \
    -O ${prefix}.snp.vqsr.vcf.gz
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
    broad = file(params.broad)
    vepCache = file(params.vepCache)
    file(params.outputDir).mkdirs()

    sampleChannel().set{ familiesSizeFile }
    
    filtered = excludeMNPs(familiesSizeFile.files)
                    .join(familiesSizeFile.sizes)
                    .map{familyId, files, size -> tuple( groupKey(familyId, size), files)}
                    .groupTuple()
                    .map{ familyId, files -> tuple(familyId, files.flatten())}

    filtered | view
    importGVCF(filtered, referenceGenome) | view
    genotypeGVCF(importGVCF.out, referenceGenome) | view 
    v = variantRecalibratorSNP(genotypeGVCF.out, referenceGenome, broad).join(genotypeGVCF.out)
    v | view
    applyVQSRSNP(v) | view

    // splitMultiAllelics(genotypeGVCF.out, referenceGenome) | view 
    // vep(splitMultiAllelics.out, referenceGenome, vepCache) 
    // tabix(vep.out) 

}