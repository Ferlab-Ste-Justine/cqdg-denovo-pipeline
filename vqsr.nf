
process variantRecalibratorSNP {
    label 'medium'

    container 'broadinstitute/gatk'

    input:
    tuple val(prefix), path(vcf)
    path referenceGenome
    path broadResource

    output:
    tuple val(prefix), path("*.recal*"), path("*.tranches")
    
    script:
    def tranches = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0"].collect{"-tranche $it"}.join(' ')
    def annotationValues = ["QD","MQRankSum","ReadPosRankSum","FS","MQ","SOR","DP"].collect{"-an $it"}.join(' ')
    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
    $tranches \
    --trust-all-polymorphic \
    -R $referenceGenome/${params.referenceGenomeFasta} \
    -V ${vcf.first()} \
    --resource:hapmap,known=false,training=true,truth=true,prior=15 \
    ${broadResource}/hapmap_3.3.hg38.vcf.gz  \
    --resource:omni,known=false,training=true,truth=false,prior=12 \
    ${broadResource}/1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10 \
    ${broadResource}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${broadResource}/Homo_sapiens_assembly38.dbsnp138.vcf \
    $annotationValues \
    --max-gaussians 6 \
    -mode SNP -O ${prefix}.recal --tranches-file ${prefix}.tranches
    """

    

}

process variantRecalibratorIndel {
    label 'medium'

    container 'broadinstitute/gatk'

    input:
    tuple val(prefix), path(vcf)
    path referenceGenome
    path broadResource

    output:
    tuple val(prefix), path("*.recal*"), path("*.tranches")
    
    script:
    def tranches = ["100.0","99.95","99.9","99.5","99.0","97.0","96.0","95.0","94.0","93.5","93.0","92.0","91.0","90.0"].collect{"-tranche $it"}.join(' ')
    def annotationValues = ["FS","ReadPosRankSum","MQRankSum","QD","SOR","DP"].collect{"-an $it"}.join(' ')
    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
    $tranches \
    -R $referenceGenome/${params.referenceGenomeFasta} \
    -V ${vcf.first()} \
    --trust-all-polymorphic \
    --resource:mills,known=false,training=true,truth=true,prior=12 ${broadResource}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
    --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${broadResource}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${broadResource}/Homo_sapiens_assembly38.dbsnp138.vcf \
    $annotationValues \
    --max-gaussians 4 \
    -mode INDEL -O ${prefix}.recal --tranches-file ${prefix}.tranches
    """

    

}

process applyVQSRSNP {
    label 'medium'

    container 'broadinstitute/gatk'

    input:
    tuple val(prefix), path(recal), path(tranches), path(vcf)

    output:
    tuple val(prefix), path("*.snp.vqsr.vcf.gz*")

    script:
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def exactRecal = recal.find { it.name.endsWith("recal") }
    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xms5G -Xmx5G -XX:ParallelGCThreads=2" ApplyVQSR \
    -V ${exactVcfFile} \
    --recal-file ${exactRecal} \
    -mode SNP \
    --tranches-file ${tranches} \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -O ${prefix}.snp.vqsr.vcf.gz
    """

}


process applyVQSRIndel {
    label 'medium'

    container 'broadinstitute/gatk'
    publishDir "${params.outputDir}", mode: 'copy'

    input:
    tuple val(prefix), path(recal), path(tranches), path(vcf)

    output:
    tuple val(prefix), path("*.vcf.gz*")

    script:
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def exactRecal = recal.find { it.name.endsWith("recal") }
    """
    echo $prefix > file
    gatk --java-options "-Xms5G -Xmx5G -XX:ParallelGCThreads=2" ApplyVQSR \
    -V ${exactVcfFile} \
    --recal-file ${exactRecal} \
    -mode INDEL \
    --tranches-file ${tranches} \
    --truth-sensitivity-filter-level 99.7 \
    --create-output-variant-index true \
    -O ${prefix}.vcf.gz
    """

}


workflow {
    referenceGenome = file(params.referenceGenome)
    broad = file(params.broad)
    Channel.fromFilePairs("$params.inputDir/*.vcf.gz{,.tbi}")
        .tap{ vcfs}


    snp = variantRecalibratorSNP(vcfs, referenceGenome, broad)
        .join(vcfs)

    indel = variantRecalibratorIndel(vcfs, referenceGenome, broad)      

    snp | view
        
    asnp = applyVQSRSNP(snp) 
    aindel = applyVQSRIndel(indel.join(asnp))


}