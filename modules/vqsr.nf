
/**
Build a recalibration model to score SNP variant quality for filtering purposes

Note: pre-requisite step for applyVQSRSNP
*/
process variantRecalibratorSNP {
    label 'medium'


    input:
    tuple val(prefix), path(vcf)
    path referenceGenome
    path broadResource

    output:
    tuple val(prefix), path("*.recal*"), path("*.tranches")
    
    script:
    def args = task.ext.args ?: ''
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def tranches = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0"].collect{"-tranche $it"}.join(' ')
    def annotationValues = ["QD","MQRankSum","ReadPosRankSum","FS","MQ","SOR","DP"].collect{"-an $it"}.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantRecalibrator \\
        $tranches \\ 
        --trust-all-polymorphic \\
        -R $referenceGenome/${params.referenceGenomeFasta} \\
        -V ${exactVcfFile} \\
        --resource:hapmap,known=false,training=true,truth=true,prior=15 \\
        ${broadResource}/hapmap_3.3.hg38.vcf.gz  \\
        --resource:omni,known=false,training=true,truth=false,prior=12 \\
        ${broadResource}/1000G_omni2.5.hg38.vcf.gz \\
        --resource:1000G,known=false,training=true,truth=false,prior=10 \\
        ${broadResource}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=7 ${broadResource}/Homo_sapiens_assembly38.dbsnp138.vcf \\
        $annotationValues \\
        --max-gaussians 6 \\
        -mode SNP -O ${prefix}.recal --tranches-file ${prefix}.tranches \\
        $args
    """

    stub:
    """
    touch ${prefix}.recal
    touch ${prefix}.tranches
    """
}

/**
Build a recalibration model to score Indel variant quality for filtering purposes

Note: pre-requisite step for applyVQSRIndel
*/
process variantRecalibratorIndel {
    label 'medium'


    input:
    tuple val(prefix), path(vcf)
    path referenceGenome
    path broadResource

    output:
    tuple val(prefix), path("*.recal*"), path("*.tranches")
    
    script:
    def args = task.ext.args ?: ''
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def tranches = ["100.0","99.95","99.9","99.5","99.0","97.0","96.0","95.0","94.0"].collect{"-tranche $it"}.join(' ')
    def annotationValues = ["FS","ReadPosRankSum","MQRankSum","QD","SOR","DP"].collect{"-an $it"}.join(' ')

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantRecalibrator \\
        $tranches \\
        -R $referenceGenome/${params.referenceGenomeFasta} \\
        -V ${exactVcfFile} \\
        --trust-all-polymorphic \\
        --resource:mills,known=false,training=true,truth=true,prior=12 ${broadResource}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \\
        --resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${broadResource}/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \\
        --resource:dbsnp,known=true,training=false,truth=false,prior=2 ${broadResource}/Homo_sapiens_assembly38.dbsnp138.vcf \\
        $annotationValues \\
        --max-gaussians 4 \\
        -mode INDEL \\
        -O ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        $args
    """

    stub:
    """
    touch ${prefix}.recal
    touch ${prefix}.tranches
    """
}

/*
Apply a score cutoff to filter SNP variants based on a recalibration table
*/
process applyVQSRSNP {
    label 'medium'


    input:
    tuple val(prefix), path(recal), path(tranches), path(vcf)

    output:
    tuple val(prefix), path("*.snp.vqsr_${params.TSfilterSNP}.vcf.gz*")

    script:
    def args = task.ext.args ?: ''
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def exactRecal = recal.find { it.name.endsWith("recal") }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyVQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyVQSR \\
        -V ${exactVcfFile} \\
        --recal-file ${exactRecal} \\
        -mode SNP \\
        --tranches-file ${tranches} \\
        --truth-sensitivity-filter-level ${params.TSfilterSNP} \\
        --create-output-variant-index true \\
        -O ${prefix}.snp.vqsr_${params.TSfilterSNP}.vcf.gz \\
        $args
    """
    stub:
    """
    touch ${prefix}.snp.vqsr_${params.TSfilterSNP}.vcf.gz
    """

}

/*
Apply a score cutoff to filter Indel variants based on a recalibration table
*/
process applyVQSRIndel {
    label 'medium'

    container 'broadinstitute/gatk'

    input:
    tuple val(prefix), path(recal), path(tranches), path(vcf)

    output:
    tuple val(prefix), path("*.snpindel.vqsr_${params.TSfilterINDEL}.vcf.gz*")

    script:
    def args = task.ext.args ?: ''
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def exactRecal = recal.find { it.name.endsWith("recal") }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK ApplyVQSR] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }
    """
    set -e
    echo $prefix > file
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        ApplyVQSR \\
        -V ${exactVcfFile} \\
        --recal-file ${exactRecal} \\
        -mode INDEL \\
        --tranches-file ${tranches} \\
        --truth-sensitivity-filter-level ${params.TSfilterINDEL} \\
        --create-output-variant-index true \\
        -O ${prefix}.snpindel.vqsr_${params.TSfilterINDEL}.vcf.gz \\
        $args
    """
    stub:
    """
    touch ${prefix}.snpindel.vqsr_${params.TSfilterINDEL}.vcf.gz
    """
}