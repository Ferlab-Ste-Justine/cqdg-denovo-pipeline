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
    bcftools annotate -x FORMAT/PRI ${exactVcfFile} | bcftools norm -c w -m -any -f $referenceGenome/${params.referenceGenomeFasta} --old-rec-tag OLD_RECORD --output-type z --output ${familyId}.normed.vcf.gz  
    bcftools view --min-ac 1 --output-type z --output ${familyId}.splitted.vcf.gz ${familyId}.normed.vcf.gz
    bcftools index -t ${familyId}.splitted.vcf.gz
    """

    stub:
    """
    touch ${familyId}.splitted.vcf.gz
    touch ${familyId}.splitted.vcf.gz.tbi
    """
}

process vep {
    label 'vep'

    publishDir "${params.outputDir}", mode: 'copy'

    input:
    tuple val(familyId), path(vcfFile)
    path referenceGenome
    path vepCache

    output:
    tuple val(familyId), path("*vep.vcf.gz")

    script:
    def args = task.ext.args ?: ''
    def exactVcfFile = vcfFile.find { it.name.endsWith("vcf.gz") }
    """
    vep \\
        --fork ${params.vepCpu} \\
        --dir ${vepCache} \\
        --offline \\
        --cache \\
        --fasta $referenceGenome/${params.referenceGenomeFasta} \\
        --input_file $exactVcfFile \\
        --format vcf \\
        --vcf \\
        --output_file variants.${familyId}.vep.vcf.gz \\
        --compress_output bgzip \\
        --xref_refseq \\
        --variant_class \\
        --numbers \\
        --hgvs \\
        --hgvsg \\
        --canonical \\
        --symbol \\
        --flag_pick \\
        --fields "Allele,Consequence,IMPACT,SYMBOL,Feature_type,Gene,PICK,Feature,EXON,BIOTYPE,INTRON,HGVSc,HGVSp,STRAND,CDS_position,cDNA_position,Protein_position,Amino_acids,Codons,VARIANT_CLASS,HGVSg,CANONICAL,RefSeq" \\
        --no_stats \\
        $args
    """

    stub:
    """
    touch variants.${familyId}.vep.vcf.gz
    """
}

process tabix {
    label 'tiny'

    publishDir "${params.outputDir}", mode: 'copy'

    input:
    tuple val(familyId), path(vcfFile)

    output:
    path "*.tbi"

    script:
    def args = task.ext.args ?: ''

    """
    tabix \\
        $vcfFile \\
        $args
    """
    script:
    """
    touch ${vcfFile}.tbi
    """

} 
