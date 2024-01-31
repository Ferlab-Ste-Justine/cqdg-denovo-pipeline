process combineGVCF {
    label 'medium'

    container 'broadinstitute/gatk'
    input:
    val familyId
    path gvcfFiles
    path referenceGenome

    output:
    path "*combined.gvcf.gz*"

    script:
    def exactGvcfFiles = gvcfFiles.findAll { it.name.endsWith("gvcf.gz") }.collect { "--variant $it" }.join(' ')

    """
    gatk CombineGVCFs -R $referenceGenome/Homo_sapiens_assembly38.fasta $exactGvcfFiles -O ${familyId}.combined.gvcf.gz
    """    

}    

process genotypeGVCF {
    label 'geno'

    container 'broadinstitute/gatk'

    input:
    val familyId
    path gvcfFile
    path referenceGenome

    output:
    path "*genotyped.vcf.gz*"

    script:
    def exactGvcfFile = gvcfFile.find { it.name.endsWith("gvcf.gz") }
    """
    gatk --java-options "-Xmx24g" GenotypeGVCFs -R $referenceGenome/Homo_sapiens_assembly38.fasta -V $exactGvcfFile -O ${familyId}.genotyped.vcf.gz
    """

}

process splitMultiAllelics{
    label 'medium'

    container 'broadinstitute/gatk:4.1.4.1'
    
    input:
    val familyId
    path vcfFile
    path referenceGenome

    output:
    path "*splitted.vcf*"

    script:
    def exactVcfFile = vcfFile.find { it.name.endsWith("vcf.gz") }
    """
    gatk LeftAlignAndTrimVariants -R $referenceGenome/Homo_sapiens_assembly38.fasta -V $exactVcfFile -O ${familyId}.splitted.vcf.gz --split-multi-allelics
    """
}

process vep {
    label 'vep'

    container 'ensemblorg/ensembl-vep'
    
    input:
    val familyId
    path vcfFile
    path referenceGenome
    path vepCache

    output:
    path "*vep.vcf.gz"

    script:
    def exactVcfFile = vcfFile.find { it.name.endsWith("vcf.gz") }
    """
    vep \
    --fork ${params.vep_cpu} \
    --dir ${vepCache} \
    --offline \
    --cache \
    --fasta $referenceGenome/Homo_sapiens_assembly38.fasta \
    --input_file $exactVcfFile \
    --format vcf \
    --vcf \
    --output_file ${familyId}.vep.vcf.gz \
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
    label 'small'

    container 'staphb/htslib'

    input:
    path vcfFile

    output:
    path "*.tbi"

    script:
    """
    tabix $vcfFile
    """

}

process copyFinalDestination {
    label 'small'

    input:
    path destination
    path vcfFile
    path tbiFile

    output:
    stdout

    script:
    """
    cp -f $vcfFile $destination/
    cp -f $tbiFile $destination/
    """

}    



workflow {

    families = Channel.fromPath(file("$params.inputDir/samples.txt"))
               .splitCsv(sep: '\t')
               .map { it[0]}

    gvcfs = Channel.fromPath(file("$params.inputDir/samples.txt"))
               .splitCsv(sep: '\t')
               .map { it.tail().collect{ e -> file("${params.inputDir}/${e}*")}.flatten()}       

    referenceGenome = file(params.referenceGenome)
    vepCache = file("s3://cqdg-prod-file-import/nextflow/vep")
    finalDestination = file("s3://cqdg-prod-file-import/nextflow/output")

    finalDestination.mkdirs()
    families | view
    gvcfs | view
    combineGVCF(families, gvcfs, referenceGenome) | view
    genotypeGVCF(families, combineGVCF.out, referenceGenome) | view
    splitMultiAllelics(families, genotypeGVCF.out, referenceGenome) | view
    vep(families, splitMultiAllelics.out, referenceGenome, vepCache) | view
    tabix(vep.out) | view

    copyFinalDestination(finalDestination, vep.out, tabix.out)


}