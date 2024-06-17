process hardFiltering {
    label 'medium'

    input:
    tuple val(prefixId), path(vcf)
    val(filters)

    output:
    tuple val(prefixId), path("*.hardfilter.vcf.gz*")

    script:
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def filterOptions = filters.collect{ "-filter \"${it.expression}\" --filter-name \"${it.name}\"" }.join(" ")
    """
    set -e

    gatk VariantFiltration \
     -V ${exactVcfFile} \
    ${filterOptions} \
     -O  ${prefixId}.hardfilter.vcf.gz
    """
    stub:
    """
    touch ${prefixId}.hardfilter.vcf.gz
    """
}