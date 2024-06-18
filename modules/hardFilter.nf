process hardFiltering {
    label 'medium'

    input:
    tuple val(prefixId), path(vcf)
    val(filters)

    output:
    tuple val(prefixId), path("*.hardfilter.vcf.gz*")

    script:
    def args = task.ext.args ?: ''
    def argsjava = task.ext.args ?: ''
    def exactVcfFile = vcf.find { it.name.endsWith("vcf.gz") }
    def filterOptions = filters.collect{ "-filter \"${it.expression}\" --filter-name \"${it.name}\"" }.join(" ")

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK VariantFiltration] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    set -e

    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData $argsjava" \\
        VariantFiltration \\
        -V ${exactVcfFile} \\
        ${filterOptions} \\
        -O  ${prefixId}.hardfilter.vcf.gz \\
        $args
    """
    stub:
    """
    touch ${prefixId}.hardfilter.vcf.gz
    """
}