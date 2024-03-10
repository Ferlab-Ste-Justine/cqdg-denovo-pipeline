process denovoCNN{
    label 'dnn'

    container 'gelana/denovocnn:1.0'
    publishDir "s3://cqdg-prod-app-datalake/datamart/denovo_dee/denovoCNN/", mode: 'copy'
    
    input:
    tuple val(familyId), path(variantList), path(childCram), path(fatherCram), path(motherCram)
    path referenceGenome

    output:
    tuple val(familyId), path("*predictions.csv")

    script:
    def exactChildCram = childCram.find { it.name.endsWith(".cram") }
    def exactFatherCram = fatherCram.find { it.name.endsWith(".cram") }
    def exactMotherCram = motherCram.find { it.name.endsWith(".cram") }
    """
    set -e
    echo $familyId > file
    /app/apply_denovocnn.sh\
        --workdir=\$(pwd) \
        --variant-list=$variantList \
        --child-bam=$exactChildCram \
        --father-bam=$exactFatherCram \
        --mother-bam=$exactMotherCram \
        --snp-model=/app/models/snp \
        --in-model=/app/models/ins \
        --del-model=/app/models/del \
        --genome=$referenceGenome/${params.referenceGenomeFasta} \
        --output=${familyId}.predictions.csv    
    """
}

workflow {
    referenceGenome = file(params.referenceGenome)
    //variantList = file("s3://cqdg-prod-app-datalake/datamart/denovo_dee/denovoCNN/variant_list/chr1/part-00000-797be2de-9892-400f-a7a3-4532a4490bd4-c000.csv")
    //variantList = file("s3://cqdg-prod-app-datalake/datamart/denovo_dee/denovoCNN/variant_list/hail_denovo.tsv")
    variantList = file("s3://cqdg-prod-app-datalake/datamart/denovo_dee/denovoCNN/variant_list/rqdm_denovo/part-00000-814a22aa-623a-465b-9e9a-6035d6244e3f-c000.csv")
    childCram = file("s3://cqdg-prod-file-import/jmichaud/DEE/dataset_default/2219/S16943.cram*")
    fatherCram = file("s3://cqdg-prod-file-import/jmichaud/DEE/dataset_default/2219/S16944.cram*")
    motherCram = file("s3://cqdg-prod-file-import/jmichaud/DEE/dataset_default/2219/S16945.cram*")
    
    denovoCNN(tuple("HSC0047_rqdm", variantList, childCram, fatherCram, motherCram), referenceGenome) | view

    // splitMultiAllelics(genotypeGVCF.out, referenceGenome) | view 
    // vep(splitMultiAllelics.out, referenceGenome, vepCache) 
    // tabix(vep.out) 

}