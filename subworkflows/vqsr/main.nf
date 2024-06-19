include { variantRecalibratorIndel; variantRecalibratorSNP; applyVQSRIndel; applyVQSRSNP} from '../../modules/vqsr'

/**
Filter out probable artifacts from the callset using the Variant Quality Score Recalibration (VQSR) procedure

The input and output formats are the same:
    Input: (prefixId,  [some.file.vcf.gz, some.file.vcf.gz.tbi])

All output files will be prefixed with the given prefixId.
*/
workflow VQSR {
    take:
        input // channel: (val(prefixId),  [.vcf.gz, .vcf.gz.tbi])

    main:
        referenceGenome = file(params.referenceGenome)
        broad = file(params.broad)
        

        outputSNP = variantRecalibratorSNP(input, referenceGenome, broad)
            | join(input)
            | applyVQSRSNP 

        output = variantRecalibratorIndel(input, referenceGenome, broad)
            | join(outputSNP)
            | applyVQSRIndel
  
    emit:
        output // channel: (val(prefixId),  [.vcf.gz, .vcf.gz.tbi])
}