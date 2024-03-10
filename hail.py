import hail as hl
gnomad = hl.import_table('s3a://cqdg-prod-app-datalake/datamart/denovo_dee/gnomad_v3/*')
gnomad = gnomad.filter(gnomad.af != '""')
gnomad = gnomad.select(variant='chr' + gnomad.chromosome + ':' + gnomad.start + ':' + gnomad.reference + ':' + gnomad.alternate, af=hl.float64(gnomad.af))
priors = gnomad.transmute(**hl.parse_variant(gnomad.variant, reference_genome='GRCh38')).key_by('locus', 'alleles')
priors.write('gnomad.ht')
priors = hl.read_table('gnomad.ht')
v = hl.import_vcf("variants.HSC0047.vep.vcf.gz", reference_genome='GRCh38', force_bgz=True, array_elements_required=False)
trio = hl.Trio(s='S16943', fam_id='S16943', is_female=True, mat_id='S16945', pat_id='S16944')
ped = hl.Pedigree([trio])
de_novo_results = hl.de_novo(v, ped, pop_frequency_prior=priors[v.row_key].af)
de_novo_results.aggregate(hl.agg.counter(de_novo_results.confidence))

 