params.Ref_Abbr = "DV10"
params.rawVCF = "${baseDir}/VCF/Samples_minicore_1-135_whole_genome.DV10.raw.vcf.gz"
params.referencePrefix = "DarmorV10_Chromosomes_Only"
params.VCF = "${baseDir}/VCF"
params.GoodCoord = "${baseDir}/data/VQSR_training_good.txt"
params.MediumCoord = "${baseDir}/data/VQSR_training_medium.txt"
params.ref = "${baseDir}/Genome/Reference"
params.KnownSites = "${baseDir}/data/Truth_Set.DV10.SNPs.vcf.gz"
params.outputVCFprefix = "Samples_minicore_1-135_whole_genome_called"

process extract_training_sets {
    container 'oras://community.wave.seqera.io/library/bcftools:1.15--6a7920e74c21902b'
    publishDir "${params.VCF}", mode: 'copy'

    input:
    path(good_coord)
    path(medium_coord)


    output:
    tuple path("VQSR_Training_Good.${params.Ref_Abbr}.vcf"),
          path("VQSR_Training_Medium.${params.Ref_Abbr}.vcf")

    script:
    """
    echo ""
    echo "\$(date) - Extracting Good SNPs"
    echo ""

    bcftools view -T ${good_coord} "${params.rawVCF}" -Oz -o "VQSR_Training_Good.${params.Ref_Abbr}.vcf"
    
    echo ""
    echo "\$(date) - DONE!"
    echo ""

    echo ""
    echo "\$(date) - Extracting Medium SNPs"
    echo ""

    bcftools view -T ${medium_coord} "${params.rawVCF}" -Oz -o "VQSR_Training_Medium.${params.Ref_Abbr}.vcf"

    echo ""
    echo "\$(date) - DONE!"
    echo ""
    """
}

process run_VQSR {
    container 'community.wave.seqera.io/library/bwa_fastp_gatk4_samtools:b552bdcea7a3515f'
    publishDir "${params.VCF}", mode: 'copy'

    input:
    tuple path(raw_SNPs), path(good_SNPs), path(medium_SNPs), path(ref), path(known)
    path vcf_index
    path ref_index

    output:
    tuple path("output.SNP.${params.Ref_Abbr}.recal"),
          path("output.SNP.${params.Ref_Abbr}.tranches"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.5.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.100.0.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.9.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.0.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.98.5.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.98.0.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.95.0.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.90.0.tranche.vcf.gz"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.5.tranche.vcf.gz.tbi"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.100.0.tranche.vcf.gz.tbi"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.9.tranche.vcf.gz.tbi"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.0.tranche.vcf.gz.tbi"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.98.5.tranche.vcf.gz.tbi"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.98.0.tranche.vcf.gz.tbi"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.95.0.tranche.vcf.gz.tbi"),
          path("${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.90.0.tranche.vcf.gz.tbi")

    script:
    """
    echo "Prepareing the VCFs"
    gatk IndexFeatureFile -I ${known}
    gatk IndexFeatureFile -I ${good_SNPs}
    gatk IndexFeatureFile -I ${medium_SNPs}

    echo "\$(date) - Recalibrating Variants"

    gatk --java-options "-Xmx100g -Xms100g" VariantRecalibrator -R "${ref}" -V "${raw_SNPs}" --trust-all-polymorphic -mode SNP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum --resource:gmm_good,known=false,training=true,truth=true,prior=15.0 "${good_SNPs}" --resource:gmm_medium,known=false,training=true,truth=false,prior=10.0 "${medium_SNPs}" --resource:known_sites,known=true,training=false,truth=false,prior=2.0 "${known}" -O "output.SNP.${params.Ref_Abbr}.recal" --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --maximum-training-variants 2000000 --truth-sensitivity-tranche 100.0 --truth-sensitivity-tranche 99.9 --truth-sensitivity-tranche 99.5 --truth-sensitivity-tranche 99.0 --truth-sensitivity-tranche 98.5 --truth-sensitivity-tranche 98.0 --truth-sensitivity-tranche 95.0 --truth-sensitivity-tranche 90.0


    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"
    echo ""
    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 99.5"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.5.tranche.vcf.gz" --truth-sensitivity-filter-level 99.5 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"
    echo ""
  
    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 100.0"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.100.0.tranche.vcf.gz" --truth-sensitivity-filter-level 100.0 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"

    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 99.9"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.9.tranche.vcf.gz" --truth-sensitivity-filter-level 99.9 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"

    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 99.0"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.99.0.tranche.vcf.gz" --truth-sensitivity-filter-level 99.0 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"

    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 98.5"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.98.5.tranche.vcf.gz" --truth-sensitivity-filter-level 98.5 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"

    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 98.0"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.98.0.tranche.vcf.gz" --truth-sensitivity-filter-level 98.0 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"

    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 95.0"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.95.0.tranche.vcf.gz" --truth-sensitivity-filter-level 95.0 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Recalibrating Variants: DONE!"

    echo "\$(date) - Apply VQSR - Truth Sensitivity Threshold: 90.0"
    echo ""

    gatk --java-options "-Xmx100g -Xms100g" ApplyVQSR -R "${ref}" -V "${raw_SNPs}" -O "${params.outputVCFprefix}.${params.Ref_Abbr}.SNPs.90.0.tranche.vcf.gz" --truth-sensitivity-filter-level 90.0 --tranches-file "output.SNP.${params.Ref_Abbr}.tranches" --recal-file "output.SNP.${params.Ref_Abbr}.recal" -mode SNP

    echo ""
    echo "\$(date) - Variant Score Recalibration Applied!"
    echo ""
    """
}

workflow {
    training_data = extract_training_sets("${params.GoodCoord}", "${params.MediumCoord}")
    vcf = file("${params.rawVCF}")
    vcf_index = file("${params.rawVCF}.tbi")
    ref = file("${params.ref}/${params.referencePrefix}.fa") 
    ref_index = Channel
        .fromPath("${params.ref}/${params.referencePrefix}.*{fai,dict}")
        .collect() 
    known = file("${params.KnownSites}")
    VQSR_inputs = training_data.map { good, medium ->
        tuple(vcf, good, medium, ref, known)
    }
    VQRS_outputs = run_VQSR(VQSR_inputs, vcf_index, ref_index) 
}
