#!/usr/bin/env nextflow

// QA AND PREPROCESSING
include { BWA_INDEX } from './modules/bwa_index'
include { SAMTOOLS_QUICKCHECK } from './modules/samtools_quickcheck'
include { SAMTOOLS_CRAMTOFASTQ } from './modules/samtools_cramtofastq'
include { XENGSORT_INDEX } from './modules/xengsort_index'
include { XENGSORT_CLASSIFY } from './modules/xengsort_classify'
include { SEQKIT_STATS } from './modules/seqkit_stats'
include { BBMAP_REPAIR } from './modules/bbmap_repair'
include { TRIMMOMATIC_PE } from './modules/trimmomatic_pe'
include { BWA_MEM } from './modules/bwa_mem'
include { SAMTOOLS_FLAGSTAT; SAMTOOLS_COVERAGE } from './modules/samtools_assess'
include { SAMTOOLS_SORTBYNAME } from './modules/samtools_sortbyname'
include { SAMTOOLS_FIXMATE } from './modules/samtools_fixmate'
include { SAMTOOLS_SORTBYCOORD } from './modules/samtools_sortbycoord'
include { SAMTOOLS_MARKDUP } from './modules/samtools_markdup'
include { SAMTOOLS_SUBSETBAM } from './modules/samtools_subsetbam'
include { SAMTOOLS_INDEXSUBSET } from './modules/samtools_indexsubset'
include { SAMTOOLS_INDEX } from './modules/samtools_index'
// STRUCTURAL VARIANT DETECTION
include { MANTA_SVDETECTION } from './modules/manta_svdetection'
include { DELLY_SVDETECTION } from './modules/delly_svdetection'
include { PICARD_CREATESEQUENCEDICT } from './modules/picard_createseqeuncedict' 
include { SVABA_SVDETECTION } from './modules/svaba_svdetection'
include { SVABA_REWRITE_INSERTIONS } from './modules2/svaba_fixvcf_no_chmod'
include { SURVIVOR_SVINTEGRATION; SURVIVOR_INSINTEGRATION } from './modules2/survivor_svintegration_workdir'
include { TABIX_INDEXVCF; TABIX_INDEXDELLYVCF; TABIX_INDEXSVABAVCF} from './modules/tabix_indexvcf'
include { BCFTOOLS_VIEW } from './modules/bcftools_view'
include { BCFTOOLS_FILTERMANTA; BCFTOOLS_FILTERDELLY; BCFTOOLS_FILTERSVABA; BCFTOOLS_FILTERFORINS; BCFTOOLS_CONCAT; BCFTOOLS_UNZIP; BCFTOOLS_REMOVE_COMMONSV; BCFTOOLS_FILTER_CONTIGS; BCFTOOLS_FILTERSNV; BCFTOOLS_FILTERSIZE } from './modules2/bcftools_filter_WG'
include { FILTER_COMMON_SVS } from './modules2/scr_find_overlapping_SV_samplename'
include { BEDTOOLS_REMOVE_REPEATS } from './modules2/bedtools_remove_samplename'
include { ANNOTSV_SVANNOTATION; ANNOTSV_FORREPEATS } from './modules2/annotsv_svannotation_samplename'
// COPY NUMBER VARIANT ANALYSIS
include { GATK_PREPROCESS_INTERVALS; GATK_ANNOTATE_INTERVALS; GATK_COLLECT_READ_COUNTS; GATK_CREATE_PON } from './modules/gatk_createPoN'
include { PICARD_ADDREADGROUPS } from './modules2/picard_addreadgroups_1input/main.nf'
include { GATK_MODEL_SEGMENTS_NO_PON } from './modules2/gatk_cnvdetection_nopon'
include { GATK_COLLECT_READ_COUNTS_TUMOR; DENOISE_READ_COUNTS; COLLECT_ALLELIC_COUNTS_TUMOR; GATK_MODEL_SEGMENTS_TUMOR_ONLY; GATK_CALL_COPYRATIO_SEGMENTS; GATK_PLOT_DENOISED_COPYRATIOS; GATK_PLOT_MODELED_SEGMENTS } from './modules2/gatk_cnvdetection_samplename'
// SNV ANALYSIS
include { GATK_MUTECT2; GATK_FILTERMUTECTCALLS; GATK_VARIANTANNOTATION } from './modules2/gatk_mutect2_WG'
// include { GATK_MUTECT2_PURECN } from './modules2/gatk_mutect2_purecn' (PURE_CN ANALYSIS IS BEST WITH WES NOT WGS)
include { ANNOVAR } from './modules2/annovar_WG'
include { ANNOVAR_GNOMAD } from './modules2/annovar_gnomad_WG'


workflow {

    // ASSIGN CRAM FILES TO CHANNEL: sample_ch

    Channel
      .fromPath(params.samplesheet)
      .splitCsv(header: true)
      .map { row -> tuple(row.name, file(row.path)) }
      .set { sample_ch }

    // INDEX THE HUMAN GENOME FASTA FILE

    BWA_INDEX(params.ref)

    // PASS THE SAMPLES TO SAMTOOLS QUICK CHECK TO ASSESS FOR FILE CORRUPTION
    // THEN CONVERT THE CRAM FILES TO PAIRED END FASTQ FILES

    SAMTOOLS_QUICKCHECK(sample_ch)

    SAMTOOLS_CRAMTOFASTQ(sample_ch, params.ref)

    // REMOVE MOUSE CONTMAINATION:
    // CREATE XENGSORT INDEX FILES FOR BOTH HUMAN AND MOUSE GENOMES
    // PASS THE SAMPLES AND INDEX FILES TO XENGSORT CLASSIFY TO SORT THE MOUSE AND HUMAN READS

    XENGSORT_INDEX(params.mouse_ref, params.human_ref)

    XENGSORT_CLASSIFY(XENGSORT_INDEX.out.hash, XENGSORT_INDEX.out.info, SAMTOOLS_CRAMTOFASTQ.out.fastq_pair)

    // QUALITY CONTROL:
    // PASS THE RESULTING HUMAN CLASSIFIED READS TO SEQKIT STATS TO ASSESS READ QUALITY
    // BBMAP REPAIR RE PAIR UNPAIRED READS AND REMOVE SINGLE
  
    // TRIMMOMATIC PE CLEAN PAIRED END READS AND REMOVES ADAPTER SEQUENCES

    SEQKIT_STATS(XENGSORT_CLASSIFY.out.combined_pair)
    
    BBMAP_REPAIR(XENGSORT_CLASSIFY.out.combined_pair)

    TRIMMOMATIC_PE(BBMAP_REPAIR.out.bbmap_pair)

    // ALIGN THE PAIRED END FASTQ FILES TO HUMAN GENOME AND CREATE BAM FILE

    BWA_MEM(params.ref, BWA_INDEX.out.bwa_index, TRIMMOMATIC_PE.out.readp)

    // PREPROCESS BAM FILE:
    // FLAGSTAT ASSESSES THE ALIGNMENT QUALITY OF THE SAMPLES
    // SORTBYNAME SORTS THE READS BY NAME
    // FIXMATE FIXES AND MATES READS USING COORDINATES
    // SORTBYCOORD SORTS READS BY NUMERICAL COORDINATES
    // MARKDUP MARKS DUPLICATE READS
    // INDEX INDEXES THE BAM FILES

    bwa_bam_ch = BWA_MEM.out.bam
      .flatten()
      .filter { !it.name.contains('.tmp.') }

    SAMTOOLS_FLAGSTAT(bwa_bam_ch)

    SAMTOOLS_SORTBYNAME(bwa_bam_ch)

    SAMTOOLS_FIXMATE(SAMTOOLS_SORTBYNAME.out.sorted_bam)
    
    SAMTOOLS_SORTBYCOORD(SAMTOOLS_FIXMATE.out.fixmate_bam)

    SAMTOOLS_MARKDUP(SAMTOOLS_SORTBYCOORD.out.sorted_bam)

    SAMTOOLS_INDEX(SAMTOOLS_MARKDUP.out.markdup_bam)

    // ZOOM IN ON P53/TCAB1 LOCUS
    // CREATE SUBSETTED BAM FILES WITH SUBSETBAM
    // CREATE AN INDEX FOR THE SUBSETTED BAM FILE

    SAMTOOLS_SUBSETBAM(SAMTOOLS_INDEX.out.bam_index)

    SAMTOOLS_INDEXSUBSET(SAMTOOLS_SUBSETBAM.out.bam)

    // STRUCTURAL VARIANT ANALYSIS:
    // MANTA
    // DELLY
    // SVABA

    // // MANTA:
    // // DETECT STRUCTURAL VARIANT USING MANTA
    // // INDEX THE MANTA VCF FILE
    // // FILTER SV CALLS FOR >= 15 SUPPORTING READS

    MANTA_SVDETECTION(SAMTOOLS_INDEX.out.bam_index, params.ref, params.ref_index)

    // // HARD PATH TO FILES
    // // create a channel path to the manta vcf, will need to index it

    TABIX_INDEXVCF(MANTA_SVDETECTION.out.dipvcf)
    
    BCFTOOLS_FILTERMANTA(MANTA_SVDETECTION.out.dipvcf)

    // DELLY:
    // DETECT STRUCTURAL VARIANT USING DELLY
    // INDEX THE DELLY VCF FILE
    // FILTER SV CALLS FOR >= 15 SUPPORTING READS

    DELLY_SVDETECTION(SAMTOOLS_INDEX.out.bam_index, params.ref, params.ref_index)

    TABIX_INDEXDELLYVCF(DELLY_SVDETECTION.out.delly_vcf)

    BCFTOOLS_FILTERDELLY(DELLY_SVDETECTION.out.delly_vcf)

    // // SVABA:
    // // PICARD CREATES A DICTIONARY FOR THE HUMAN GENOME REFERENCE FILE
    // // DETECT STRUCTURAL VARIANT USING SVABA
    // // INDEX THE SVABA VCF FILE
    // // FILTER SV CALLS FOR >= 15 SUPPORTING READS

    PICARD_CREATESEQUENCEDICT(params.ref)

    SVABA_SVDETECTION(SAMTOOLS_INDEX.out.bam_index, params.ref, params.ref_index, BWA_INDEX.out.bwa_index, PICARD_CREATESEQUENCEDICT.out.dict)

    TABIX_INDEXSVABAVCF(SVABA_SVDETECTION.out.svaba_sv)

    BCFTOOLS_FILTERSVABA(SVABA_SVDETECTION.out.svaba_sv)

    // MERGE INSERTION CALLS IN SVABA SO IT IS COMPATIBLE IN SURVIVOR

    SVABA_REWRITE_INSERTIONS(BCFTOOLS_FILTERSVABA.out.filtered_vcf, params.fix_svaba_ins)

    // // MERGE CALLS USING 2 METHODS:
    // // CREATE A CHANNEL COMBINING CALLS FOR EVERY SAMPLE FROM THE 3 DETECTORS
    // // METHOD 1: PASS THE CHANNEL TO SURVIVOR TO MERGE THE CALLS IN EACH SAMPLE FOR TRA, DEL, DUP AND INV
    // // METHOD 2: PASS THE CHANNEL TO SURVIVOR TO MERGE ALL CALLS IN EVERY SAMPLE TO MAINTAIN INSERTIONS
    // ////// FILTER METHOD 2 FOR INS ONLY
    // // CONCATENATE BOTH METHODS INTO ONE MERGED VCF
    
    // Build sv_channel from BCFTOOLS filter outputs directly
    BCFTOOLS_FILTERDELLY.out.filtered_vcf
        .map { vcf -> 
            def sample = vcf.getName().replace('-combined_1P_target.delly.filtered.vcf', '')
            tuple(sample, vcf)
        }
        .join(
            BCFTOOLS_FILTERMANTA.out.filtered_vcf
            .map { vcf -> 
                def sample = vcf.getName().replace('-combined_1P_target.manta.filtered.vcf', '')
                tuple(sample, vcf)
            }
        )
        .join(
            BCFTOOLS_FILTERSVABA.out.filtered_vcf
            .map { vcf ->
                def sample = vcf.getName().replace('-combined_1P.svaba.filtered.vcf', '')
                tuple(sample, vcf)
            }
        )
        .map { sample, delly, manta, svaba -> tuple(sample, delly, manta, svaba) }
        .set { sv_channel }

    SURVIVOR_SVINTEGRATION(sv_channel, SVABA_REWRITE_INSERTIONS.out.svaba_ins_fixed)

    SURVIVOR_INSINTEGRATION(sv_channel, SVABA_REWRITE_INSERTIONS.out.svaba_ins_fixed)

    BCFTOOLS_FILTERFORINS(SURVIVOR_INSINTEGRATION.out.survivor_vcf)

    SURVIVOR_SVINTEGRATION.out.survivor_vcf
      .join(BCFTOOLS_FILTERFORINS.out.filtered_vcf)
      .map { sample, sv_vcf, ins_vcf -> tuple(sample, sv_vcf, ins_vcf) }
      .set { survivor_vcfs_ch }

    BCFTOOLS_CONCAT(survivor_vcfs_ch)

    // /// SUBSET THE SV CALLS FOR THE P53/TCAB1 REGION

    BCFTOOLS_VIEW(BCFTOOLS_CONCAT.out.vcf_allSV)

    // // // THE COMMON SV VCF FILE IS ZIPPED, UNZIP IT USING THE INDEX FILE

    // BCFTOOLS_UNZIP(params.common_var, params.common_var_index)

    // // Create an AnnotSV tsv for the samples
    // // The purpose of creating an annotsv tsv prior to final processing is to provide an intermediate step with a flexible data type

    ANNOTSV_FORREPEATS(BCFTOOLS_CONCAT.out.concat_vcf)

    // // create a bed file storing the coordinates of sample SV that overlap with common SVs

    FILTER_COMMON_SVS(ANNOTSV_FORREPEATS.out.annotsv, params.common_sv_bed, params.find_overlapping)

    BCFTOOLS_CONCAT.out.concat_vcf
      .join(FILTER_COMMON_SVS.out.overlap)
      .map { sample, vcf, bed -> tuple(sample, vcf, bed) }
      .set { common_sv_ch }

    // // REMOVE COMMON SVs IN THE GREATER POPULATION

    BCFTOOLS_REMOVE_COMMONSV(common_sv_ch)

    // // REMOVE SV CALLS FOUND WITHIN REPETITIVE REGIONS OF THE GENOME AS ANNOTATED BY REPEATMASKER

    BEDTOOLS_REMOVE_REPEATS(BCFTOOLS_REMOVE_COMMONSV.out.common_filtered_vcf, params.repeat_regions)

    // // FILTER SV CALLS STRICTLY FOR CANONICAL CHROMOSOMES

    BCFTOOLS_FILTER_CONTIGS(BEDTOOLS_REMOVE_REPEATS.out.repetitive_filtered_vcf)

    // // FILTER SV CALLS for 1Mb size limit

    BCFTOOLS_FILTERSIZE(BCFTOOLS_FILTER_CONTIGS.out.final_vcf)

    // // ANNOTATE THE FINAL SET OF STRUCTURAL VARIANTS USING ANNOTSV

    ANNOTSV_SVANNOTATION(BCFTOOLS_FILTERSIZE.out.filtered_vcf)

    // // CREATE A CHANNEL TO PAIR EACH NORMAL BAM WITH THE CORRESPONDING BAI
    
    Channel
    .fromPath(params.matched_normals, checkIfExists: true)
    .splitCsv(header: true, sep: ',', quote: '"', strip: true)
    .map { row ->
        def bam = file(row.Path)         // header is "Path"
        def bai = file(row.Path_Index)   // header is "Path_Index"
        if( !bam.exists() ) error "BAM not found: ${bam}"
        if( !bai.exists() ) error "BAI not found: ${bai}"
        tuple(bam, bai)
    }
    .set { normals_ch }

    // // CREATE A CHANNEL TO STORE CANONICAL CONTIGS 

    Channel
    .value(['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20', 'chr21','chr22'])
    .set { ch_primary_contigs }

    // // COPY NUMBER VARIANT ANALYSIS WITH GATK:
    // // BUILD A PANEL OF NORMALS
    // // PERFORM COPY NUMBER VARIANT ANALYSIS

    // //// BUILD A PANEL OF NORMALS:
    // //// PREPOCESS INTERVALS BUILD BINNED WINDOWS OF THE GENOME
    // //// COLLECT READ COUNTS FOR THE NORMAL SAMPLES
    // //// CREATE A CHANNEL TO COMBINE THE READ COUNTS INTO ONE OBJECT
    // //// CREATE PANEL OF NORMALS (PON) FROM THE COMBINED OBJECT

    GATK_PREPROCESS_INTERVALS(params.ref, params.dict, params.ref_index, ch_primary_contigs)

    GATK_COLLECT_READ_COUNTS(normals_ch, params.ref, params.ref_index, params.dict, GATK_PREPROCESS_INTERVALS.out.intervals)

    ch_pon_counts = GATK_COLLECT_READ_COUNTS.out.counts_hdf5
                      .flatten()
                      .collect()

    GATK_CREATE_PON(ch_pon_counts)

    //// PERFORM COPY NUMBER VARIANT ANALYSIS:
    //// CREATE A CHANNEL TO PAIR EACH TUMOR BAM FILE WITH THE CORRESPONDING BAI INDEX
    //// ADD READ GROUP SPECIFICATION INTO THE BAM HEADERS
    //// COLLECT READ COUNTS FOR THE TUMOR SAMPLES
    //// DENOISE THE TUMOR READ COUNTS USING THE PON
    //// COLLECT ALLELIC COUNTS IN TUMOR SAMPLES TO IMPROVE READ COUNT ESTIMATION
    //// CREATE A CHANEL FOR DENOISED AND ALLELIC COUNTS AND PAIR BY SAMPLE
    //// FEED DENOISED AND ALLELIC COUNTS TO MODEL SEGMENTS 
    //// MODEL SEGMENTS IS TUMOR ONLY GIVEN WE DO NOT HAVE MATCHED NORMAL SAMPLES
    //// CALL COPY RATIO SEGMENTS DETECTS DEVIATION AWAY FROM EXPECTED COPY NUMBER
    //// PLOT DENOISED COPY RATIO CREATES COPY RATIO FIGURES
    //// CREATE CHANNELS FOR THE MODEL AND THEN ANOTHER TO PAIR IT WITH DENOISED DATA
    //// PLOT MODELLED SEGMENTS ALLOWS FOR FUTHER ALLELIC INSIGHT
    
    SAMTOOLS_INDEX.out.bam_index
    .map { bam, bai ->
        def sample = bam.getName().replace('-combined_1P.markdup.bam', '')
        tuple(sample, bam, bai)
    }
    .set { indexed_bam_ch }

    PICARD_ADDREADGROUPS(indexed_bam_ch)

    GATK_COLLECT_READ_COUNTS_TUMOR(PICARD_ADDREADGROUPS.out.bam_rg, params.ref, params.ref_index, params.dict, GATK_PREPROCESS_INTERVALS.out.intervals)

    DENOISE_READ_COUNTS(GATK_COLLECT_READ_COUNTS_TUMOR.out.counts_hdf5, GATK_CREATE_PON.out.pon)

    COLLECT_ALLELIC_COUNTS_TUMOR(PICARD_ADDREADGROUPS.out.bam_rg, params.ref, params.ref_index, params.dict, params.snps, params.snps_index)

    DENOISE_READ_COUNTS.out.denoised
      .join(COLLECT_ALLELIC_COUNTS_TUMOR.out.allelic_counts)
      .map { sample, den, ac -> tuple(sample, den, ac) }
      .set { CH_DENOISED_ALLELIC }

    GATK_MODEL_SEGMENTS_TUMOR_ONLY(CH_DENOISED_ALLELIC)

    GATK_CALL_COPYRATIO_SEGMENTS(GATK_MODEL_SEGMENTS_TUMOR_ONLY.out.cr_seg)

    GATK_PLOT_DENOISED_COPYRATIOS(DENOISE_READ_COUNTS.out.standard, DENOISE_READ_COUNTS.out.denoised, params.dict)

    GATK_MODEL_SEGMENTS_TUMOR_ONLY.out.model_segments_out
      .join(DENOISE_READ_COUNTS.out.denoised)
      .map { sample, mod, den -> tuple(sample, den, mod) }
      .set { CH_DENOISED_MODEL }

    GATK_PLOT_MODELED_SEGMENTS(CH_DENOISED_MODEL, params.dict)

    // SNV Analysis with GATK Mutect2:
    // A tumor-only protocol in order to accomodate PDX samples without matched normals

    // using the samplesheet with all samples, we assign the sample name, bam and bam index file to a channel

    Channel
    .fromPath("${params.samplesheet_PDX}")
    .splitCsv(header: true)
    .map { row ->
        tuple(row.Final_Identifier, file(row.Path), file(row.Path_Index))
    }
    .set { PDX_samples_ch }

    // Call the channel, reference genome and indexes, and gnomad germline reference for Mutect2 SNV analysis using the tumor-only mode
    
    GATK_MUTECT2(PDX_samples_ch, params.ref, params.ref_index, params.dict, params.af_only_gnomad, params.af_only_gnomad_index, params.pon, params.pon_index)

    // GATK_MUTECT2_PURECN(PDX_samples_ch, params.ref, params.ref_index, params.dict, params.af_only_gnomad, params.af_only_gnomad_index, params.pon, params.pon_index)

    GATK_FILTERMUTECTCALLS(GATK_MUTECT2.out.mutect2_vcf, GATK_MUTECT2.out.mutect2_tbi, GATK_MUTECT2.out.mutect2_stats, params.ref, params.ref_index, params.dict)

    BCFTOOLS_FILTERSNV(GATK_FILTERMUTECTCALLS.out.filtered_mutect_vcf)

    // Load PDX BAM files from samplesheet, keyed by Final_Identifier

    Channel
      .fromPath("${params.samplesheet_PDX}")
      .splitCsv(header: true)
      .map { row ->
          def id  = (row.Final_Identifier as String).trim()
          def bam = file(row.Path)
          def bai = file(row.Path_Index)
          tuple(id, bam, bai)
      }
      .set { bams_ch }

    // Join filtered SNV VCFs with their corresponding BAM files by sample ID
    // to prepare inputs for variant annotation
    
    BCFTOOLS_FILTERSNV.out.filtered_snv
    .join(bams_ch)
    .map { id, vcf, bam, bai -> tuple(id, vcf, bam, bai) }
    .set { varannot_ch }

    GATK_VARIANTANNOTATION(varannot_ch, params.dbsnp, params.dbsnp_idx, params.ref, params.ref_index, params.dict)

    ANNOVAR(GATK_VARIANTANNOTATION.out.variantannotated_vcf, params.humandb_dir, params.table_annovar, params.convert2annovar, params.annotate_variant, params.coding_change, params.retrieve_seq_from_fasta, params.variants_reduction)

    ANNOVAR_GNOMAD(GATK_VARIANTANNOTATION.out.variantannotated_vcf, params.humandb_dir, params.table_annovar, params.convert2annovar, params.annotate_variant, params.coding_change, params.retrieve_seq_from_fasta, params.variants_reduction)


}