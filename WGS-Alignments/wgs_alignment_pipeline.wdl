version 1.0

workflow WGSpipelines {
    input {
        String halb
        Array[String] extra_halbs = []
        Array[File] extra_fastqs = []
        String? title
        String? description
        Boolean use_atlas = true
        File reference_fa
        File reference_fa_alt = "~{reference_fa}.alt"
        File reference_fa_amb = "~{reference_fa}.amb"
        File reference_fa_ann = "~{reference_fa}.ann"
        File reference_fa_bwt = "~{reference_fa}.bwt"
        File reference_fa_fai = "~{reference_fa}.fai"
        File reference_fa_pac = "~{reference_fa}.pac"
        File reference_fa_sa = "~{reference_fa}.sa"
        String license_server
        File dbsnp = "/cluster/home/jtaylor/reference_files/cromwell_data_files/dbsnp_146.hg38.vcf.gz"
        File dbsnp_tbi = "/cluster/home/jtaylor/reference_files/cromwell_data_files/dbsnp_146.hg38.vcf.gz.tbi"
        File mills = "/cluster/home/jtaylor/reference_files/cromwell_data_files/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
        File mills_tbi = "/cluster/home/jtaylor/reference_files/cromwell_data_files/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"
        File strelka_call_contigs = "/cluster/home/jtaylor/reference_files/cromwell_data_files/hg38_contigs.bed.gz"
        File strelka_call_contigs_tbi = "/cluster/home/jtaylor/reference_files/cromwell_data_files/hg38_contigs.bed.gz.tbi"

        Boolean write_recal = false
        Boolean run_strelka = false
        Boolean run_gatk = false
        Boolean write_cram = false
		
        Boolean remove_dups = true
        String read_filter = "-q "
        String genotype_model = "multinomial"

        String read_1_identifier = "_1_"
        String read_2_identifier = "_2_"
    }
    call LocateFastqs as lf {
        input:
            halb=halb,
            use_atlas=use_atlas,
            extra_halbs=extra_halbs,
            extra_fastqs=extra_fastqs,
            read_1_identifier=read_1_identifier,
            read_2_identifier=read_2_identifier
    }
    scatter (fqPair in lf.fastq_pairs) {
        call SentieonPostalt as alignment {
            input:
                halb=halb,
                fqPair=fqPair,
                reference_fa=reference_fa,
                reference_fa_alt=reference_fa_alt,
                reference_fa_amb=reference_fa_amb,
                reference_fa_ann=reference_fa_ann,
                reference_fa_bwt=reference_fa_bwt,
                reference_fa_fai=reference_fa_fai,
                reference_fa_pac=reference_fa_pac,
                reference_fa_sa=reference_fa_sa,
                license_server=license_server
        }
        call Fastp as fastp {
            input:
                halb=halb,
                fqPair=fqPair
        }
    }
    call SentieonDedup as dedup {
        input:
            halb=halb,
            bams=alignment.bam,
            bais=alignment.bai,
            license_server=license_server,
            remove_dups=remove_dups
    }
    call SentieonRecal as recal {
        input:
            halb=halb,
            dedup_bam=dedup.bam,
            dedup_bai=dedup.bai,
            dbsnp=dbsnp,
            dbsnp_tbi=dbsnp_tbi,
            mills=mills,
            mills_tbi=mills_tbi,
            reference_fa=reference_fa,
            reference_fa_fai=reference_fa_fai,
            license_server=license_server
    }
    call WgsStats as stats {
        input:
           halb=halb,
           reference_fa=reference_fa,
           reference_fa_fai=reference_fa_fai,
           dedup_bam=dedup.bam,
           dedup_bai=dedup.bai,
           license_server=license_server,
           read_filter=read_filter,
           recal_table=recal.recal_table
    }
    call CompileFastp as compile_fastp {
        input:
            halb=halb,
            fastp_jsons=fastp.fastpJson
    }
    if (run_strelka || write_recal) {
        call WriteRecal as save_recal {
            input:
                halb=halb,
                dedup_bam=dedup.bam,
                dedup_bai=dedup.bai,
                recal_table=recal.recal_table,
                read_filter=read_filter,
                reference_fa=reference_fa,
                reference_fa_fai=reference_fa_fai,
                license_server=license_server
        }
        if (run_strelka) {
            call Strelka2 as strelka2 {
                input:
                    bam=save_recal.bam,
                    bai=save_recal.bai,
                    reference_fa=reference_fa,
                    reference_fa_fai=reference_fa_fai,
                    call_contigs=strelka_call_contigs,
                    call_contigs_tbi=strelka_call_contigs_tbi
            }
            call AnnotateDBSNP as annotate {
                input:
                    halb=halb,
                    vcf=strelka2.vcf,
                    vcf_tbi=strelka2.vcf_tbi,
                    dbsnp=dbsnp,
                    dbsnp_tbi=dbsnp_tbi
            }
        }
    }
    if (run_gatk) {
        call SentieonGATK4 as gatk4 {
            input:
                halb=halb,
                bam=dedup.bam,
                bai=dedup.bai,
                recal_table=recal.recal_table,
                read_filter=read_filter,
                reference_fa=reference_fa,
                reference_fa_fai=reference_fa_fai,
                license_server=license_server,
                genotype_model=genotype_model,
                dbsnp=dbsnp,
                dbsnp_tbi=dbsnp_tbi
        }
    }
    if(write_cram) {
        call BamToCram as bam2cram {
            input:
                halb=halb,
                dedup_bam=dedup.bam,
                dedup_bai=dedup.bai,
                recal_table=recal.recal_table,
                read_filter=read_filter,
                reference_fa=reference_fa,
                reference_fa_fai=reference_fa_fai,
                license_server=license_server
        }
    }
    output {
        File? bam=save_recal.bam
        File? bai=save_recal.bai
        File? vcf=annotate.annot_vcf
        File? vcf_tbi=annotate.annot_tbi
        File? strelka_gvcf=strelka2.gvcf
        File? strelka_gvcf_tbi=strelka2.gvcf_tbi
        File? gatk4_gvcf=gatk4.gvcf
        File? gatk4_gvcf_tbi=gatk4.gvcf_tbi
        File? gatk4_vcf=gatk4.vcf
        File? gatk4_vcf_tbi=gatk4.vcf_tbi
        File? cram=bam2cram.cram
        File? crai=bam2cram.crai
        File dedup_metrics=dedup.metrics
        File dedup_readGroup=dedup.readGroup
        File aligment_stat = stats.alignment_stat
        File insert_stat = stats.insert_stat
        File quality_yield_stat = stats.quality_yield_stat
        File quality_distro_stat = stats.quality_distro_stat
        File gc_bias_stat = stats.gc_bias_stat
        File gc_bias_all = stats.gc_bias_all
        File wgs_stat = stats.wgs_stat
        File fastpJson = compile_fastp.fastpJson
    }
}

task LocateFastqs {
    meta {
        volatile: true
    }
    input {
        String halb
        Boolean use_atlas
        Array[String] extra_halbs
        Array[File] extra_fastqs
        String read_1_identifier
        String read_2_identifier
    }
    command {
        if [ ~{use_atlas} == "true" ]; then
            atlas fmf init
            atlas fmf add ~{halb}
            for L in $(cat ~{write_lines(extra_halbs)}); do
                atlas fmf add $L
            done
            atlas fmf download
        fi

        mkdir HALB_extra_fastqs
        for F in $(cat ~{write_lines(extra_fastqs)}); do
            ln -s $F HALB_extra_fastqs/
        done
    }
    output {
        Array[Pair[File, File]] fastq_pairs = zip(glob("HALB*/*~{read_1_identifier}*fastq.gz"), glob("HALB*/*~{read_2_identifier}*fastq.gz"))
    }
    runtime {
        cpu: 1
        memory: '1 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: "" 
    }
}

task SentieonPostalt {
    input {
        String halb
        Pair[File, File] fqPair
        File reference_fa
        File reference_fa_alt
        File reference_fa_amb
        File reference_fa_ann
        File reference_fa_bwt
        File reference_fa_fai
        File reference_fa_pac
        File reference_fa_sa
        String license_server
        String tempParams = ""
        Int threads = 16
    }
    String bamFN = "~{halb}.aligned.bam"
    String baiFN = "~{halb}.aligned.bam.bai"
    command {
        RG_ID=$(basename ~{fqPair.left} | cut -f1-2 -d_)
        RG_OPTIONS="@RG\tID:$RG_ID\tSM:~{halb}\tLB:~{halb}\tPL:ILLUMINA"
        export SENTIEON_LICENSE="${license_server}"

        sentieon \
            bwa mem -M \
            -R "$RG_OPTIONS" \
            -t ${threads} \
            -K 10000000 \
            ${reference_fa} \
            ${fqPair.left} ${fqPair.right} | \
        k8 \
            /opt/bwakit/bwa-postalt.js \
            ${reference_fa_alt} | \
        sentieon util sort ${tempParams} \
            --bam_compression 1 \
            -r ${reference_fa} \
            -o ${bamFN} \
            -t ${threads} \
            --sam2bam \
            -i - 
    }
    output {
        File bam = "${bamFN}"
        File bai = "${baiFN}"
    }
    runtime { 
        cpu: threads
        memory: '24 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}
task Fastp {
    input {
        String halb
        Pair[File, File] fqPair
    }
    String htmlFN = "${basename(fqPair.left, '.fast.gz')}.fastp.html"
    String jsonFN = "${basename(fqPair.left, '.fast.gz')}.fastp.json"
    command {
        fastp -i ${fqPair.left} -I ${fqPair.right} -h ${htmlFN} -j ${jsonFN}
    }
    output {
        File fastpJson = "${jsonFN}"
    }
    runtime {
        cpu: 8
        memory: '24 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}

task SentieonDedup {
    input {
        String halb
        Array[File] bams
        Array[File] bais
        String license_server
        String tempParams = ""
        Boolean remove_dups
        Int threads = 16
    }
    String bamFN = "~{halb}.dedup.bam"
    String baiFN = "~{halb}.dedup.bam.bai"
    String scoreFN = "~{halb}.score.txt"
    String metricsFN = "~{halb}.dedupmetrics.txt"
    String readGroupFN = "~{halb}.readGroup.txt"
    String dupOption = if remove_dups then "--rmdup" else ""
    command {
        export SENTIEON_LICENSE="${license_server}"
        sentieon driver \
            -t ${threads} \
            -i ${sep=" -i " bams} \
            --algo LocusCollector \
            --fun score_info \
            ${scoreFN} && \
        sentieon driver ${tempParams} \
            -t ${threads} \
            -i ${sep=" -i " bams} \
            --algo Dedup \
            ${dupOption} \
            --score_info ${scoreFN} \
            --metrics ${metricsFN} \
            --bam_compression 1 \
            ${bamFN} && \
        echo "${length(bams)}" > ${readGroupFN}
    }
    output {
        File bam = "${bamFN}"
        File bai = "${baiFN}"
        File metrics = "~{metricsFN}"
        File readGroup = "~{readGroupFN}"
    }
    runtime {
        cpu: threads
        memory: '128 GB'
        disks: "local-disk 500 SSD"
        container: ""
        docker: ""
    }
}

task SentieonRecal {
    input {
        String halb
        File dedup_bam
        File dedup_bai
        File dbsnp
        File dbsnp_tbi
        File mills
        File mills_tbi
        File reference_fa
        File reference_fa_fai
        String license_server
        String tempParams = ""
        Int threads = 16
    }
    String recal_table_fn = "~{halb}.recal.table"
    String recal_bam = "~{halb}.recal.bam"
    String recal_bai = "~{halb}.recal.bam.bai"
    command {
        export SENTIEON_LICENSE="${license_server}"

        sentieon driver \
            -r ${reference_fa} \
            -t ${threads} \
            -i ${dedup_bam} \
            --algo QualCal \
            -k ${dbsnp} \
            -k ${mills} \
            ${recal_table_fn}
    }
    output {
        File recal_table = "~{recal_table_fn}"
    }
    runtime {
        cpu: threads
        memory: '64 GB'
        disks: "local-disk 500 SSD"
        container: ""
        docker: ""
    }
}

task CompileFastp {
    input {
        String halb
        Array[File] fastp_jsons
    }
    String halbJsonFN = "~{halb}.fastp.json"
    command {
        python /usr/local/bin/compile_fastp.py ${sep=" " fastp_jsons} > ${halbJsonFN}
    }
    output {
        File fastpJson = "${halbJsonFN}"
    }
    runtime {
        cpu: 1
        memory: '8 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}

task WgsStats {
    input {
       String halb
       String license_server
       File reference_fa
       File reference_fa_fai
       File dedup_bam
       File dedup_bai
       File recal_table
       String read_filter
       Int threads = 16
    }

    String alignment_stat_fn = "~{halb}.alignment_metrics.txt"
    String insert_stat_fn = "~{halb}.insert_metrics.txt"
    String quality_yield_stat_fn = "~{halb}.quality_yield_metrics.txt"
    String quality_distro_stat_fn = "~{halb}.quality_distro_metrics.txt"
    String gc_bias_stat_fn = "~{halb}.gcbias_metrics.txt"
    String gc_bias_all_fn = "~{halb}.gcbias_all.txt"
    String wgs_stat_fn = "~{halb}.output_metrics.txt"
    command {
        export SENTIEON_LICENSE="${license_server}"

        sentieon driver -t ~{threads} -r "~{reference_fa}" ~{read_filter}~{recal_table} -i ~{dedup_bam} \
          --algo AlignmentStat ~{alignment_stat_fn} \
          --algo InsertSizeMetricAlgo ~{insert_stat_fn} \
          --algo QualityYield ~{quality_yield_stat_fn} \
          --algo QualDistribution ~{quality_distro_stat_fn} \
          --algo GCBias --summary ~{gc_bias_stat_fn} ~{gc_bias_all_fn} \
          --algo WgsMetricsAlgo --base_qual_histogram true ~{wgs_stat_fn}
    }

    output {
        File alignment_stat ="~{alignment_stat_fn}"
        File insert_stat = "~{insert_stat_fn}"
        File quality_yield_stat = "~{quality_yield_stat_fn}"
        File quality_distro_stat = "~{quality_distro_stat_fn}"
        File gc_bias_stat = "~{gc_bias_stat_fn}"
        File gc_bias_all = "~{gc_bias_all_fn}"
        File wgs_stat = "~{wgs_stat_fn}"
    }
    runtime {
        cpu: threads
        memory: '64 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}

task WriteRecal {
    input {
        String halb
        File dedup_bam
        File dedup_bai
        File reference_fa
        File reference_fa_fai
        String license_server
        File recal_table
        String read_filter
        String tempParams = ""
        Int threads = 16
    }
    String recal_bam = "~{halb}.recal.bam"
    String recal_bai = "~{halb}.recal.bam.bai"
    command {
        export SENTIEON_LICENSE="~{license_server}"

        sentieon driver ~{tempParams} \
            -r ~{reference_fa} \
            -t ~{threads} \
            -i ~{dedup_bam} \
            ~{read_filter}~{recal_table} \
            --algo ReadWriter \
            ~{recal_bam}
    }
    output {
        File bam = "~{recal_bam}"
        File bai = "~{recal_bai}"
    }
    runtime {
        cpu: threads
        memory: '16 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}

task Strelka2 {
    input {
        File bam
        File bai
        File call_contigs
        File call_contigs_tbi
        File reference_fa
        File reference_fa_fai
        Int threads = 16
    }
    String run_dir="strelka2_rundir"
    String mem_GB="16"
    command {
        /usr/local/bin/configureStrelkaGermlineWorkflow.py \
            --bam ${bam} \
            --referenceFasta ${reference_fa} \
            --callRegions ${call_contigs} \
            --runDir ${run_dir} && \
        ${run_dir}/runWorkflow.py \
            -m local \
            -j ${threads} \
            -g ${mem_GB}
    }
    output {
        File vcf="${run_dir}/results/variants/variants.vcf.gz"
        File vcf_tbi="${run_dir}/results/variants/variants.vcf.gz.tbi"
        File gvcf="${run_dir}/results/variants/genome.vcf.gz"
        File gvcf_tbi="${run_dir}/results/variants/genome.vcf.gz.tbi"
    }
    runtime {
        cpu: threads
        memory: '16 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}

task AnnotateDBSNP {
    input {
        String halb
        File vcf
        File vcf_tbi
        File dbsnp
        File dbsnp_tbi
    }
    String out_fn="~{halb}_annotated.vcf.gz"
    command {
        bcftools annotate \
            -a ${dbsnp} \
            -c ID \
            -O z \
            -o ${out_fn} \
            ${vcf} && \
        tabix ${out_fn}
    }
    output {
        File annot_vcf="${out_fn}"
        File annot_tbi="${out_fn}.tbi"
    }
    runtime {
        cpu: 1
        memory: '4 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}

task SentieonGATK4 {
    input {
        String halb
        File bam
        File bai
        File reference_fa
        File reference_fa_fai
        File recal_table
        String read_filter
        String license_server
        String genotype_model
        File dbsnp
        File dbsnp_tbi
        Int threads = 16
    }
    String gvcf_fn="~{halb}.gvcf.gz"
    String gtbi_fn="~{halb}.gvcf.gz.tbi"
    String vcf_fn="~{halb}.vcf.gz"
    String tbi_fn="~{halb}.vcf.gz.tbi"
    command {
        export SENTIEON_LICENSE="${license_server}"

        sentieon driver \
            -r ${reference_fa} \
            -t ${threads} \
            -i ${bam} \
            ~{read_filter}~{recal_table} \
            --algo Haplotyper \
            -d ${dbsnp} \
            --genotype_model ~{genotype_model} \
            --emit_mode gvcf \
            ${gvcf_fn}

        sentieon driver \
        -r ${reference_fa} \
        -t ${threads} \
        -i ${bam} \
        ~{read_filter}~{recal_table} \
        --algo Haplotyper \
        --genotype_model ${genotype_model} \
        -d ${dbsnp} \
        ${vcf_fn}
    }
    output {
        File gvcf="${gvcf_fn}"
        File gvcf_tbi="${gtbi_fn}"
        File vcf="${vcf_fn}"
        File vcf_tbi="${tbi_fn}"
    }
    runtime {
        cpu: threads
        memory: '64 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}

task BamToCram {
    input {
        String halb
        File dedup_bam
        File dedup_bai
        File reference_fa
        File reference_fa_fai
        String license_server
        File recal_table
        String read_filter
        String tempParams = ""
        Int threads = 16
    }
    String recal_cram = "~{halb}.recal.cram"
    String recal_crai = "~{halb}.recal.cram.crai"
    command {
        export SENTIEON_LICENSE="~{license_server}"

        sentieon driver ~{tempParams} \
            -r ~{reference_fa} \
            -t ~{threads} \
            -i ~{dedup_bam} \
            ~{read_filter}~{recal_table} \
            --algo ReadWriter \
            ~{recal_cram}
    }
    output {
        File cram = "~{recal_cram}"
        File crai = "~{recal_crai}"
    }
    runtime {
        cpu: threads
        memory: '64 GB'
        disks: "local-disk 300 SSD"
        container: ""
        docker: ""
    }
}
