import glob

#sample_pattern = "data/AM1686*.bam"
sample_pattern = "data/*.bam"
sample_files = glob.glob(sample_pattern)
sample_ids = [sample[5:-4] for sample in sample_files]

#print(sample_ids)

sampleLane_pattern = "data/trimmed_Lan/*_R1_001.fastq"
sampleLane_files = glob.glob(sampleLane_pattern)
sampleLane_ids = [sampleLane[17:-13] for sampleLane in sampleLane_files]


def get_reps(wildcards):
    #bamRep_pattern = data/mapped_RG/{wildcards.sample}_*.bam
    #bamRep_files = glob.glob(bamRep_pattern)
    bamRep_files=glob.glob('data/mapped_RG/' + wildcards.sample + '_*.bam')
    return bamRep_files[0]
    #return bamRep_files[0]


## I am running the alt-aware alignment tutorial from: https://gatkforums.broadinstitute.org/gatk/discussion/8017/how-to-map-reads-to-a-reference-with-alternate-contigs-like-grch38
## I will add notes for roles required changes to match this tutorial

rule all:
    input:
        #expand("data/shuf/{sample}_shuf.bam", sample=sample_ids),
        expand("data/fastq/{sample}_R_001.pe.fq", sample=sample_ids),
        expand("qc/fastqc/fastq/{sample}_R_001.pe_fastqc.{ext}", sample=sample_ids, ext=['html', 'zip']),
        "qc/multiqc/fastq/multiqc_report.html",
        expand("data/fastqPE/{sample}_R1_001.fq", sample=sample_ids), expand("data/fastqPE/{sample}_R2_001.fq", sample=sample_ids),
        expand('data/trimmed/{sample}_{R}_001.fastq.gz', sample=sample_ids, R=['R1', 'R2']),
        expand("data/trimmed_Lan/{sample}.laneReport", sample=sample_ids),
        expand("data/uBAM/{sampleLane}.bam", sampleLane=sampleLane_ids),
        "refGenome/hs38DH.fa",
        expand("refGenome/hs38DH.fa.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        "refGenome/gatkIndex/genome.fa", "refGenome/gatkIndex/genome.fa.fai", "refGenome/gatkIndex/genome.dict",
        expand("data/mapped_reads/{sampleLane}.aln.bam", sampleLane=sampleLane_ids), 
        expand("data/mapped_reads/{sampleLane}.doneReport", sampleLane=sampleLane_ids),
##        expand("data/mapped_reads/{sample}.aln.bam", sample=sample_ids),
        expand("data/mapped_RG/{sampleLane}.bam", sampleLane=sampleLane_ids),
        expand("data/dedup/{sample}.bam", sample=sample_ids),
        expand("data/clean/{sample}.bam", sample=sample_ids),
        "knowVar/Homo_sapiens_assembly38.known_indels.vcf.gz","knowVar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        expand("data/recalib/{sample}.txt", sample=sample_ids),
        expand("data/recalib/{sample}.bam", sample=sample_ids),
        expand("vc/hapCaller_single/{sample}.g.vcf", sample=sample_ids),
        "vc/combined_patch.g.vcf",
        "vc/hapCaller_noped_raw.vcf", "vc/hapCaller_ped_raw.vcf",
        expand("vc/hapCaller_{fam}_ExcessHet.vcf", fam=['ped','noped']),
        "vc/hapCaller_ped_sitesOnly.vcf",
        "vc/ADHD_ped.indels.recal", "vc/ADHD_ped.indels.tranches",
        "vc/ADHD_ped.snps.recal", "vc/ADHD_ped.snps.tranches", 
        "vc/hapCaller_ped_ExcessHet_indelRecal.vcf",
        "vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf",
        "vc/ADHD_ped.variant_calling_detail_metrics", "vc/ADHD_ped.variant_calling_summary_metrics",
        expand("vc/hapCaller_ped_ExcessHet_Recal_{var}.vcf", var=['SNP','INDEL']),

rule retrieve_fastq:
    input:
        "data/{sample}.bam"
    output:
        output_pe="data/fastq/{sample}_R_001.pe.fq",
        #output_se="data/fastq/{sample}_R_001.se.fq"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    threads:1
    shell:
        '''
        mkdir -p data/tmp
        module load SAMtools/1.9
        samtools bamshuf -uOn 64 {input} data/tmp/tmp_{wildcards.sample} > data/shuf/{wildcards.sample}_shuf.bam
        samtools bam2fq data/shuf/{wildcards.sample}_shuf.bam -s data/fastq/{wildcards.sample}_R_001.se.fq > {output.output_pe}
        '''

rule fastqc_pre:
    input:
        "data/fastq/{sample}_R_001.pe.fq"
    output:
        html="qc/fastqc/fastq/{sample}_R_001.pe_fastqc.html",
        zip="qc/fastqc/fastq/{sample}_R_001.pe_fastqc.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    threads:1
    shell:
        '''
        #module load FastQC/0.11.5
        module load FastQC/0.11.7-Java-1.8.0_162
        fastqc -t {threads} -f fastq -noextract -o qc/fastqc/fastq {input}
        '''

rule multiQC_pre:
    input:
        html=expand("qc/fastqc/fastq/{sample}_R_001.pe_fastqc.html", sample=sample_ids)
    output:
        "qc/multiqc/fastq/multiqc_report.html",
        "qc/multiqc/fastq/multiqc_data.zip"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 32
    shell:
        '''
        source activate multiQC
        multiqc -z -o qc/multiqc/fastq qc/fastqc/fastq
        '''

rule split_fastq:
    input:
        "data/fastq/{sample}_R_001.pe.fq"
    output:
        R1="data/fastqPE/{sample}_R1_001.fq",
        R2="data/fastqPE/{sample}_R2_001.fq"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    shell:
        '''
        paste - - - - - - - - < {input} \
            | tee >(cut -f 1-4 | tr "\t" "\n" > {output.R1}) \
            | cut -f 5-8 | tr "\t" "\n" >  {output.R2}
        '''

rule trimmomatic_pe:
    input:
        r1="data/fastqPE/{sample}_R1_001.fq",
        r2="data/fastqPE/{sample}_R2_001.fq",
    output:
        r1="data/trimmed/{sample}_R1_001.fastq.gz",
        r2="data/trimmed/{sample}_R2_001.fastq.gz",
        sr1="data/trimmed/{sample}_R1.unpaired_001.fastq.gz",
        sr2="data/trimmed/{sample}_R2.unpaired_001.fastq.gz",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 4
    threads:1
    shell:
        '''
        module load Trimmomatic/0.38-Java-1.8.0_162
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.38.jar PE -threads {threads} -phred33 {input.r1} {input.r2} {output.r1} {output.sr1} {output.r2} {output.sr2} ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/TruSeq3-PE.fa:2:30:10:1 SLIDINGWINDOW:4:2 MINLEN:20
        '''

rule split_lanes:
    input:
        r1="data/trimmed/{sample}_R1_001.fastq.gz",
        r2="data/trimmed/{sample}_R2_001.fastq.gz",
    output:
        "data/trimmed_Lan/{sample}.laneReport",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    threads:1
    shell:
        '''
        zcat {input.r1} | awk -v id={wildcards.sample} 'BEGIN {{FS = ":"}} {{lane=$2 ; print > "data/trimmed_Lan/"id"_"lane"_R1_001.fastq" ; for (i = 1; i <= 3; i++) {{getline ; print > "data/trimmed_Lan/"id"_"lane"_R1_001.fastq"}}}}'
        zcat {input.r2} | awk -v id={wildcards.sample} 'BEGIN {{FS = ":"}} {{lane=$2 ; print > "data/trimmed_Lan/"id"_"lane"_R2_001.fastq" ; for (i = 1; i <= 3; i++) {{getline ; print > "data/trimmed_Lan/"id"_"lane"_R2_001.fastq"}}}}'
        ls data/trimmed_Lan/{wildcards.sample}_* > {output}
        '''

rule FastqToSam:
    input:
        r1="data/trimmed_Lan/{sampleLane}_R1_001.fastq",
        r2="data/trimmed_Lan/{sampleLane}_R2_001.fastq"
    output:
        "data/uBAM/{sampleLane}.bam"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        name=$(basename {input.r1})
        SM=$(echo $name | cut -d "_" -f1)
        LB=$SM                                                    # I assume each sample has one library prep
        PL="Illumina"
        RGID=$(echo $name | cut -d "_" -f1,2)

        #module load Java/jdk1.8.0
        module load Java/1.8.0_192
        source activate gatk
        #java -jar picard.jar FastqToSam
        gatk --java-options "-Xmx6G" FastqToSam \
        -F1={input.r1} \
        -F2={input.r2} \
        -O={output} \
        -SM=$SM \
        -LB=$LB \
        -PL=$PL \
        -RG=$RGID \
        --TMP_DIR="tmp/{wildcards.sampleLane}"
        '''

## In this pipeline, I will skip "MarkIlluminaAdapters" because I trimmed the fastq files with trimmomatic already


## Using "bwa.kit" to download and index the ref to follow the alt-aware alignment tutorial
rule download_ref:
    output:
        "refGenome/hs38DH.fa",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    shell:
        '''
        mkdir -p refGenome && cd refGenome
        ../bwa.kit/run-gen-ref hs38DH
        '''

rule bwa_index:
    input:
        "refGenome/hs38DH.fa"
    output:
        expand("refGenome/hs38DH.fa.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
    shell:
        '''
        bwa.kit/bwa index refGenome/hs38DH.fa
        '''

## I do separet bwa mapping to use "bwa.kit" to follow the alt-aware alignment tutorial. This mapping is done from fastq instead of reversing uBAM
rule bwa_map:
    input:
        expand("refGenome/hs38DH.fa.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
        r1="data/trimmed_Lan/{sampleLane}_R1_001.fastq",
        r2="data/trimmed_Lan/{sampleLane}_R2_001.fastq"
    output:
        bam="data/mapped_reads/{sampleLane}.aln.bam",
        log="data/mapped_reads/{sampleLane}.doneReport"  
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
    threads:4
    shell:
        '''
        bwa.kit/run-bwamem -t {threads} -o data/mapped_reads/{wildcards.sampleLane} refGenome/hs38DH.fa {input.r1} {input.r2} | sh
        echo {output.bam} > {output.log}
        '''

## clean_failed_map.sh: I used this script to clean failing bwa_map jobs before adding the idea of log files 

#rule bwa_map:
#    input:
#        expand("refGenome/hs38DH.fa.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
#        r1="data/trimmed/{sample}_R1_001.fastq.gz",
#        r2="data/trimmed/{sample}_R2_001.fastq.gz",
#    output:
#        "data/mapped_reads/{sample}.aln.bam",   
#    resources:
#        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
#        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
#    threads:4
#    shell:
#        '''
#        name=$(basename {input.r1})
#        SM=$(echo $name | cut -d "_" -f1)
#        LB=$SM                                                                # I am assuming each sample has one prep
#        PL="Illumina"
#        RGID=$(head -n1 <(zcat {input.r1}) | sed 's/:/_/g' |cut -d "_" -f1,2) # FASTQ from multiple lanes are not correct
#        PU=$RGID.$LB  
#        echo RGID $RGID LB $LB PL $PL PU $PU SM $SM >> test.txt
#
#        bwa.kit/run-bwamem -t {threads} -o data/mapped_reads/{wildcards.sample} -R "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU" refGenome/hs38DH.fa {input.r1} {input.r2} | sh
#        ls {output}
#        '''

rule GATK_index:
    input:
        "refGenome/hs38DH.fa"
    output:
        ref="refGenome/gatkIndex/genome.fa",
        index="refGenome/gatkIndex/genome.fa.fai",
        dict="refGenome/gatkIndex/genome.dict",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    shell:
        '''
        if [ ! -f refGenome/gatkIndex/genome.fa ];then ln -s ../hs38DH.fa refGenome/gatkIndex/genome.fa;fi
        #module load SAMTools/1.5
        #module load picardTools/1.89
        #samtools faidx "refGenome/gatkIndex/genome.fa"
        #java -Xmx4g -jar $PICARD/CreateSequenceDictionary.jar R= {input} O= {output.dict}
        module load SAMtools/1.9
        samtools faidx "refGenome/gatkIndex/genome.fa"
        module load Java/1.8.0_192
        source activate gatk
        gatk --java-options "-Xmx6G" CreateSequenceDictionary -R= {input} -O= {output.dict}
        '''


## To use the new MarkDuplicates technique, I am changing the SORT_ORDER (and thus removing CREATE_INDEX) unlike the pipeline in dogTumors.
## I am adding "--ATTRIBUTES_TO_RETAIN=XA" to follow the command in the alt-aware alignment tutorial
## Also, I am using BAM file generated already (by the bwa.kit) instead of reversing uBAM and mapping on the fly 
rule mergeBamAlignment:
    input:
        gatk_ref="refGenome/gatkIndex/genome.fa",
        ubam="data/uBAM/{sampleLane}.bam",
        bam="data/mapped_reads/{sampleLane}.aln.bam"
    output:
        bam="data/mapped_RG/{sampleLane}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        #module load Java/jdk1.8.0
        module load Java/1.8.0_192
        source activate gatk
        gatk --java-options "-Xmx6G" MergeBamAlignment \
        --ALIGNED_BAM={input.bam} \
        --UNMAPPED_BAM={input.ubam} \
        --OUTPUT={output.bam} \
        -R={input.gatk_ref} --SORT_ORDER=unsorted --ADD_MATE_CIGAR=true \
        --CLIP_ADAPTERS=false --CLIP_OVERLAPPING_READS=true \
        --INCLUDE_SECONDARY_ALIGNMENTS=true --MAX_INSERTIONS_OR_DELETIONS=-1 \
        --PRIMARY_ALIGNMENT_STRATEGY=MostDistant --ATTRIBUTES_TO_RETAIN=XS --ATTRIBUTES_TO_RETAIN=XA \
        --TMP_DIR="tmp2/{wildcards.sampleLane}"
        '''

# the new MarkDuplicates technique is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
# This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
# While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
rule MarkDuplicates:
    input:
        get_reps,
    output:
        bam="data/dedup/{sample}.bam",
        metrics="data/dedup/{sample}.metrics.txt",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 4,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        #inputs_bam = [$"--INPUT="bamRep for bamRep in {input}]
        while read bam;do
            input_bam=$"--INPUT="$bam;
            inputs_bam+=($input_bam);
        done < <(find data/mapped_RG -name "{wildcards.sample}_*.bam")

        #module load Java/jdk1.8.0
        module load Java/1.8.0_192
        source activate gatk
        gatk --java-options "-Xmx6G" MarkDuplicates \
        ${{inputs_bam[*]}} \
        --OUTPUT={output.bam} \
        --METRICS_FILE={output.metrics} \
        --VALIDATION_STRINGENCY SILENT \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        --ASSUME_SORT_ORDER "queryname" \
        --TMP_DIR="tmp3/{wildcards.sample}"
        '''

rule sort_fix_index:
    input:
        bam="data/dedup/{sample}.bam",
        gatk_ref="refGenome/gatkIndex/genome.fa",
    output:
        bam="data/clean/{sample}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        #module load Java/jdk1.8.0
        module load Java/1.8.0_192
        source activate gatk
        set -o pipefail
        gatk --java-options "-Xmx6G" SortSam \
        --INPUT={input.bam} \
        --OUTPUT=/dev/stdout \
        --SORT_ORDER=coordinate | \
        gatk --java-options "-Xmx6G" SetNmAndUqTags \
        --INPUT=/dev/stdin \
        --OUTPUT={output.bam} \
        --CREATE_INDEX=true \
        -R={input.gatk_ref}\
        --TMP_DIR="tmp4/{wildcards.sample}"
        '''

rule download_knowVar:
    output:
        dbSNP_vcf="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbSNP_vcf_index="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        knownIndels=["knowVar/Homo_sapiens_assembly38.known_indels.vcf.gz","knowVar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],
        knownIndels_indices=["knowVar/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi","knowVar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"],
        axiomPoly_resource_vcf="knowVar/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        axiomPoly_resource_vcf_index="knowVar/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
        hapmap_resource_vcf="knowVar/hapmap_3.3.hg38.vcf.gz",
        hapmap_resource_vcf_index="knowVar/hapmap_3.3.hg38.vcf.gz.tbi",
        omni_resource_vcf="knowVar/1000G_omni2.5.hg38.vcf.gz",
        omni_resource_vcf_index="knowVar/1000G_omni2.5.hg38.vcf.gz.tbi",
        one_thousand_genomes_resource_vcf="knowVar/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        one_thousand_genomes_resource_vcf_index="knowVar/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000
    shell:
        '''
        mkdir -p knowVar
        source activate gsutil-env
        gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz knowVar/.
        gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi knowVar/.
        gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz knowVar/.
        gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi knowVar/.
        gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf knowVar/.
        gsutil cp gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx knowVar/.
        gsutil cp gs://broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz knowVar/.
        gsutil cp gs://broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi knowVar/.
        gsutil cp gs://broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz knowVar/.
        gsutil cp gs://broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi knowVar/.
        gsutil cp gs://broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz knowVar/.
        gsutil cp gs://broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi knowVar/.
        gsutil cp gs://broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz knowVar/.
        gsutil cp gs://broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi knowVar/.
       
         #module load Java/jdk1.8.0
        #source activate gatk
        #gatk IndexFeatureFile -F knowVar/canis_familiaris_SNPs.vcf
        #gatk IndexFeatureFile -F knowVar/canis_familiaris_indels.vcf
        '''

rule BaseRecalib:
    input:
        #intervals="refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed",
        bam="data/clean/{sample}.bam",
        ref="refGenome/gatkIndex/genome.fa",
        dbSNP_vcf="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf",
        dbSNP_vcf_index="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
        knownIndels=["knowVar/Homo_sapiens_assembly38.known_indels.vcf.gz","knowVar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"],
        knownIndels_indices=["knowVar/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi","knowVar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi"]
    output:
        report="data/recalib/{sample}.txt",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        #module load Java/jdk1.8.0
        module load Java/1.8.0_192
        source activate gatk
        gatk --java-options "-Xmx7G" BaseRecalibrator \
        -R {input.ref} \
        -I {input.bam} \
        --use-original-qualities \
        -O {output.report} \
        --known-sites {input.dbSNP_vcf} \
        --known-sites {input.knownIndels[0]} \
        --known-sites {input.knownIndels[1]} 
        '''

        #-L {input.intervals} \
        #-ip 150

rule ApplyBQSR:
    input:
        #intervals="refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed",
        bam="data/clean/{sample}.bam",
        ref="refGenome/gatkIndex/genome.fa",
        report="data/recalib/{sample}.txt",
    output:
        bam="data/recalib/{sample}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        #module load Java/jdk1.8.0
        module load Java/1.8.0_192
        source activate gatk
        gatk --java-options "-Xmx7G" ApplyBQSR \
        -R {input.ref} \
        -I {input.bam} \
        -O {output.bam} \
        -bqsr {input.report} \
        --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 60 \
        --add-output-sam-program-record \
        --use-original-qualities 
        '''

        #-L {input.intervals} \
        #-ip 150

rule HaplotypeCaller_single:
    input:
        "refGenome/gatkIndex/genome.fa.fai", "refGenome/gatkIndex/genome.dict",
        ref="refGenome/gatkIndex/genome.fa",
        bam="data/recalib/{sample}.bam",
        intervals="refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed",
    output:
        gvcf="vc/hapCaller_single/{sample}.g.vcf",
        idx="vc/hapCaller_single/{sample}.g.vcf.idx",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 10,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 16
    threads:1
    shell:
        '''
        module load Java/1.8.0_192
        source activate gatk
        gatk --java-options "-Xmx15G" HaplotypeCaller \
        -R {input.ref} \
        -I {input.bam} \
        --emit-ref-confidence GVCF \
        --pcr-indel-model NONE \
        --read-filter OverclippedReadFilter \
        -L {input.intervals} \
        -ip 150 \
        -O {output.gvcf}
        '''


rule combineGVCFs:
    input:
        gvcfs=expand("vc/hapCaller_single/{sample}.g.vcf", sample=sample_ids),
        ref="refGenome/gatkIndex/genome.fa",
        intervals="refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed",
    output:
        "vc/combined_patch.g.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 24 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:1
    shell:
        '''
        inputs=()
        for gvcf in {input.gvcfs};do
            input_gvcf=$" -V "$gvcf;
            inputs+=($input_gvcf);
        done

        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx62G" CombineGVCFs \
        -R {input.ref} \
        ${{inputs[*]}} \
        -L {input.intervals} \
        -ip 150 \
        -O {output}
        '''

rule genotypeGVCFs:
    input:
        gvcf="vc/combined_patch.g.vcf",
        ref="refGenome/gatkIndex/genome.fa",
        dbSNP_vcf="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf",
        intervals="refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed",
    output:
        "vc/hapCaller_noped_raw.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx62G" GenotypeGVCFs \
        -R {input.ref} \
        -V {input.gvcf} \
        -D {input.dbSNP_vcf} \
        -G StandardAnnotation \
        --use-new-qual-calculator \
        --max-alternate-alleles 6 \
        -L {input.intervals} \
        -ip 150 \
        -O {output}
        '''


rule genotypeGVCFs_ped:
    input:
        gvcf="vc/combined_patch.g.vcf",
        ref="refGenome/gatkIndex/genome.fa",
        dbSNP_vcf="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf",
        intervals="refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed",
        ped="refGenome/ped_Paisas_seqIDs.ped"
    output:
        "vc/hapCaller_ped_raw.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 64
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx62G" GenotypeGVCFs \
        -R {input.ref} \
        -V {input.gvcf} \
        -D {input.dbSNP_vcf} \
        -G StandardAnnotation \
        --use-new-qual-calculator \
        --max-alternate-alleles 6 \
        -L {input.intervals} \
        -ip 150 \
        -ped {input.ped} \
        -O {output}
        '''

# ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
# than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
# Inbreeding Coefficient & family effect: https://software.broadinstitute.org/gatk/documentation/article?id=8032
# https://gatkforums.broadinstitute.org/gatk/discussion/2450/what-does-inbreedingcoeff-mean
# https://gatkforums.broadinstitute.org/gatk/discussion/2945/pedigree-analysis
# https://gatkforums.broadinstitute.org/gatk/discussion/7696/pedigree-ped-files
# https://software.broadinstitute.org/gatk/documentation/article?id=11016
rule HardFilterExcessHet:
    input:
        vcf="vc/hapCaller_{fam}_raw.vcf",
    output:
        "vc/hapCaller_{fam}_ExcessHet.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 2
    threads:1
    shell:
        '''
        excess_het_threshold="54.69"
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx3G" VariantFiltration \
        -V {input.vcf} \
        --filter-expression "ExcessHet > ${{excess_het_threshold}}" \
        --filter-name ExcessHet \
        -O {output}
        '''

rule MakeSitesOnlyVcf:
    input:
        vcf="vc/hapCaller_ped_ExcessHet.vcf",
    output:
        "vc/hapCaller_ped_sitesOnly.vcf",
    resources: 
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 2
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx3G" MakeSitesOnlyVcf \
        -I {input.vcf} \
        -O {output}
        '''


# Evaluating the quality of a variant callset: https://software.broadinstitute.org/gatk/documentation/article?id=6308
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.9.0/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.php
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.beta.6/org_broadinstitute_hellbender_tools_walkers_vqsr_VariantRecalibrator.php
# https://software.broadinstitute.org/gatk/documentation/article.php?id=39
# https://software.broadinstitute.org/gatk/documentation/article.php?id=1259
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
# error in the documentation: https://gatkforums.broadinstitute.org/gatk/discussion/12440/help-me-gatk4-vqsr-error
rule IndelsVariantRecalibrator:
    input:
        vcf="vc/hapCaller_ped_sitesOnly.vcf",
        mills_resource_vcf="knowVar/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
        axiomPoly_resource_vcf="knowVar/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
        dbSNP_vcf="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf",
    output:
        recal="vc/ADHD_ped.indels.recal",
        tranches="vc/ADHD_ped.indels.tranches",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 10,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 28
    threads:1
    shell:
        '''
        indel_recal_tranche_values=("100.0" "99.95" "99.9" "99.5" "99.0" "97.0" "96.0" "95.0" "94.0" "93.5" "93.0" "92.0" "91.0" "90.0")
        input_tranches=()
        for tranche in ${{indel_recal_tranche_values[*]}};do
            input_tranche=$" -tranche "$tranche;
            input_tranches+=($input_tranche);
        done

        indel_recal_ann_values=("FS" "ReadPosRankSum" "MQRankSum" "QD" "SOR")  ## "DP"
        input_anns=()
        for ann in ${{indel_recal_ann_values[*]}};do
            input_ann=$" -an "$ann;
            input_anns+=($input_ann);
        done

        module load Java/1.8.0_172
        source activate gatk
        module load R
        gatk --java-options "-Xmx24G" VariantRecalibrator \
        -V {input.vcf} \
        -O {output.recal} \
        --tranches-file {output.tranches} \
        --trust-all-polymorphic \
        ${{input_tranches[*]}} \
        ${{input_anns[*]}} \
        -mode INDEL \
        --max-gaussians 4 \
        --resource mills,known=false,training=true,truth=true,prior=12.0:{input.mills_resource_vcf} \
        --resource axiomPoly,known=false,training=true,truth=false,prior=10.0:{input.axiomPoly_resource_vcf} \
        --resource dbsnp,known=true,training=false,truth=false,prior=2.0:{input.dbSNP_vcf} \
        --rscript-file output.plots.indels.R
        '''

rule SNPsVariantRecalibrator:
    input:
        vcf="vc/hapCaller_ped_sitesOnly.vcf",
        hapmap_resource_vcf="knowVar/hapmap_3.3.hg38.vcf.gz",
        omni_resource_vcf="knowVar/1000G_omni2.5.hg38.vcf.gz",
        one_thousand_genomes_resource_vcf="knowVar/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbSNP_vcf="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf",
    output:
        recal="vc/ADHD_ped.snps.recal",
        tranches="vc/ADHD_ped.snps.tranches",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 10,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 28
    threads:1
    shell:
        '''
        snp_recal_tranche_values=("100.0" "99.95" "99.9" "99.8" "99.6" "99.5" "99.4" "99.3" "99.0" "98.0" "97.0" "90.0")
        input_tranches=()
        for tranche in ${{snp_recal_tranche_values[*]}};do
            input_tranche=$" -tranche "$tranche;
            input_tranches+=($input_tranche);
        done

        snp_recal_ann_values=("QD" "MQRankSum" "ReadPosRankSum" "FS" "MQ" "SOR")  ## "DP"
        input_anns=()
        for ann in ${{snp_recal_ann_values[*]}};do
            input_ann=$" -an "$ann;
            input_anns+=($input_ann);
        done

        module load Java/1.8.0_172
        source activate gatk
        module load R
        gatk --java-options "-Xmx24G" VariantRecalibrator \
        -V {input.vcf} \
        -O {output.recal} \
        --tranches-file {output.tranches} \
        --trust-all-polymorphic \
        ${{input_tranches[*]}} \
        ${{input_anns[*]}} \
        -mode SNP \
        --max-gaussians 6 \
        --resource hapmap,known=false,training=true,truth=true,prior=15:{input.hapmap_resource_vcf} \
        --resource omni,known=false,training=true,truth=true,prior=12:{input.omni_resource_vcf} \
        --resource 1000G,known=false,training=true,truth=false,prior=10:{input.one_thousand_genomes_resource_vcf} \
        --resource dbsnp,known=true,training=false,truth=false,prior=7:{input.dbSNP_vcf} \
        --rscript-file output.plots.snps.R
        '''

rule IndelsApplyRecalibration:
    input:
        vcf="vc/hapCaller_ped_ExcessHet.vcf",
        recal="vc/ADHD_ped.indels.recal",
        tranches="vc/ADHD_ped.indels.tranches",
    output:
        "vc/hapCaller_ped_ExcessHet_indelRecal.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        indel_filter_level=99.7
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx7G" ApplyVQSR \
        -V {input.vcf} \
        -O {output} \
        --recal-file {input.recal} \
        --tranches-file {input.tranches} \
        --truth-sensitivity-filter-level ${{indel_filter_level}} \
        --create-output-variant-index true \
        -mode INDEL
        '''

rule SNPsApplyRecalibration:
    input:
        vcf="vc/hapCaller_ped_ExcessHet_indelRecal.vcf",
        recal="vc/ADHD_ped.snps.recal",
        tranches="vc/ADHD_ped.snps.tranches",
    output:
        "vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        snp_filter_level=99.7
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx7G" ApplyVQSR \
        -V {input.vcf} \
        -O {output} \
        --recal-file {input.recal} \
        --tranches-file {input.tranches} \
        --truth-sensitivity-filter-level ${{snp_filter_level}} \
        --create-output-variant-index true \
        -mode SNP
        '''

rule CollectVariantCallingMetrics:
    input:
        vcf="vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf",
        dbSNP_vcf="knowVar/Homo_sapiens_assembly38.dbsnp138.vcf",
        ref_dict="refGenome/gatkIndex/genome.dict",
    output:
        "vc/ADHD_ped.variant_calling_detail_metrics", "vc/ADHD_ped.variant_calling_summary_metrics",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        metrics_filename_prefix="vc/ADHD_ped"
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx7G" CollectVariantCallingMetrics \
        --INPUT {input.vcf} \
        --DBSNP {input.dbSNP_vcf} \
        --SEQUENCE_DICTIONARY {input.ref_dict} \
        --OUTPUT ${{metrics_filename_prefix}} \
        --THREAD_COUNT 1 
        '''

rule split_vcf:
    input:
        vcf="vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf",
        gatk_ref="refGenome/gatkIndex/genome.fa",
    output:
        "vc/hapCaller_ped_ExcessHet_Recal_{var}.vcf",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 15,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000 * 8
    threads:1
    shell:
        '''
        module load Java/1.8.0_172
        source activate gatk
        gatk --java-options "-Xmx7G" SelectVariants \
        -R {input.gatk_ref} \
        -V {input.vcf} \
        -select-type {wildcards.var} \
        -O {output}
        '''

