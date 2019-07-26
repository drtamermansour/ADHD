#!/bin/sh

#  ADAH.sh
#  
#
#  Created by Tamer Mansour on 6/10/19.
#  

## Litrature:
#Arcos-Burgos 2004. Attention-Deficit/Hyperactivity Disorder in a Population Isolate: Linkage to Loci at 4q13.2, 5q33.3, 11q22, and 17p11
#Arcos-Burgos 2010. A common variant of the latrophilin 3 gene, LPHN3, confers susceptibility to ADHD and predicts effectiveness of stimulant medication
#Arcos-Burgos 2012. Review: A common genetic network underlies substance use disorders and disruptive or externalizing disorders.
#Mastronardi 2016. Linkage and association analysis of ADHD endophenotypes in extended and multigenerational pedigrees from a genetic isolate
#Acosta 2016. ADGRL3 (LPHN3) variants are associated with a refined phenotype of ADHD in the MTA study

## Data downlowd
mkdir -p $SCRATCH/Tamer2/ADAH
cd $SCRATCH/Tamer2/ADAH
source activate workEnv1 ## Python 2.7.15
rclone copy "remote_box:ADHD Flagship BAM files" $SCRATCH/Tamer2/ADAH/data

## the bed file of the target sequence intervals: ADHD_capture_regions.nr.bed  (hg19 coordinates)
# The header line was deleted then
# coordinates were converted to hg38 using UCSC web tool: http://genome.ucsc.edu/cgi-bin/hgLiftOver then
# manually merge the small intervals into 5 large intervals

## The pedigree file ped_Paisas.xls
# The empty cells were filled with (?) then all (?) were replaced with NA
# file was saved in csv format & upload to the refGenome directory
cat refGenome/ped_Paisas.csv | tr '\r' '\n' | tr -s '\n' > refGenome/ped_Paisas_unix.csv
ls data/*.bam | awk -F"[/.]" '{print $2}' | sed 's/AM//' > refGenome/seq_samples
awk -F',' 'NR==FNR{c[$1]++;next};c[$8]' refGenome/seq_samples refGenome/ped_Paisas_unix.csv > refGenome/ped_Paisas_seq.csv
sed -i 's/NA/0/g' refGenome/ped_Paisas_seq.csv
#cat refGenome/ped_Paisas_seq.csv | awk -F"," '{print $1,$2,$3,$4,$5+1,$6}' > refGenome/ped_Paisas_seq.ped ## (0=male,1=female) changed to (1=male, 2=female) && (1=unaffected, 2=affected) kept as its
cat refGenome/ped_Paisas_seq.csv | awk -F"," '{print $1+"1000"$2+100,"AM"$8}' > refGenome/samples_ids
cat refGenome/ped_Paisas_seq.csv | awk -F"," '{print $1+"1000"$3+100, $3}' | sort | uniq > refGenome/fathers_ids
cat refGenome/ped_Paisas_seq.csv | awk -F"," '{print $1+"1000"$4+100, $4}' | sort | uniq > refGenome/mothers_ids
cat refGenome/fathers_ids refGenome/mothers_ids | sort | uniq > refGenome/parents_ids
awk '{print $1}' refGenome/samples_ids | grep -vFwf - refGenome/parents_ids > refGenome/unseq.parents_ids
cat refGenome/unseq.parents_ids >> refGenome/samples_ids
#awk 'NR==FNR{c[$1]=$2;next};{print $1,c[$1+"1000"$2+100],c[$1+"1000"$3+100],c[$1+"1000"$4+100],$5,$6}' refGenome/samples_ids refGenome/ped_Paisas_seq.ped > refGenome/ped_Paisas_seqIDs.ped
awk -F'[ ,]' 'NR==FNR{c[$1]=$2;next};{print $1,c[$1+"1000"$2+100],c[$1+"1000"$3+100],c[$1+"1000"$4+100],$5+1,$6}' refGenome/samples_ids refGenome/ped_Paisas_seq.csv > refGenome/ped_Paisas_seqIDs.ped  ## update ids and change (0=male,1=female) to (1=male, 2=female)

## preliminary QC
source activate workEnv1 ## Python 2.7.15
conda install -c bioconda qualimap
qualimap bamqc -bam data/AM5425.bam -outdir qualimap_results
## Check the header info
module load SAMtools/0.1.19
samtools view -H data/AM5425.bam > example_bam_header

## downlaod the software you need on the server
## multiQC (for Slurm configuration)
conda create -n multiQC multiqc==1.6

## bwa.kit
# Download the bwa-0.7.15 binary package
wget -O - https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download | tar xjf -
# update the download link of the reference genome
sed -i 's|url38=.*|url38="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"|' bwa.kit/run-gen-ref

## GATK
# https://software.broadinstitute.org/gatk/documentation/quickstart
# 1. Requirements
module load Java/jdk1.8.0 ## java -version ## java version "1.8.0_102"
# 2. Get GATK
# 2A. download: https://software.broadinstitute.org/gatk/download/
wget https://github.com/broadinstitute/gatk/releases/download/4.0.9.0/gatk-4.0.9.0.zip
unzip gatk-4.0.9.0.zip
cd gatk-4.0.9.0
# 2B. Add the path to the .bashrc file
echo "export PATH=$PATH:$(pwd)/gatk-4.0.9.0/gatk" >> $HOME/.bashrc
source $HOME/.bashrc
# 2B. run the conda env (I have installed miniconda already & I have miniconda3/bin in the PATH)
# https://software.broadinstitute.org/gatk/documentation/article?id=12836
conda update conda  ## conda --version ## conda 4.5.11
conda env create -n gatk -f gatkcondaenv.yml

## gsutil
conda create -n gsutil python=2.7
source activate gsutil
conda install -c flyem-forge gsutil-env
source activate gsutil-env
pip install oauth2client
###############
## Install Snakemake using conda & deploy PBS or slurm profile
bash snakemake_setup.sh

## Run Snakemake on slurm
# 1. add conda and the hpcc folder to you path
export PATH=$HOME/miniconda3/bin:$(pwd)/slurm:$PATH
# 2. turn on the environment
source activate snakemake
# 3. run by the submoit script
. slurm/submit-slurm.sh

## Note:
# I got error: /mnt/home/mansourt/miniconda3/etc/profile.d/conda.sh: line 55: PS1: unbound variable
# based on: https://github.com/conda/conda/issues/8186
# I changed line 55 in /mnt/home/mansourt/miniconda3/etc/profile.d/conda.sh
# old: ask_conda="$(PS1="$PS1" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix "$cmd" "$@")" || \return $?
# new: ask_conda="$(PS1="${PS1:-}" "$CONDA_EXE" $_CE_M $_CE_CONDA shell.posix "$cmd" "$@")" || \return $?

##################
## some temp code to find finised and failed jobs after Snakmake runs
mkdir error
grep error slurm-* | awk -F":" '{print $1}' | while read out;do mv $out error/.;done
grep output error/* | awk '{print $7}' | awk -F"," '{print $1}' | while read bam;do mv $bam error/.;done
grep output error/* | awk '{print $7}' | grep vcf | awk -F"," '{print $1}' | while read bam;do mv $bam error/.;done

mkdir failed
grep Failed slurm-* | awk -F":" '{print $1}' | while read out;do mv $out failed/.;done
grep output failed/* | awk '{print $7}' | while read bam;do mv $bam failed/.;done

mkdir logs/HaplotypeCaller_single
grep "Finished job" slurm-* | awk -F":" '{print $1}' | while read out;do mv $out logs/HaplotypeCaller_single/.;done

grep Failed slurm-* | awk -F":" '{print $1}' | awk -F"[-.]" '{print $2}' | while read job;do scancel $job;done
mv failed failed3
mkdir failed
grep Failed slurm-* | awk -F":" '{print $1}' | while read out;do mv $out failed/.;done
grep output failed/*.out | awk '{print $7}' | while read bam;do mv $bam failed/.;done
grep output failed/*.out | awk '{print $7}' | awk -F"[./]" '{print $3}' | while read temp;do rm -rf tmp2/$temp;done
############################
## Assess the coverage of the target regions
for bam in data/recalib/*.bam;do echo $bam; samtools view -c -F 2308 -L refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed $bam;done > target_coverage
for bam in data/recalib/*.bam;do echo $bam; samtools view -c -F 2308 $bam;done > total_coverage
cat target_coverage | paste - - > target_coverage2
cat total_coverage | paste - - > total_coverage2
paste <(head -n40 target_coverage2) <(head -n40 total_coverage2) | awk '{print $0,$2/$4}' > percentage_coverage
paste <(cat target_coverage2) <(cat total_coverage2) | awk '{print $0,$2/$4}' > percentage_coverage
cat percentage_coverage | awk '{if($5<0.7)print}' ## 87 sample
cat percentage_coverage | awk '{if($5<0.6)print}' ## 65 sample
cat percentage_coverage | awk '{if($5<0.5)print}' ## 36 sample
##############################
#stats
grep -v "^#" vc/hapCaller_raw_v1.vcf | wc -l
#52380
grep -v "^#" vc/hapCaller_noped_raw.vcf | wc -l  ## -D {input.dbSNP_vcf} -G StandardAnnotation --use-new-qual-calculator
grep -v "^#" vc/hapCaller_ped_raw.vcf | wc -l   ## add the ped file
#57509
grep -v "^#" vc/hapCaller_noped_raw.vcf | awk '{print $3}' | grep ^rs | wc -l
grep -v "^#" vc/hapCaller_ped_raw.vcf | awk '{print $3}' | grep ^rs | wc -l
#27365
grep -v "^#" vc/hapCaller_noped_ExcessHet.vcf | grep -v PASS | awk '{print $3,$7}' | wc -l
#2118
grep -v "^#" vc/hapCaller_ped_ExcessHet.vcf | grep -v PASS | awk '{print $3,$7}' | wc -l
#1374
grep -v "^#" vc/hapCaller_ped_ExcessHet_indelRecal.vcf | grep -v PASS | awk '{print $3,$7}' | wc -l
#1783
grep -v "^#" vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf | grep -v PASS | awk '{print $3,$7}' | wc -l
#17587
grep -v "^#" vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf | grep PASS | awk '{print $3,$7}' | wc -l
#39922
grep -v "^#" vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf | grep PASS | awk '{print $3}'| grep ^rs | wc -l
#26126
##################
mkdir -p vc/GWAS
GWAS=$(pwd)/vc/GWAS
snps="hapCaller_ped_ExcessHet_Recal_SNP"
indels="hapCaller_ped_ExcessHet_Recal_INDEL"
## a) prepare VCF files for analysis
## INDELS: select bi-alleleic indels, replace with ACTG character,and remove chrUn from VCF
awk '/#/{print;next}{if($5 !~ /,/){print}}' vc/$indels.vcf > $GWAS/$indels.monoAllel.vcf
awk 'BEGIN{FS="\t";OFS="\t"}/#/{print;next}{if(length($4)>1){$5!="A"?$4="A":$4="T";};if(length($5)>1){$4!="A"?$5="A":$5="T";};print;}' $GWAS/$indels.monoAllel.vcf > $GWAS/$indels.monoAllel_edit.vcf

## b) "new" prepare Plink input & create binary inputs & create file of alternative alleles
cd $GWAS
## merge filtered SNPs and indels in the fake snp format
module load VCFtools/0.1.15-Perl-5.26.1
vcf-concat ../$snps.vcf $indels.monoAllel_edit.vcf | vcf-sort > allSnp.vcf
## prepare Plink input
#vcftools --vcf allSnp.vcf --plink --out allSnp  ##input: 57509  & output: 51894 biallelic variants (The difference is multiallelic snps)
## create binary inputs
module load PLINK/1.9b_4.1-x86_64  #module load plink/1.9
#plink --file allSnp --make-bed --out allSnp.vcftools  ## 49759 variants and 422 people pass filters and QC.
#plink --vcf allSnp.vcf --biallelic-only 'strict' 'list' --make-bed --set-missing-var-ids "@:#" --out allSnp.plink ## 49759 variants and 422 people pass filters and QC. ## exactly perform as using the input of VCFtools except when all samples are homozygous mutant, the A1 will be 0 instead of the minor allele
plink --vcf allSnp.vcf --biallelic-only 'strict' 'list' --make-bed --set-missing-var-ids "@:#" --keep-allele-order --out allSnp.plinkA2 ## 49759 variants and 422 people pass filters and QC. # --keep-allele-order keep the refernce allele as A2
## create file of alternative alleles
#cat allSnp.vcf | awk 'BEGIN{FS="\t";OFS="\t";}/#/{next;}{{if($3==".")$3=$1":"$2;}print $3,$5;}'  > alt_alleles
ped="../../refGenome/ped_Paisas_seqIDs.ped"
cat $ped | awk '{print $1,$2}' > new_ids
cat allSnp.plinkA2.fam | awk '{print $1,$2,$2}' > old_ids
join -j 2 <(sort -k 2b,2 old_ids) <(sort -k 2b,2 new_ids) | awk '{print $1,$2,$4,$3}' > id_convert
plink -bfile allSnp.plinkA2 --update-ids id_convert --make-bed -out allSnp.plinkA2_FID
plink -bfile allSnp.plinkA2_FID --update-parents $ped --update-sex $ped 3 --make-bed -out allSnp.final # 135 founders and 287 nonfounders present.
rm allSnp.plinkA2.*
rm allSnp.plinkA2_FID.*
rm new_ids old_ids id_convert
module load R/3.5.1-X11-20180131
Rscript -e "install.packages('qqman', lib='~/R/v3.5.1/library', contriburl=contrib.url('http://cran.r-project.org/'))"

cat $ped | awk '{print $1,$2,$6}' > ADHD_list  ## 196 unaffected + 226 affected

## case control association analysis
#plink --bfile allSnp.final --make-pheno ADHD_list "2" --assoc --geno "0.5" --maf "0.01" --a2-allele allSnp.final.bim --adjust gc qq-plot --out allsnp.ADHD.asc
#Rscript -e 'library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
#'data=read.table("allsnp.ADHD.asc.assoc", header=TRUE); data=data[!is.na(data$P),];'\
#'pdf(paste("allsnp.ADHD.asc","pdf",sep="."), width=20, height=10);'\
#'manhattan(data, p = "P", col = c("blue4", "orange3"), cex = 0.6);'\
#'graphics.off();'



## Updated case control association analysis with more filters and appropriate test for family stratification
#--geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed, while --mind does the same for samples.
#--maf filters out all variants with minor allele frequency below the provided threshold (default 0.01),
#--hwe filters out all variants which have Hardy-Weinberg equilibrium exact test p-value below the provided threshold.
#--me filters out variants and samples/trios with Mendel error rates exceeding the given thresholds
# Given a case/control phenotype and a set of clusters defined via --within/--family, --mh computes a weighted average of the per-stratum odds ratios for each variant, along with a 1df chi-square statistic and p-value (for the null hypothesis that odds ratios for all strata are equal to 1)
# Variants are pruned based on the variance inflation factor (VIF), which recursively removes SNPs within a sliding window of 50 SNPs and 5 SNPs is used to shift the window at each step with VIF threathold equals 2 (The parameters for --indep). To read more about VIF: https://en.wikipedia.org/wiki/Variance_inflation_factor
plink --bfile allSnp.final --geno "0.1" --mind "0.1" --maf "0.05" --hwe 0.001 'midp' --me 0.05 0.1 'var-first' --a2-allele allSnp.final.bim --make-bed --keep-allele-order --out filtered.allSnp --indep 50 5 2
#1488 variants removed due to missing genotype data (--geno).
#--hwe: 301 variants removed due to Hardy-Weinberg exact test.
#37977 variants removed due to minor allele threshold(s)
#--me/--mendel: 11467 Mendel errors detected.
#46 variants and 14 people excluded.
#9947 variants and 408 people pass filters and QC.
#Pruning complete.  8727 of 9947 variants removed. ===> keep 1220 haplotype

plink --bfile filtered.allSnp --make-pheno ADHD_list "2" --family --mh --a2-allele filtered.allSnp.bim --adjust qq-plot --out filtered.allsnp.mh
Rscript -e 'library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
'data=read.table("filtered.allsnp.mh.cmh", header=TRUE); data=data[!is.na(data$P),];'\
'data$P=data$P * 1220;'\
'pdf(paste("filtered.allsnp.mh.cmh","pdf",sep="."), width=20, height=10);'\
'manhattan(data, p = "P", col = c("blue4", "orange3"), suggestiveline = -log(0.05,10), genomewideline = -log(0.01,10), cex = 0.6);'\
'graphics.off();'


## permutation and set-test
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
#ann=gencode.v31.annotation
#gunzip $ann.gtf.gz
#awk '{if($3=="gene")print;}' $ann.gtf | awk -F"\"" 'BEGIN{OFS="\t"}{print $1,$2}' | awk 'BEGIN{OFS="\t"}{print $1,$4,$5,$10}' > $ann.genes
#while read chr start end;do
#awk -v chr="$chr" -v start="$start" -v end="$end" '{if($1==chr && (($2 >= start && $2 <= end) || ($3 >= start && $3 <= end)))print;}' $ann.genes
#done < ../../refGenome/ADHD_capture_regions.nr.noheader.hg38_UCSCcovert_summerized.bed > target_gene_sets
#mkdir ../../plink_update && cd ../../plink_update
#wget http://s3.amazonaws.com/plink1-assets/dev/plink_linux_x86_64.zip
#unzip plink_linux_x86_64.zip
#cd $GWAS
#../../plink_update/plink --bfile filtered.allSnp --make-pheno ADHD_list "2" --make-set target_gene_sets --make-set-border 5 --family --mh 'perm' 'set-test' --a2-allele filtered.allSnp.bim --out filtered.allsnp.mh_set

plink --bfile filtered.allSnp --extract filtered.allSnp.prune.in --out pruned.filtered.allSnp --make-bed --keep-allele-order
plink --bfile pruned.filtered.allSnp --make-pheno ADHD_list "2" --family --mh --a2-allele pruned.filtered.allSnp.bim --adjust qq-plot --out pruned.filtered.allSnp.mh
# --adjust: Genomic inflation est. lambda (based on median chisq) = 1.
Rscript -e 'library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
'data=read.table("pruned.filtered.allSnp.mh", header=TRUE); data=data[!is.na(data$P),];'\
'pdf(paste("pruned.filtered.allSnp.mh","pdf",sep="."), width=20, height=10);'\
'manhattan(data, p = "P", col = c("blue4", "orange3"), cex = 0.6);'\
'graphics.off();'

## model association analysis
plink --bfile filtered.allSnp --make-pheno ADHD_list "2" --family --model --a2-allele filtered.allSnp.bim --out filtered.allsnp.mod
head -n1 filtered.allsnp.mod.model > filtered.allsnp.mod.model.sorted
cat filtered.allsnp.mod.model | grep -v NA | sort -k 10,10g >> filtered.allsnp.mod.model.sorted

## Chr17
awk '{if(NR==1)print;}{if($1==17)print;}' filtered.allsnp.mh.cmh | sort -k8,8g > chr17.filtered.allsnp.mh.cmh
awk '{if(NR==1)print;}{if($1==17 && $5=="ALLELIC")print;}' filtered.allsnp.mod.model.sorted > chr17.filtered.allsnp.mod.model.sorted.allel
awk '{if(NR==1)print;}{if($1==17 && $5=="GENO")print;}' filtered.allsnp.mod.model.sorted > chr17.filtered.allsnp.mod.model.sorted.geno

Rscript -e 'library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
'data=read.table("filtered.allsnp.mh.cmh", header=TRUE); data=data[!is.na(data$P),];'\
'data$P=data$P * 1220;'\
'pdf(paste("chr17","filtered.allsnp.mh.cmh","pdf",sep="."), width=20, height=10);'\
'manhattan(subset(data, CHR == "17"), p = "P", suggestiveline = -log(0.05,10), genomewideline = -log(0.01,10), annotatePval="0.05", annotateTop=TRUE, xlim = c(4160870,4660688));'\
'graphics.off();'

## The peak mutation
cd ../../
grep "^#CHR" vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf > rs4500776.vcf
grep -w rs4500776 vc/hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf >> rs4500776.vcf
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < rs4500776.vcf > rs4500776.tvcf
#tail -n+10 rs4500776.tvcf | awk -F":" '{print $1}' | sort > new

join -1 2 -2 1 <(sort -k 2b,2 $ped) <(tail -n+10 rs4500776.tvcf | sort -k 1b,1) | awk '{print $2,$1,$6, $7}' | sort -k1b,1 -k3b,3 > rs4500776.ped
cat rs4500776.ped | awk '{if($3==1)print $4}' | awk -F"[/:]" '{mut+=$1+$2;if($1!=".")c+=2}END{print mut, c-mut, mut/c}' # unaffected 138 250
cat rs4500776.ped | awk '{if($3==2)print $4}' | awk -F"[/:]" '{mut+=$1+$2;if($1!=".")c+=2}END{print mut, c-mut, mut/c}' # affected   230 220
## remember: in chr17.filtered.allsnp.mod.model.sorted.allel, rs4500776 was 222/212        135/241

cat rs4500776.ped | awk -F"[ /:]" '{a[$1]+=$4+$5;b[$1]+=2;}END{for(i in a)print i, a[i], b[i]}' > rs4500776.famFreq
cat rs4500776.famFreq | awk '{if($2==0){fam+=1;sum+=$3;print $0;}}END{print fam,sum}' ## 8 families with 28 individuals
cat rs4500776.ped | awk -F"[ /:]" '{a[$1]+=$4+$5;if($4!=".")b[$1]+=2;}END{for(i in a)print i, a[i], b[i]}' > rs4500776.famFreq2
cat rs4500776.famFreq2 | awk '{if($2==0){fam+=1;sum+=$3;print $0;}}END{print fam,sum}' ## 8 families with 27 genotyped individuals

cp rs4500776.ped temp
cat rs4500776.famFreq | awk '{if($2==0)print $1}' | while read famID;do grep -v ^$famID temp > temp2; cp temp2 temp;done
mv temp rs4500776_non0.ped
cat rs4500776_non0.ped | awk '{if($3==1)print $4}' | awk -F"[/:]" '{mut+=$1+$2;if($1!=".")c+=2}END{print c-mut, mut, mut/c}' # 228, 138
cat rs4500776_non0.ped | awk '{if($3==2)print $4}' | awk -F"[/:]" '{mut+=$1+$2;if($1!=".")c+=2}END{print c-mut, mut, mut/c}' # 188, 230
# https://www.socscistatistics.com/tests/chisquare/default2.aspx # p-value is .000001

#### check the variants in Chr4
awk '{if(NR==1)print;}{if($1==4)print;}' filtered.allsnp.mh.cmh | sort -k8,8g > chr4.filtered.allsnp.mh.cmh
grep "^#CHROM" ../hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf | awk '{print $1,$2,$3,$4,$5,$6,$7}' > chr4_rs.vcf
cat chr4_rs.list | grep -Fwf - ../hapCaller_ped_ExcessHet_indelRecal_snpRecal.vcf | awk '{print $1,$2,$3,$4,$5,$6,$7}' >> chr4_rs.vcf
head -n1 chr4.filtered.allsnp.mh.cmh > chr4_rs.filtered.allsnp.mh.cmh
cat chr4_rs.list | grep -Fwf - chr4.filtered.allsnp.mh.cmh >> chr4_rs.filtered.allsnp.mh.cmh
awk '{if(NR==1)print;}{if($1==4 && $5=="ALLELIC")print;}' filtered.allsnp.mod.model.sorted > chr4.filtered.allsnp.mod.model.sorted.allel
cat chr4.filtered.allsnp.mod.model.sorted.allel | grep -Fwf chr4_rs.list > chr4_rs.filtered.allsnp.mod.model.sorted.allel


head -n4 chr4.filtered.allsnp.mh.cmh
head -n4 chr4.filtered.allsnp.mh.cmh | awk '{print $2}' | grep -Fwf - chr4.filtered.allsnp.mod.model.sorted.allel
cat chr4_rs.vcf
cat chr4_rs.filtered.allsnp.mh.cmh
cat chr4_rs.filtered.allsnp.mod.model.sorted.allel
#############
## Association anlaysis using the original VCF 
mkdir -p vc/GWAS2
GWAS2=$(pwd)/vc/GWAS
cd $GWAS2 # copy the files from the local derive
gunzip div.vcf.gz snv.vcf.gz
## a) prepare VCF files for analysis
## INDELS: select bi-alleleic indels, replace with ACTG character,and remove chrUn from VCF
awk '/#/{print;next}{if($5 !~ /,/){print}}' div.vcf > div.monoAllel.vcf
awk 'BEGIN{FS="\t";OFS="\t"}/#/{print;next}{if(length($4)!=length($5))print;}' div.monoAllel.vcf > div.monoAllel_noSNPs.vcf
awk 'BEGIN{FS="\t";OFS="\t"}/#/{print;next}{if(length($4)>1){$5!="A"?$4="A":$4="T";};if(length($5)>1){$4!="A"?$5="A":$5="T";};print;}' div.monoAllel_noSNPs.vcf > div.monoAllel_edit.vcf
## b) "new" prepare Plink input & create binary inputs & create file of alternative alleles
## merge filtered SNPs and indels in the fake snp format
module load VCFtools/0.1.15-Perl-5.26.1
vcf-concat snv.vcf div.monoAllel_edit.vcf | vcf-sort > allSnp.vcf

## create binary inputs
module load PLINK/1.9b_4.1-x86_64  #module load plink/1.9
cp snv.vcf allSnp.vcf
grep "^##" allSnp.vcf > allSnp_editIDs.vcf
grep "^#CHROM" allSnp.vcf | sed 's/AM_/AM/g' | sed 's/_2//g' >> allSnp_editIDs.vcf
grep -v "^#" allSnp.vcf | awk '{if(length($1)<=5)print;}' >> allSnp_editIDs.vcf
plink --vcf allSnp_editIDs.vcf --biallelic-only 'strict' 'list' --make-bed --set-missing-var-ids "@:#" --keep-allele-order --out allSnp.plinkA2 ## 49759 variants and 422 people pass filters and QC. # --keep-allele-order keep the refernce allele as A2

ped="../../refGenome/ped_Paisas_seqIDs.ped"
cat $ped | awk '{print $1,$2}' > new_ids
cat allSnp.plinkA2.fam | awk '{print $1,$2,$2}' > old_ids
join -j 2 <(sort -k 2b,2 old_ids) <(sort -k 2b,2 new_ids) | awk '{print $1,$2,$4,$3}' > id_convert
plink -bfile allSnp.plinkA2 --update-ids id_convert --make-bed -out allSnp.plinkA2_FID
plink -bfile allSnp.plinkA2_FID --update-parents $ped --update-sex $ped 3 --make-bed -out allSnp.final

cat $ped | awk '{print $1,$2,$6}' > ADHD_list  ## 196 unaffected + 226 affected
## case control association analysis
plink --bfile allSnp.final --make-pheno ADHD_list "2" --assoc --geno "0.5" --maf "0.01" --a2-allele allSnp.final.bim --adjust gc qq-plot --out allsnp.ADHD.asc

module load R/3.5.1-X11-20180131
Rscript -e 'library(qqman, lib.loc = "/mnt/home/mansourt/R/v3.0.1/library");'\
'data=read.table("allsnp.ADHD.asc.assoc", header=TRUE); data=data[!is.na(data$P),];'\
'pdf(paste("allsnp.ADHD.asc","pdf",sep="."), width=20, height=10);'\
'manhattan(data, p = "P", col = c("blue4", "orange3"), cex = 0.6);'\
'graphics.off();'

# peak mutation
grep "^#CHR" allSnp.vcf > chr5_133201502.vcf
grep "^chr5" allSnp.vcf | grep -w 133201502 >> chr5_133201502.vcf
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < chr5_133201502.vcf > chr5_133201502.tvcf
join -1 2 -2 1 <(sort -k 2b,2 $ped) <(tail -n+10 chr5_133201502.tvcf | sort -k 1b,1) | awk '{print $2,$1,$6, $7}' | sort -k1b,1 -k3b,3 > chr5_133201502.ped
cat chr5_133201502.ped | awk '{if($3==1)print $4}' | awk -F"[/:]" '{mut+=$1+$2;if($1!=".")c+=2}END{print c-mut, mut, mut/c}' # 250, 138
cat chr5_133201502.ped | awk '{if($3==2)print $4}' | awk -F"[/:]" '{mut+=$1+$2;if($1!=".")c+=2}END{print c-mut, mut, mut/c}' # 220, 230


#####
## local assessment of pregenerated VCF files
cd /Users/drtamermansour/Desktop/FileZilla/max_maria_youssef/ADHD_exom_data
mkdir -p work_dir
cp ADHD*/*snv.vcf.gz work_dir/snv.vcf.gz
cp ADHD*/*div.vcf.gz work_dir/div.vcf.gz
gunzip work_dir/snv.vcf.gz
gunzip work_dir/div.vcf.gz
cd work_dir
grep "^#CHR" snv.vcf > rs4500776.vcf
grep ^chr17 snv.vcf | grep -w 4107098 >> rs4500776.vcf
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < rs4500776.vcf > rs4500776.tvcf
tail -n+10 rs4500776.tvcf | awk -F":" '{print $1}' | sed 's/AM_/AM/' | sort > original ## cp the new file from the server

###### 
#Note: how to run bwa.kit/bwa-postalt.js
bwa.kit/k8 bwa.kit/bwa-postalt.js refGenome/hs38DH.fa.alt altalt.sam > altalt_postalt.sam

