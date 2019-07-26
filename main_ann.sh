##ann_v38
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/gencode.v31.annotation.gtf.gz
ann=gencode.v31.annotation
gunzip $ann.gtf.gz
awk '{if($3=="gene")print;}' $ann.gtf > $ann.gene.gtf

## ann_v37
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_31/GRCh37_mapping/gencode.v31lift37.annotation.gtf.gz
ann=gencode.v31lift37.annotation
gunzip $ann.gtf.gz
awk '{if($3=="gene")print;}' $ann.gtf > $ann.gene.gtf
listID="demontis_2019"
#listID="arcos_2004"
interval="2500000";
while read chr BP;do
awk -v chr="$chr" -v BP="$BP" -v intv=$interval '{if($1==chr && (($4 > BP-intv && $4 < BP+intv) || ($5 > BP-intv && $5 < BP+intv) || ($4 < BP && $5 > BP)))print;}' $ann.gene.gtf > $chr.$BP.$interval
cat $chr.$BP.$interval | awk -F "\"" '{print $2}' | awk -F"." '{print $1}' > $chr.$BP.$interval.locus
grep -Fwf $chr.$BP.$interval.locus GN.txt | sort -k3,3n > $chr.$BP.$interval.GN
done < $listID.txt
cat *.GN | sort -k3,3n | head > $listID.$interval.head
awk '{print $2}' $listID.$interval.head | grep -Fwf - *.GN | awk '{print $2,$1,$3,$4,$5}' | sort > $listID.peak.$interval.head
awk '{print $2}' $listID.$interval.head | grep -Fwf - $ann.gene.gtf | awk -F "[\t;]" '{print $1,$4,$5,$9,$10,$11}' | awk '{print $5,$1,$2,$3,$7,$9}' | sort > $listID.ann.$interval.head
paste $listID.peak.$interval.head $listID.ann.$interval.head | sort -k3,3n > $listID.complete.$interval.head
mkdir ${listID}_$interval;
mv *.$interval* ${listID}_$interval/.
cat ${listID}_$interval/$listID.complete.$interval.head

##########################
##########################

* match the candidate genes with the genome network list
## interval = 1M
ENSG00000125414 chr17.10647130.500000.GN:MYH2 1397 2.177 2.177  "ENSG00000125414.19_6" chr17 10424465 10453017 "protein_coding" "MYH2"
ENSG00000272736 chr17.10647130.500000.GN:RP11-799N11.1 2188 1.872 1.872 "ENSG00000272736.5_7" chr17 10286449 10441179 "lncRNA" "AC005323.2"
ENSG00000230561 chr5.126842699.500000.GN:LINC01183 2690 1.737 1.737     "ENSG00000230561.4_3" chr5 127039081 127277326 "protein_coding" "CCDC192"
ENSG00000250956 chr5.126842699.500000.GN:CTB-88F18.3 3372 1.587 1.587   "ENSG00000250956.1_5" chr5 126515448 126515592 "processed_pseudogene" "AC011416.2"
ENSG00000245937 chr5.126842699.500000.GN:LINC01184 4007 1.477 1.477     "ENSG00000245937.8_6" chr5 127274844 127418864 "lncRNA" "LINC01184"
ENSG00000250375 chr4.60021455.500000.GN:RP11-340A13.3 4293 1.427 1.427  "ENSG00000250375.2_7" chr4 59912738 59940974 "lncRNA" "AC097501.1"
ENSG00000160588 chr11.117697731.500000.GN:MPZL3 4515 1.393 1.393        "ENSG00000160588.10_3" chr11 118097405 118123080 "protein_coding" "MPZL3"
ENSG00000242692 chr17.10647130.500000.GN:RPS27AP1 6038 1.204 1.204      "ENSG00000242692.1_2" chr17 10162079 10162549 "processed_pseudogene" "RPS27AP1"
ENSG00000173926 chr5.126842699.500000.GN:MARCH3 6644 1.142 1.142        "ENSG00000173926.6_3" chr5 126203406 126366250 "protein_coding" "MARCH3"
ENSG00000206579 chr8.56130095.500000.GN:XKR4 7685 1.043 1.043   "ENSG00000206579.8_3" chr8 56014949 56454613 "protein_coding" "XKR4"

## interval = 5M
ENSG00000110395 chr11.117697731.2500000.GN:CBL 149 4.212 4.212  "ENSG00000110395.6_4" chr11 119076986 119184636 "protein_coding" "CBL"
ENSG00000224077 chr11.117697731.2500000.GN:AP000936.4 256 3.741 3.741   "ENSG00000224077.1_6" chr11 116969703 116978886 "lncRNA" "AP000936.1"
ENSG00000118058 chr11.117697731.2500000.GN:KMT2A 286 3.628 3.628        "ENSG00000118058.21_5" chr11 118307179 118397547 "protein_coding" "KMT2A"
ENSG00000023287 chr8.56130095.2500000.GN:RB1CC1 365 3.359 3.359 "ENSG00000023287.13_5" chr8 53535018 53658403 "protein_coding" "RB1CC1"
ENSG00000154096 chr11.117697731.2500000.GN:THY1 385 3.284 3.284 "ENSG00000154096.13_5" chr11 119288088 119295695 "protein_coding" "THY1"
ENSG00000186174 chr11.117697731.2500000.GN:BCL9L 527 2.991 2.991        "ENSG00000186174.12_5" chr11 118764584 118796317 "protein_coding" "BCL9L"
ENSG00000260254 chr11.117697731.2500000.GN:AP000997.2 863 2.548 2.548   "ENSG00000260254.1_3" chr11 115509281 115517680 "lncRNA" "AP000997.2"
ENSG00000230716 chr11.117697731.2500000.GN:KRT8P7 1027 2.406 2.406      "ENSG00000230716.3_2" chr11 119473587 119475018 "processed_pseudogene" "KRT8P7"
ENSG00000095139 chr11.117697731.2500000.GN:ARCN1 1062 2.38 2.38 "ENSG00000095139.14_3" chr11 118443105 118473748 "protein_coding" "ARCN1"
ENSG00000110367 chr11.117697731.2500000.GN:DDX6 1095 2.363 2.363        "ENSG00000110367.13_4" chr11 118618472 118661873 "protein_coding" "DDX6"

## interval = 10M
ENSG00000262990 chr17.10647130.5000000.GN:CTC-281F24.2 53 5.229 5.229   "ENSG00000262990.1_6" chr17 6653004 6653297 "processed_pseudogene" "AC004706.2"
ENSG00000110395 chr11.117697731.5000000.GN:CBL 149 4.212 4.212  "ENSG00000110395.6_4" chr11 119076986 119184636 "protein_coding" "CBL"
ENSG00000224077 chr11.117697731.5000000.GN:AP000936.4 256 3.741 3.741   "ENSG00000224077.1_6" chr11 116969703 116978886 "lncRNA" "AP000936.1"
ENSG00000181222 chr17.10647130.5000000.GN:POLR2A 261 3.724 3.724        "ENSG00000181222.16_5" chr17 7387685 7417935 "protein_coding" "POLR2A"
ENSG00000118058 chr11.117697731.5000000.GN:KMT2A 286 3.628 3.628        "ENSG00000118058.21_5" chr11 118307179 118397547 "protein_coding" "KMT2A"
ENSG00000023287 chr8.56130095.5000000.GN:RB1CC1 365 3.359 3.359 "ENSG00000023287.13_5" chr8 53535018 53658403 "protein_coding" "RB1CC1"
ENSG00000154096 chr11.117697731.5000000.GN:THY1 385 3.284 3.284 "ENSG00000154096.13_5" chr11 119288088 119295695 "protein_coding" "THY1"
ENSG00000261996 chr17.10647130.5000000.GN:CTC-281F24.1 395 3.253 3.253  "ENSG00000261996.1_6" chr17 6556783 6558328 "lncRNA" "AC004706.1"
ENSG00000134852 chr4.60021455.5000000.GN:CLOCK 525 2.995 2.995  "ENSG00000134852.15_4" chr4 56294070 56413076 "protein_coding" "CLOCK"
ENSG00000186174 chr11.117697731.5000000.GN:BCL9L 527 2.991 2.991        "ENSG00000186174.12_5" chr11 118764584 118796317 "protein_coding" "BCL9L"

* match the genes in demontis_2019 loci with genome network list (smaller interval)
## interval = 100k
ENSG00000213917 chr2.215181889.50000.GN:RPL5P8 2265 1.849 1.849 "ENSG00000213917.2_4" chr2 215145362 215146176 "processed_pseudogene" "RPL5P8"
ENSG00000144451 chr2.215181889.50000.GN:SPAG16 3468 1.568 1.568 "ENSG00000144451.19_3" chr2 214149103 215275225 "protein_coding" "SPAG16"
ENSG00000229444 chr1.44184192.50000.GN:RP11-184I16.4 4434 1.403 1.403   "ENSG00000229444.1_6" chr1 44175063 44193014 "lncRNA" "AL451062.1"
ENSG00000271327 chr12.89760744.50000.GN:RP11-1109F11.3 6035 1.205 1.205 "ENSG00000271327.1_6" chr12 89761584 89763078 "lncRNA" "AC010201.2"
ENSG00000126091 chr1.44184192.50000.GN:ST3GAL3 8365 0.985 0.985 "ENSG00000126091.20_6" chr1 44171495 44396837 "protein_coding" "ST3GAL3"
ENSG00000066135 chr1.44184192.50000.GN:KDM4A 8905 0.939 0.939   "ENSG00000066135.13_5" chr1 44115820 44171189 "protein_coding" "KDM4A"
ENSG00000251182 chr4.31151456.50000.GN:RP11-617I14.1 11737 0.74 0.74    "ENSG00000251182.2_6" chr4 31172635 31213300 "lncRNA" "LINC02497"
ENSG00000271259 chr12.89760744.50000.GN:RP11-1109F11.5 13678 0.624 0.624        "ENSG00000271259.1_6" chr12 89765597 89766136 "lncRNA" "AC010201.1"
ENSG00000250705 chr5.87854395.50000.GN:CTC-470C15.1 13710 0.622 0.622   "ENSG00000250705.1_5" chr5 87898598 87898935 "processed_pseudogene" "AC008526.1"
ENSG00000197585 chr2.215181889.50000.GN:AC107218.3 14391 0.587 0.587    "ENSG00000197585.10_7" chr2 215106400 215548970 "lncRNA" "AC068051.1"

## interval = 0.5M
ENSG00000257594 chr12.89760744.250000.GN:GALNT4 1429 2.158 2.158        "ENSG00000257594.4_4" chr12 89913189 89918573 "protein_coding" "GALNT4"
ENSG00000213917 chr2.215181889.250000.GN:RPL5P8 2265 1.849 1.849        "ENSG00000213917.2_4" chr2 215145362 215146176 "processed_pseudogene" "RPL5P8"
ENSG00000259768 chr16.72578131.250000.GN:RP5-991G20.1 3113 1.641 1.641  "ENSG00000259768.6_8" chr16 72699022 72856680 "lncRNA" "AC004943.2"
ENSG00000144451 chr2.215181889.250000.GN:SPAG16 3468 1.568 1.568        "ENSG00000144451.19_3" chr2 214149103 215275225 "protein_coding" "SPAG16"
ENSG00000229444 chr1.44184192.250000.GN:RP11-184I16.4 4434 1.403 1.403  "ENSG00000229444.1_6" chr1 44175063 44193014 "lncRNA" "AL451062.1"
ENSG00000271327 chr12.89760744.250000.GN:RP11-1109F11.3 6035 1.205 1.205        "ENSG00000271327.1_6" chr12 89761584 89763078 "lncRNA" "AC010201.2"
ENSG00000070961 chr12.89760744.250000.GN:ATP2B1 6507 1.156 1.156        "ENSG00000070961.15_2" chr12 89981826 90103077 "protein_coding" "ATP2B1"
ENSG00000126091 chr1.44184192.250000.GN:ST3GAL3 8365 0.985 0.985        "ENSG00000126091.20_6" chr1 44171495 44396837 "protein_coding" "ST3GAL3"
ENSG00000066135 chr1.44184192.250000.GN:KDM4A 8905 0.939 0.939  "ENSG00000066135.13_5" chr1 44115820 44171189 "protein_coding" "KDM4A"
ENSG00000081189 chr5.87854395.250000.GN:MEF2C 10149 0.843 0.843 "ENSG00000081189.15_7" chr5 88012934 88200074 "protein_coding" "MEF2C"

hg19: POC1B (chr12:89813495-89919824), GALNT4 (chr12:89913189-89918573), POC1B-GALNT4 (chr12:89913185-89920039), POC1B-AS1 (chr12:89918371-89941782)
hg38: POC1B (chr12:89419718-89526047), GALNT4 (chr12:89519412-89524796), POC1B-GALNT4 (chr12:89519408-89526262)

Bogari et al. Whole Exome Sequencing: Novel Genetic Polymorphisms in Saudi Arabian Attention Deficit Hyperactivity Disorder (ADHD) Children. Natural Science, 2019. ==> rs2230283 (C/T; chr12:89916811) is one of 6 missense variants found in 3 out of 4 patients screened for 222 candidate genes


