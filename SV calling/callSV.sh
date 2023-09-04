#!/bin/bash

# evironment need
samtools=samtools
vcftools=vcftools
minimap2=minimap2
cuteSV=bincuteSV
python3=python3
SURVIVOR=SURVIVOR


DIR=$PWD

parameters=$(getopt -o d:r:m:h -l data:,REF:,Maf:,help -n "$0" -- "$@")
[ $? != 0 ] && exit 1
eval set -- "$parameters"   # 将$parameters设置为位置参数

while true ; do             # 循环解析位置参数
    case "$1" in
        -h|--help)     echo  "
callsv.sh      usage : callsv.sh -d path.file -r ref.path

data:        -d|--data              必选参数；三列文件，第一列样本名，第二列和第三列fastq数据（加上路径）,TAB分隔符
REF:         -r|--ref               必选参数；参考基因组序列（加上路径）
Maf:        [-m|--Maf]              可选参数；提取vcf文件后的过滤参数 ( 默认 0.05 ; eg: -m 0.05 | --Maf=0.05 | ' ' )
help:       [-h|--help]             非选项参数
"                       ;exit 1;          shift 1 ;;   # 非选项参数的选项
        -d|--data)     data=$2 ;         shift 2 ;;   # 带参数的选项
        -r|--REF)      REF=$2 ;    shift 2 ;;   # 带参数的选项
        -m|--Maf)      Maf=$2;         shift 2 ;;  # 给了可选参数
        --)  break ;;       # 开始解析非选项类型的参数，break后，它们都保留在$@中
        *) echo "wrong";exit 1;;
    esac
done

# minimap2tobam_single-SVcalling(){ # Sam.name 1_clean_fastq.gz 2_clean_fastq.gz dir -- bam gt.vcf gt.flt.vcf
#     sam=$1
#     clean_q1=$2
#     clean_q2=$3
#     cd $4
#     #minimap2 比对得到排序后的bam并建索引
#     if [[ $clean_q2 == "NA" ]];then
#         $minimap2 -t 10 -ax map-ont $REF $clean_q1 | \
#         $samtools view --threads 10 -bS -q 20 - | \
#         $samtools sort --threads 10 - -o $sam.sort.mq20.bam
#     else
#         $minimap2 -t 10 -ax -p 0.0001 map-ont $REF $clean_q1 $clean_q2 | \
#         $samtools view --threads 10 -bS -q 20 - | \
#         $samtools sort --threads 10 - -o $sam.sort.mq20.bam
#     fi
#     $samtools index $sam.sort.mq20.bam
#     #cuteSV 预 call SV（群call）
#     $python3 $cuteSV $sam.sort.mq20.bam $REF $sam.cuteSV.gt.vcf ./ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -l 50 -L 1000000 --genotype -S $sam -t 10
#     #过滤
#     $vcftools --vcf $sam.cuteSV.gt.vcf --remove-filtered IMPRECISE --keep-filtered PASS --recode --recode-INFO-all --stdout > $sam.cuteSV.flt.gt.vcf

#     #统计特征
#     $samtools stats $sam.sort.mq20.bam > $sam.stats
#     grep ^SN $sam.stats| cut -f 2-|grep error|cut -f 2 > $sam.error
# }



cat $data |while read line ;do
    read sam q1 q2 <<< $(echo $line | awk '{print $1,$2,$3}')
    #生成bam和pre-VCF
    mkdir -p $DIR/$sam 

    case $? in 
        0)
            echo  "minimap2tobam(){ # Sam.name 1_clean_fastq.gz 2_clean_fastq.gz dir -- bam gt.vcf gt.flt.vcf
    sam=\$1
    clean_q1=\$2
    clean_q2=\$3
    cd \$4
    #minimap2 比对得到排序后的bam并建索引
    if [[ \$clean_q2 == \"NA\" ]];then
        $minimap2 -t 10 -ax map-ont $REF \$clean_q1 | $samtools view --threads 10 -bS -q 20 - | $samtools sort --threads 10 - -o \$sam.sort.mq20.bam
    else
        $minimap2 -t 10 -ax -p 0.0001 map-ont $REF \$clean_q1 \$clean_q2 |  $samtools view --threads 10 -bS -q 20 - | $samtools sort --threads 10 - -o \$sam.sort.mq20.bam
    fi
    $samtools index \$sam.sort.mq20.bam
    #cuteSV 预 call SV（群call）
    $python3 $cuteSV \$sam.sort.mq20.bam $REF \$sam.cuteSV.gt.vcf ./ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -l 50 -L 1000000 --genotype -S \$sam -t 10
    #过滤
    $vcftools --vcf \$sam.cuteSV.gt.vcf --remove-filtered IMPRECISE --keep-filtered PASS --recode --recode-INFO-all --stdout > \$sam.cuteSV.flt.gt.vcf

    #统计特征
    $samtools stats \$sam.sort.mq20.bam > \$sam.stats
    grep ^SN \$sam.stats| cut -f 2-|grep error|cut -f 2 > \$sam.ErRoR
}
minimap2tobam $sam $q1 NA $DIR/$sam"  >  $DIR/$sam/$sam.bam.sh
            cd $DIR/$sam && qsub -q shangLG -l nodes=1:ppn=10,mem=100g $sam.bam.sh 
    ;;

        *)
        echo "error in mkdir $DIR/$sam"
    esac
done

#pre-merger
echo  "
all_N=$(cat $data |wc -l)
cd $DIR
find -name \"*cuteSV.flt.gt.vcf\" | awk -F \"/\" 'OFS=\"/\"{\$1=\"$DIR\"; print \$0}' > \${all_N}.path
$SURVIVOR merge \${all_N}.path 1000 1 1 -1 -1 50 \${all_N}.raw.minimap2.cuteSV.vcf
" > $DIR/premerge.sh


#force-calling
echo "
awk '{print \$1}' $data |while read sam ;do
    cd $DIR
    all_N=$(cat $data |wc -l)
    echo \"
$python3 $cuteSV $DIR/\$sam/\$sam.sort.mq20.bam $REF $DIR/\$sam/\$sam.cuteSV.flt.gt.vcf $DIR/\$sam/ -Ivcf $DIR/\${all_N}.raw.minimap2.cuteSV.vcf --min_support 10 --sample \$sam\" > $DIR/\$sam/\$sam.force_calling.sh
done" > $DIR/force_calling.sh


#re-merge
echo  "
all_N=$(cat $data |wc -l)
cd $DIR
find -name \"*cuteSV.flt.gt.vcf\" | awk -F \"/\" 'OFS=\"/\"{\$1=\"$DIR\"; print \$0}' > \${all_N}.path.re
$SURVIVOR merge \${all_N}.path.re 1000 1 1 -1 -1 50 \${all_N}.force-callinged.minimap2.cuteSV.vcf
" > $DIR/remerge.sh
