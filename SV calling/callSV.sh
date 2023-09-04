#!/bin/bash

# evironment need
samtools=samtools
vcftools=vcftools
minimap2=minimap2
cuteSV=bincuteSV
python3=python3
SURVIVOR=SURVIVOR


DIR=$PWD

parameters=$(getopt -o d:r:m:h -l data:,REF:,help -n "$0" -- "$@")
[ $? != 0 ] && exit 1
eval set -- "$parameters"   

while true ; do            
    case "$1" in
        -h|--help)     echo  "
callsv.sh      usage : callsv.sh -d path.file -r ref.path

data:        -d|--data              Mandatory parameter. A file with three columns separated by TABs, containing the accession ID, the file path for the fastq_1 file, and the file path for the fastq_2 file.
REF:         -r|--ref               Mandatory parameter.  reference genome fasta.
help:       [-h|--help]             Nonoptional parameter.
"                       ;exit 1;          shift 1 ;;  
        -d|--data)     data=$2 ;         shift 2 ;;   
        -r|--REF)      REF=$2 ;    shift 2 ;;  
        --)  break ;;       
        *) echo "wrong";exit 1;;
    esac
done

cat $data |while read line ;do
    read sam q1 q2 <<< $(echo $line | awk '{print $1,$2,$3}')
    # generate BAM and pre-calling
    mkdir -p $DIR/$sam 

    case $? in 
        0)
            echo  "minimap2tobam(){ # Sam.name 1_clean_fastq.gz 2_clean_fastq.gz dir -- bam gt.vcf gt.flt.vcf
    sam=\$1
    clean_q1=\$2
    clean_q2=\$3
    cd \$4
    # use minimap2 to align and make index
    if [[ \$clean_q2 == \"NA\" ]];then
        $minimap2 -t 10 -ax map-ont $REF \$clean_q1 | $samtools view --threads 10 -bS -q 20 - | $samtools sort --threads 10 - -o \$sam.sort.mq20.bam
    else
        $minimap2 -t 10 -ax -p 0.0001 map-ont $REF \$clean_q1 \$clean_q2 |  $samtools view --threads 10 -bS -q 20 - | $samtools sort --threads 10 - -o \$sam.sort.mq20.bam
    fi
    $samtools index \$sam.sort.mq20.bam
    # cuteSV pre-calling SV
    $python3 $cuteSV \$sam.sort.mq20.bam $REF \$sam.cuteSV.gt.vcf ./ --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 -l 50 -L 1000000 --genotype -S \$sam -t 10
    # filter PRECISE variants
    $vcftools --vcf \$sam.cuteSV.gt.vcf --remove-filtered IMPRECISE --keep-filtered PASS --recode --recode-INFO-all --stdout > \$sam.cuteSV.flt.gt.vcf

    # calulate 
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
