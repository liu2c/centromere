#!/bin/bash

# environment need
gffread_soft=gffread
hisat2_path=hisat2-2.2.0
featureCounts=featureCounts
Rscript=/public/home/shangLG/software/anaconda3/envs/R4/bin/Rscript
samtools=samtools
stringtie=stringtie

DIR=$PWD

parameters=$(getopt -o r:f:g:t:h -l ref:,fqfile:,gff:,thread:,help -n "$0" -- "$@")
[ $? != 0 ] && exit 1
eval set -- "$parameters"   

while true ; do             
    case "$1" in
        -h|--help)     echo  "
FPKM.sh      usage : FPKM.sh -r NIP.fa -f file -g nip.gff

ref:            -r|--ref                Mandatory parameter. reference genome fasta.
fqfile:         -f|--fqfile             optional parameter. A file with three columns separated by TABs, containing the accession ID, the file path for the fastq_1 file, and the file path for the fastq_2 file.
gff:            -g|--gff                optional parameters. reference genome gene annotation.
thread:         [-t|--thread]           optional parameters. core nodes per task
help:           [-h|--help]             Nonoptional parameter

"                       ;exit 1;          shift 1 ;;   
        -r|--ref)     REF_fa=$2 ;         shift 2 ;;  
        -f|--fqfile)      FQ_File=$2 ;    shift 2 ;;   
        -g|--gff)      GFF=$2 ;    shift 2 ;;  
        -t|--thread)      
            case "$2" in
                "")     shift 2 ;;  
                *)      NT=$2 ;    shift 2 ;;   
            esac;;
        --)  break ;;       
        *) echo "wrong";exit 1;;
    esac
done

if [[ $NT == "" ]];then NT=5 ; fi

if [[ $REF_fa != "" && $FQ_File != "" && $GFF != "" && $quantify == "" ]];then
    # make index and transform gff to gtf
    mkdir -p $DIR/index
    PREFIX=$(echo $REF_fa | awk -F/ '{print $NF}' )
    if [[ $(ls $DIR/index/$PREFIX*ht2 |wc -l) -ne 8  ]];then
        $hisat2_path/hisat2-build -p 20 $REF_fa $DIR/index/$PREFIX
    fi
    
    if [[ ! -f  $DIR/index/$GFF.gtf ]];then
        $gffread_soft $GFF -T -o $DIR/index/$GFF.gtf
    fi

    # align
    cat $FQ_File | while read line;do

        read PREFIX1 FQ_clean_1 FQ_clean_2 <<< $line
        mkdir -p $DIR/$PREFIX1

        echo "
    if [[ ! -f  $DIR/$PREFIX1/$PREFIX1.hisat.bam ]];then

        if [[ ! -f  $DIR/$PREFIX1/$PREFIX1.hisat.sam ]];then
            $hisat2_path/hisat2 -p $NT --dta -x $DIR/index/$PREFIX -1 $FQ_clean_1 -2 $FQ_clean_2 --new-summary --summary-file $DIR/$PREFIX1/$PREFIX1.ht2.txt -S $DIR/$PREFIX1/$PREFIX1.hisat.sam
        fi

        $samtools sort -O bam -@ $NT -o $DIR/$PREFIX1/$PREFIX1.hisat.bam $DIR/$PREFIX1/$PREFIX1.hisat.sam
        ls $DIR/$PREFIX1/$PREFIX1.hisat.sam && rm $DIR/$PREFIX1/$PREFIX1.hisat.sam

    fi
    $stringtie $DIR/$PREFIX1/$PREFIX1.hisat.bam -o $DIR/$PREFIX1/$PREFIX1.gtf  -p $NT -G $DIR/index/$GFF.gtf -e -A $DIR/$PREFIX1/$PREFIX1.file
    
    "  > $DIR/$PREFIX1/$PREFIX1.hisat2.gtf.sh
    nohup sh $DIR/$PREFIX1/$PREFIX1.hisat2.gtf.sh &
    done
else
    echo "wrong parameters"
fi