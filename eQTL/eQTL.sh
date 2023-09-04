#!/bin/bash


Rscript=/public/home/anaconda3/envs/R4/bin/Rscript
R1=eQTL-genotype.r
R2=matrixeQTL.r

DIR=$PWD
parameters=$(getopt -o v:f:c:k:g:p:t:s:h -l vcf:,fpkm:,cov:,keep:,gff:,singlegene:,help,prefix:,threshold:,Maf:,MaxMissing: -n "$0" -- "$@")
[ $? != 0 ] && exit 1
eval set -- "$parameters"  

while true ; do           
    case "$1" in
        -h|--help)     echo  "
eQTL.sh      usage : eQTL.sh -v NIP.vcf -f fpkm.txt -c $COV

vcf:            -v|--vcf                Mandatory parameter. reference genome fasta.
fpkm:           -f|--fpkm               optional parameter. expression file.
cov:            -c|--cov                optional parameter. covariant file of expression.
keep:           -k|--keep               optional parameter. keep sample for analysis
gff:            -g|--gff                optional parameter. gene annotation of reference
gene:           [-s|--singlegene]       optional parameter. for single gene eQTL.
prefix:         [-p|--prefix]           optional parameter. prefix for eQTL output.
threshold:      [-t|--threshold]        optional parameter. pvalue output threshold.
Maf:            [--Maf]                 optional parameter. vcftools'parameter (default 0.05).
MaxMissing:     [--MaxMissing]          optional parameter. vcftools'parameter (default 0.8).
help:           [-h|--help]             
"                       ;exit 1;          shift 1 ;;   
        -v|--vcf)       VCF=$2 ;        shift 2 ;;  
        -f|--fpkm)      FPKM=$2 ;       shift 2 ;;   
        -c|--cov)       COV=$2 ;        shift 2 ;;  
        -k|--keep)      KEEP=$2 ;       shift 2 ;;
        -g|--gff)       GFF=$2 ;        shift 2 ;; 
        -s|--singlegene)       GENE=$2 ;        shift 2 ;;  
        -p|--prefix)       PREFIX=$2 ;        shift 2 ;; 
        -t|--threshold)       THRESHOLD=$2 ;        shift 2 ;; 
        --Maf)          Maf=$2;         shift 2 ;;  
        --MaxMissing)   MaxMissing=$2;  shift 2 ;;  
        --)  break ;;       
        *) echo "wrong";exit 1;;
    esac
done

if [[ $Maf == "" ]];then Maf=0.05; fi
if [[ $MaxMissing == "" ]];then MaxMissing=0.8; fi
if [[ $THRESHOLD == "" ]];then 
    if [[ $GENE != "" ]];then
        THRESHOLD=1e-4
    else
        THRESHOLD=1
    fi
fi

Vcf_Filter(){ # VCF KEEP Maf MaxMissing
    if [[ ! -f ${VCF}_${KEEP}${Maf}-${MaxMissing}.vcf ]];then
        vcftools --vcf $VCF --keep $KEEP --recode --recode-INFO-all --stdout --maf $Maf  --max-missing $MaxMissing > ${VCF}_${KEEP}${Maf}-${MaxMissing}.vcf
    fi
}

Check_Sample(){ # VCF FPKM COV
    # etermine whether the gene is in the annotation file or the expression quantity file, or exit script if it is not
    if [[ $GENE != "" ]];then
        test -z "$( grep $GENE $FPKM)" && echo -e "$GENE not in $FPKM" && exit 1
        test -z "$(grep $GENE $GFF)" && echo -e "$GENE not in $GFF" && exit 1
    fi
    # Original vcf Indicates the KEEP filtered VCF
    head -100 $VCF | grep "#CHROM" | tr "\t" "\n" |sed 1,9d | sort -V > title
    head -1 $COV | tr "\t" "\n" > peer.title
    head -1 $FPKM | tr "\t" "\n"  > ex.titile
    join  title peer.title | join - ex.titile > comm.title
}

# 1. Genotype location information
Genotype_Pos(){
    if [[ ! -f $VCF-genotype.pos ]];then
        { echo "SNP CHR POS"; grep -v "#"  $VCF |awk '{print $1"_"$2,$1,$2}'; } > $VCF-genotype.pos
    fi
}

# 2. Genotype information (012 format)
Matrix(){
    if [[ ! -f $VCF.matrix ]];then
        vcftools --vcf $VCF --012  --out $VCF
        $Rscript $R1 -d $VCF.012 -o $VCF.matrix
        sed -i -r 's/X[.]+CHROM[.]*/#CHROM/g' $VCF.matrix
    fi
}


# 3. PCA
PCA(){
    if [[ ! -f ${VCF}.cov ]];then
        plink --vcf $VCF --out ${VCF} --make-bed --allow-extra-chr
        plink --allow-extra-chr --threads 5 -bfile ${VCF} --pca 20 --out ${VCF}
        # combine Eigenvalue vectors (top 5) and expression covariates (20) 
        cut -f $(for f in `cat comm.title`;do grep -n $f  peer.title ;done | awk -F: '{printf "%s,",$1}'|sed 's/,$//g') $COV | awk -v i=1 'NR==1{print "peer""\t"$0}NR!=1{print "V"i"\t"$0; let i++}' > $VCF.peer

        { head -1 $VCF.peer; cut -d " " --complement -f1-2 ${VCF}.eigenvec | awk '{for(i=1;i<=NF;i++){a[FNR,i]=$i}}END{for(i=1;i<=NF;i++){for(j=1;j<=FNR;j++){printf a[j,i]" "}print ""}}' | tr " " "\t"|head -5|awk -v i=1 '{print "PC"i"\t"$0;let i++}';  sed -n '2,21p' $VCF.peer; } | sed 's/\t$//g' > ${VCF}.cov
    fi
}

# 4. Expression data
Expression_Data(){
    if [[ $GENE == "" ]];then
        if [[ ! -f $VCF.fpkm ]];then
            cut -f 1,$(for f in `cat comm.title`;do grep -n $f  ex.titile ;done | awk -F: '{printf "%s,",$1}'|sed 's/,$//g') $FPKM > $VCF.fpkm
        fi
    else
        if [[ ! -f $VCF.$GENE.fpkm ]];then
            { head -1 $FPKM; grep $GENE $FPKM; }| cut -f 1,$(for f in `cat comm.title`;do grep -n $f  ex.titile ;done | awk -F: '{printf "%s,",$1}'|sed 's/,$//g') > $VCF.$GENE.fpkm
        fi
    fi
}

# 5. Gene location file
Gene_Position(){
    if [[ $GENE == "" ]];then
        if [[ ! -f $VCF.pos.gff ]];then
            #grep -v "#" $GFF| grep -w "gene"| awk -F "[=;\t]" 'BEGIN{print "geneid""\t""chr""\t""left""\t""right"}{print $10"\t"$1"\t"$4"\t"$5}' > $VCF.pos.gff
            grep -v "#" $GFF| grep -w "gene"| awk -F "[=;\t]" '{print $10"\t"$1"\t"$4"\t"$5}' | sort -V -k 1,1 > $FPKM.pos
            { echo -e "geneid\tchr\tleft\tright"; cut -f1 $FPKM | sed 1d | sort -V -k 1,1|join -t $'\t' - $FPKM.pos -1 1 -2 1; } > $VCF.pos.gff
        fi

    else
        if [[ ! -f $VCF.$GENE.pos.gff ]];then
            #grep -v "#" $GFF| grep -w "gene"| awk -F "[=;\t]" 'BEGIN{print "geneid""\t""chr""\t""left""\t""right"}{print $10"\t"$1"\t"$4"\t"$5}' > $VCF.pos.gff
            grep -v "#" $GFF| grep $GENE | grep -w "gene"| awk -F "[=;\t]" 'BEGIN{print "geneid""\t""chr""\t""left""\t""right"}{print $10"\t"$1"\t"$4"\t"$5}'  > $VCF.$GENE.pos.gff
        fi
    fi
}

main(){
    cd $DIR
    Vcf_Filter
    TMP=$(echo $VCF)
    VCF=$( echo ${TMP}_${KEEP}${Maf}-${MaxMissing}.vcf )
    Check_Sample
    Genotype_Pos
    Matrix
    PCA
    Expression_Data
    Gene_Position

    if [[ $PREFIX == "" ]];then
        if [[ $GENE == "" ]];then
            $Rscript $R2 -v $VCF.matrix -e $VCF.fpkm -s $VCF-genotype.pos -g $VCF.pos.gff -c $VCF.cov -t $THRESHOLD
        else
            $Rscript $R2 -v $VCF.matrix -e $VCF.$GENE.fpkm -s $VCF-genotype.pos -g $VCF.$GENE.pos.gff -c $VCF.cov -t $THRESHOLD
        fi
    else
        if [[ $GENE == "" ]];then
            $Rscript $R2 -v $VCF.matrix -e $VCF.fpkm -s $VCF-genotype.pos -g $VCF.pos.gff -c $VCF.cov -p $PREFIX -t $THRESHOLD
        else
            $Rscript $R2 -v $VCF.matrix -e $VCF.$GENE.fpkm -s $VCF-genotype.pos -g $VCF.$GENE.pos.gff -c $VCF.cov -p $PREFIX -t $THRESHOLD
        fi
    fi

}


if [[ $VCF != "" && $FPKM != "" && $COV != "" && $KEEP != "" && $GFF != "" ]];then
    main 
else
    echo "wrong parameters"
fi
