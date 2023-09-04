#!/public/home/miniconda3/envs/R3.5/bin/Rscript
suppressMessages(library(MatrixEQTL))
suppressMessages(library(data.table))
suppressMessages(library(getopt))

command=matrix(c(
  'help', 'h', 0, 'logic', 'help information',
  'vcf', 'v', 1, 'character', 'vcf file',
  'expresion_data', 'e', 1 ,'character', 'FPKM data',
  'snps_location', 's', 1, 'character', 'snps location file',
  'gene_location', 'g', 1, 'character', 'gene location file',
  'covariates_file', 'c', 1, 'character', 'covariates file',
  'output_prefix', 'p', 1, 'character', 'output prefix',
  'output_pvalue_threshold', 't', 1, 'numeric', 'output pvalue threshold'

),byrow = T, ncol = 5)

args = getopt(command)
## help information
if (!is.null(args$help)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}


useModel = modelLINEAR

#5 FILES preparation
#1 SNP_genetype_file
SNP_file_name=args$vcf
#2 Gene_file
expression_file_name=args$expresion_data
#3 SNP_loc_file
snps_location_file_name=args$snps_location
#4 Gene_loc_file
gene_location_file_name=args$gene_location
#5 Covariates_file
covariates_file_name=args$covariates_file

#6 threshold value
pvs=args$output_pvalue_threshold

# SNP_file_name="NIP_T2T.SV.vcf_Os0.05-0.8.vcf.matrix"
# expression_file_name="NIP_T2T.SV.vcf_Os0.05-0.8.vcf.fpkm"
# snps_location_file_name="NIP_T2T.SV.vcf_Os0.05-0.8.vcf-genotype.pos"
# gene_location_file_name="NIP.T2T.pos.gff"
# covariates_file_name="NIP_T2T.SV.vcf_Os0.05-0.8.vcf.cov"

prefix=args$output_prefix
#output_file_name = tempfile()
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

#pvOutputThreshold_cis = 2e-2;
#pvOutputThreshold_tra = 1e-2;
pvOutputThreshold_cis = pvs
pvOutputThreshold_tra = pvs

#pvOutputThreshold = 1e-20
errorCovariance = numeric();

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

cisDist = 2e4;

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 10000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 10000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels

if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}


me = Matrix_eQTL_main(
snps = snps, 
gene = gene, 
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel, 
errorCovariance = errorCovariance, 
verbose = TRUE, 
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos, 
genepos = genepos,
cisDist = cisDist,
pvalue.hist = "qqplot",
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);


#noFDRsaveMemory = FALSE);
data.table::fwrite(me$cis$eqtls, file=paste0(prefix,"_cis.txt"),sep="\t")
data.table::fwrite(me$tran$eqtls, file=paste0(prefix,"_trans.txt"),sep="\t")
