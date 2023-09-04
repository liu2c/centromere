suppressMessages(library(data.table))
suppressMessages(library(getopt))

options(scipen=200) # nolint
options(stringsAsFactors = F)

command=matrix(c(
  'help',  'h', 0, 'logic', 'help information',
  'data',  'd', 1, 'character', 'inputfile: indispensable , 012 vcf and current dir has .pos file and .indv file Or vcf in pca',
  'output','o', 2, 'character', 'inputfile: optional , the name of output file'
),byrow = T, ncol = 5)

args = getopt(command)
## help information
if (!is.null(args$help)) {
  cat(paste(getopt(command, usage = T), "\n"))
  q(status=1)
}

x <- args$data
out.name = args$output


# transform to 012 matrix
data <- data.table::fread(x, header = F, sep = "\t")
indv = data.table::fread(paste0(x,".indv"), header = F, sep = "\t")
pos = data.table::fread(paste0(x,".pos"), header = F, sep = "\t")
pos$info = paste0(pos$V1,"_",pos$V2)

data = t(as.matrix(data))
data = data[-1,]

add = data.frame("`#CHROM`"=paste0(pos$V1,"_",pos$V2))

colnames(data) = indv$V1
data = cbind(add,data)

if( is.character(out.name) ){
  data.table::fwrite(data,out.name,sep="\t")
}else{
  data.table::fwrite(data, paste0(x,".matrix"),sep="\t")
}
