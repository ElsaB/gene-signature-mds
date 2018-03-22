# -------------------------------------- #
# REDUCE DIMENSION OF EXPRESSION MATRIX
# GENTLE FILTERING
# -------------------------------------- #
# filter low express genes
# filter very low variance genes

# CD34
cd34 = read.table("../../data/CD34_log2cpm.txt", sep="\t", stringsAsFactors=F, header=T)
cd34 = cd34[,-ncol(cd34)]
# BMMNC
bm = read.table("../../data/BMMNC_log2cpm.txt", sep="\t", stringsAsFactors=F, header=T)
colnames(bm) = as.vector(sapply(colnames(bm), function(x) strsplit(x,split="_BMMNC")[[1]][1]))
bm = bm[,-ncol(bm)]


FilterLog2cpmMatrix <- function(X, minsample0=5, valuemax=4, variancemin=0.1) {

    # filter very low express genes
    iX = apply( X, 1, function(x) ( sum(x>0)>minsample0 & max(x)>valuemax ) )
    X = X[iX,]
    # filter the very low variance genes
    myvar = apply(X,1,var)
    X = X[myvar>variancemin,]

    return(X)
}

cd34  = FilterLog2cpmMatrix(cd34)
bm  = FilterLog2cpmMatrix(bm)

# SAVE
system("mkdir -p ../../data/tmp_reduce")
write.table(cd34, sep="\t", quote=F, file="../../data/tmp_reduce/CD34_log2cpm.txt")
write.table(bm, sep="\t", quote=F, file="../../data/tmp_reduce/BMMNC_log2cpm.txt")
