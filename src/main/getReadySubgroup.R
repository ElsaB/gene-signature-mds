genes.paper = c("HOXB6","HOXB3","GPSM1","NYNRIN","SLMO1","ZNF662","SPINK2","MAPK12","HOXA3","AMN",
"ROBO4","MYO5C","PRODH","UGCG","OLFM4","LCN2","CAST","WARS","ELL2",
"CYSTM1","CYBB",
"JHDM1D","OLR1","CAMP","CRISP3")

FilterLog2cpmMatrix <- function(X, minsample0=5, valuemax=4, variancemin=0.1) {

    # filter very low express genes
    iX = apply( X, 1, function(x) ( sum(x>0)>minsample0 & max(x)>valuemax ) )
    X = X[iX,]
    # filter the very low variance genes
    myvar = apply(X,1,var)
    X = X[myvar>variancemin,]

    return(X)
}

# -------------------------- #
# DATA
# -------------------------- #

bm = read.table("../../data/BMMNC_log2cpm.txt", sep="\t", stringsAsFactors=F)
colnames(bm) = as.vector(sapply(colnames(bm), function(x) strsplit(x,split="_BMMNC")[[1]][1]))
bm = bm[,-ncol(bm)]
cd34 = read.table("../../data/CD34_log2cpm.txt", sep="\t", stringsAsFactors=F)
cd34 = cd34[,-ncol(cd34)]

clinical = read.csv("../../data/Clinical_info.csv", stringsAsFactors=F)
clinical$ID = paste0("PV",clinical$ID)

# training set = subgroup created with cd34 and should be predicted with bmmnc
training.samples = intersect(colnames(bm),colnames(cd34))

print(paste("We have",length(training.samples),"training samples"))

#cd34  = FilterLog2cpmMatrix(cd34, minsample0=5, valuemax=2, variancemin=0.1)
bm  = FilterLog2cpmMatrix(bm, minsample0=5, valuemax=2, variancemin=0.1)

design = t(bm[ , colnames(bm) %in% training.samples ])
response = clinical$subgroup[match(rownames(design),clinical$ID)]

rm(cd34)
rm(bm)

print(paste("There are",sum(!genes.paper %in% colnames(design)),"genes found in the paper but not even in the design"))

print(dim(design))
print(table(response))
