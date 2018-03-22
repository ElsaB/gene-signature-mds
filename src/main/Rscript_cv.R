# ------------------- #
# 2-class prediction
# ------------------- #
args <- commandArgs(trailingOnly = TRUE)

source("../corefunctions/predictors.R")
source("../corefunctions/evaluate.R")

## import the training set
print("Get the training set ready")
source("./getReadySubgroup.R")
print("Done")

# ------------------------- #
# PARAMETERS
# ------------------------- #

# CV PARAMETERS
nfolds=args[1]
nrepeats=args[2]
mc.cores=args[3]
mypred=args[4] # RF SVM L1 L2 EN

print(paste("number of cores used for parallel CV is", mc.cores))

nfolds=as.integer(nfolds)
nrepeats=as.integer(nrepeats)
mc.cores=as.integer(mc.cores)
mypred=as.character(mypred)

seed=234

# OUTPUT
out.path=paste0(mypred,"_",nfolds,"x",nrepeats)
system("mkdir -p results_cv")
out.path = paste("results_cv",out.path,sep="/")
system(paste("mkdir -p", out.path))

# ------------------ #
# RUN CV
# ------------------ #
print(paste("CV start"))

if (mypred=="RF") {
rescv = runCV(mypredictor = predictorRF,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores)
}
if (mypred=="SVM") {
rescv = runCV(mypredictor = predictorSVM,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               intnfolds=3
               )
}
if (mypred=="L1") {
rescv = runCV(mypredictor = predictorLogistic,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               intnfolds=3,
               alpha=1
               )
}
if (mypred=="L2") {
rescv = runCV(mypredictor = predictorLogistic,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               intnfolds=3,
               alpha=0
               )
}
if (mypred=="EN") {
rescv = runCV(mypredictor = predictorLogistic,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               intnfolds=3,
               alpha=0.5
               )
}
if (mypred=="constant") {
rescv = runCV(mypredictor = predictorConstant,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores
               )
}

print(paste("CV done"))

# ------------------ #
# SAVE OBJECT
# ------------------ #

print(paste("save in",out.path))
save(rescv, file=paste0(out.path,"/resCV.Rdata"))

auc.cv = unlist(slot(evaluateCVwithROCR(rescv, measure="auc", x.measure="cutoff"), "y.values"))
print(paste("the mean AUC over", length(rescv), "CV experiments is", round(mean(auc.cv),3)))
acc.obj = evaluateCVwithROCR(rescv, measure="acc", x.measure="cutoff")
acc.max = sapply(slot(acc.obj, "y.values"), max)
fscore.obj = evaluateCVwithROCR(rescv, measure="f", x.measure="cutoff")
fscore.max = sapply(slot(fscore.obj, "y.values"), function(y) max(y, na.rm=T))
ddperf = data.frame(auc=auc.cv, acc_max=acc.max, fscore_max=fscore.max, predictor=mypred)
write.table(ddperf, file=paste0(out.path,"/resCV.txt"), sep="\t", quote=F, row.names=F)
