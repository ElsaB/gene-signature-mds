# ------------------- #
# 2-class prediction
# ------------------- #

# model with restricted numbers of features [Q]

args <- commandArgs(trailingOnly = TRUE)

source("../corefunctions/predictors.R")
source("../corefunctions/Qpredictors.R")
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
mypred=args[4] # RF SVM L1 EN
q = args[5]
nbootstrap = args[6]

print(paste("number of cores used for parallel CV is", mc.cores))
print(paste("number of restricted features is", q))
print(paste("number of nbootstrap for SS is", nbootstrap))


nfolds=as.integer(nfolds)
nrepeats=as.integer(nrepeats)
mc.cores=as.integer(mc.cores)
mypred=as.character(mypred)
q = as.integer(q)
nbootstrap = as.integer(nbootstrap)

seed=234

# OUTPUT
out.path=paste0(mypred,"_",nfolds,"x",nrepeats,"_Q",q,"_B",nbootstrap)
system("mkdir -p results_cv_q")
out.path = paste("results_cv",out.path,sep="/")
system(paste("mkdir -p", out.path))

# ------------------ #
# RUN CV
# ------------------ #
print(paste("CV start"))

if (mypred=="RF") {
rescvq = runCV(mypredictor = QpredictorRF,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               q = q
               )
}
if (mypred=="SVM") {
rescvq = runCV(mypredictor = QpredictorSVM,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               q = q,
               intnfolds=3,
               nbootstrap=nbootstrap,
               beta = 0.2,
               nsamplemin=80,nsamplemax=90
               )
}
if (mypred=="L1") {
rescvq = runCV(mypredictor = QpredictorLogistic,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               intnfolds=3,
               alpha=1,
               nbootstrap=nbootstrap,
               beta = 0.2,
               nsamplemin=80,nsamplemax=90,
               measure="deviance"
               )
}
if (mypred=="EN") {
rescvq = runCV(mypredictor = QpredictorLogistic,
               labelTrain = response,
               patientTrain = design,
               nfolds = nfolds, nrepeats = nrepeats,
               seed = seed,
               mc.cores = mc.cores,
               intnfolds=3,
               alpha=0.5,
               nbootstrap=nbootstrap,
               beta = 0.2,
               nsamplemin=80,nsamplemax=90
               )
}
#if (mypred=="constant") {
#rescvq = runCV(mypredictor = predictorConstant,
#               labelTrain = response,
#               patientTrain = design,
#               nfolds = nfolds, nrepeats = nrepeats,
#               seed = seed,
#               mc.cores = mc.cores
#               )
#}

print(paste("CV done"))

# ------------------ #
# SAVE OBJECT
# ------------------ #

print(paste("save in",out.path))
save(rescvq, file=paste0(out.path,"/resCV.Rdata"))

z=sapply(rescvq, function(x) length(unique(x$ref)))
ishort = which(z==1)
if (length(ishort)>0) {
    rescvq.trick = rescvq[-ishort]
} else {
    rescvq.trick = rescvq
}

auc.cv = unlist(slot(evaluateCVwithROCR(rescvq.trick, measure="auc", x.measure="cutoff"), "y.values"))
print(paste("the mean AUC over", length(rescvq.trick), "CV experiments is", round(mean(auc.cv),3)))
acc.obj = evaluateCVwithROCR(rescvq.trick, measure="acc", x.measure="cutoff")
acc.max = sapply(slot(acc.obj, "y.values"), max)
fscore.obj = evaluateCVwithROCR(rescvq.trick, measure="f", x.measure="cutoff")
fscore.max = sapply(slot(fscore.obj, "y.values"), function(y) max(y, na.rm=T))
ddperf = data.frame(auc=auc.cv, acc_max=acc.max, fscore_max=fscore.max, predictor=mypred)
write.table(ddperf, file=paste0(out.path,"/resCV.txt"), sep="\t", quote=F, row.names=F)
