# ------------------- #
# 2-class prediction
# stability selection
# ------------------- #
args <- commandArgs(trailingOnly = TRUE)

# stability selection on the entire training set
# to have an idea of the robustness/stability of the feature overall

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
nbootstrap = args[1]
beta = args[2]
mc.cores=args[3]
mypred=args[4] # SVM L1 EN

print(paste("number of cores used for parallel SS is", mc.cores))

nbootstrap=as.integer(nbootstrap)
beta=as.numeric(beta)
mc.cores=as.integer(mc.cores)
mypred=as.character(mypred)

seed=434

# OUTPUT
out.path=paste0(mypred,"_",nbootstrap,"x",beta)
system("mkdir -p results_ss")
out.path = paste("results_ss",out.path,sep="/")
system(paste("mkdir -p", out.path))

# ------------------ #
# RUN CV
# ------------------ #
print(paste("SS start"))

if (mypred=="SVM") {
res.ss = runSS(mypredictor = predictorSVM,
               labelTrain = response,
               patientTrain = design,
               nbootstrap = nbootstrap,
               beta = beta,
               mc.cores = mc.cores,
               intnfolds=3
               )
}

if (mypred=="L1") {
res.ss = runSS(mypredictor = predictorLogistic,
               labelTrain = response,
               patientTrain = design,
               nbootstrap = nbootstrap,
               beta = beta,
               mc.cores = mc.cores,
               intnfolds=3,
               alpha=1
               )
}
if (mypred=="EN") {
res.ss = runSS(mypredictor = predictorLogistic,
               labelTrain = response,
               patientTrain = design,
               nbootstrap = nbootstrap,
               beta = beta,
               mc.cores = mc.cores,
               intnfolds=3,
               alpha=0.5
               )
}


print(paste("SS done"))

# ------------------ #
# SAVE OBJECT
# ------------------ #

res.ss$predictor=mypred
print(paste("save in",out.path))
save(res.ss, file=paste0(out.path,"/resSS.Rdata"))

write.table(res.ss, file=paste0(out.path,"/resSS.txt"), sep="\t", quote=F, row.names=F)
