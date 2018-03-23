# ------------------- #
# 2-class prediction
# stability selection
# ------------------- #

#source("../corefunctions/predictors.R")
#source("../corefunctions/evaluate.R")

## import the training set
print("Get the training set ready")
source("./getReadySubgroup.R")
print("Done")

# ------------------------- #
# PARAMETERS
# ------------------------- #
mypred = "RF"
mypred=as.character(mypred)

# OUTPUT
out.path = "results_impRF"
system(paste("mkdir -p", out.path))

# ------------------ #
# RUN IMPORTANCE FOR RF
# ------------------ #

library(randomForest)

if (mypred=="RF") {

    model.rf = randomForest(design, factor(response), ntree = 1000, importance=TRUE)

    dd.imp.rf = as.data.frame( importance(model.rf) )

}


# ------------------ #
# SAVE OBJECT
# ------------------ #

dd.imp.rf$predictor=mypred

print(paste("save in",out.path))

write.table(dd.imp.rf, file=paste0(out.path,"/impRF.txt"), sep="\t", quote=F, row.names=F)
