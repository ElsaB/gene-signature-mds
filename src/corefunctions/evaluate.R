library('parallel')
library('ROCR')

# CROSS-VALIDATION
runCV <- function(mypredictor, labelTrain, patientTrain, nfolds=5, nrepeats=10, seed=2396, mc.cores=1, ...) {
    # function that run mypredictor on a CV setting
    # output a list of size the number of CV experiments (eg 50)
    # each entry of the list is itself with 2 entries: "proba" and "ref"
    # "proba" contains the vector of predicted probabilities
    # "ref" contains the labels of the fold test set

    # Set random number generator seed
    set.seed(seed)

    # Make folds
    n = nrow(patientTrain)
    folds <- list()
    for (i in seq(nrepeats)) {
        folds <- c(folds,split(sample(seq(n)), rep(1:nfolds, length = n)))
    }
    nexp = length(folds) # the total number CV of experiments

    print("start CV")
    rescv = mclapply(seq(nexp),
                   FUN=function(iexp) {
                       cat(".")
                       vTrain = patientTrain[-folds[[iexp]],,drop=F]
                       vTest = patientTrain[folds[[iexp]],,drop=F]
                       lTrain = labelTrain[-folds[[iexp]]]
                       lTest = labelTrain[folds[[iexp]]]

                       pred = mypredictor(patientTrain=vTrain, patientTest=vTest, labelTrain=lTrain, ...)
                       
                       res.fold = list(class=pred$class, proba=pred$proba, gene=pred$gene, ref=lTest)
                       return(res.fold)
                   },
                   mc.cores=mc.cores
                   )

    return(rescv)

}

# EVALUATE CV
evaluateCVwithROCR <- function(resCV, measure="acc", x.measure="cutoff") {
    # example of usage: 
    #    - measure="tpr", x.measure="fpr" ---> ROC curve
    #    - measure="auc", x.measure="cutoff" ---> AUC
    list.pred = lapply(resCV, function(c) c$proba) # TO MAKE SURE
    #list.pred = lapply(resCV, function(c) 1-c$proba)
    list.ref = lapply(resCV, function(c) c$ref)
    pred.obj = prediction(list.pred, list.ref)
    perf.obj = performance(pred.obj, measure=measure, x.measure=x.measure)
    return(perf.obj)
}


# LEAVE-ONE-OUT
LeaveOneOut <- function(mypredictor, labelTrain, patientTrain, mc.cores=1, ...) {

    res.proba = mclapply(1:length(labelTrain),
                   FUN=function(j) {
                       cat(".")
                       mypred = mypredictor(patientTrain=patientTrain[-j,,drop=F], patientTest=patientTrain[j,,drop=F], labelTrain=labelTrain[-j], ...) 
                       return(mypred)
                   },
                   mc.cores=mc.cores
                   )

    return(unlist(res.proba))

}
