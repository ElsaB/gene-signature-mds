# ---------------------------------- #
# 2-CLASS CLASSIFICATION PREDICTORS
# ---------------------------------- #
#
# In this file we have predictor functions, that take at least 3 arguments:
#	- patientTrain : n*p matrix of descriptors for n training patients
#	- patientTest : m*p matrix of descriptors for m test patients
#	- labelTrain : n given labels on the n training patients
# and that output:
#       - class : the vector of predicted classes [1 or 2]
#       - proba : the vector of predicted probabilities for class 1
#       - genes [optional] : the names of the selected genes for feature-selection algorithm

# NB: probability = probability of belonging to the positive class (eg. "driver")

# Qpredictors are predictors where you select Q features on the training set
# and you retrain a model using those Q features
# and then you test


# ---------------------------------- #
# random classification
# ---------------------------------- #
predictorConstant <- function(patientTrain, patientTest, labelTrain, positiveClass="1") {

    class.random = factor(rep(positiveClass, nrow(patientTest)))
    proba.random = rep(0, nrow(patientTest))

    return(list(class=class.random, proba=proba.random, gene=NULL))
}


# ---------------------------------- #
# random forest classification
# ---------------------------------- #
library(randomForest)
QpredictorRF <- function(patientTrain, patientTest, labelTrain, positiveClass="1", ntree=500, q=100) {

    # Train a random forest model
    model.rf = randomForest(patientTrain, factor(labelTrain), ntree = ntree, importance=TRUE)
    # Consider importance
    imp.rf = importance(model.rf)[,"MeanDecreaseGini"]
    # Select the Q most important feature
    isort.rf = sort(imp.rf, decreasing=T, index.r=T)$ix[1:q]
    # Re-train with limited features
    model.rf.q = randomForest(patientTrain[,isort.rf], factor(labelTrain), ntree = ntree)
    # Test
    pred.rf = 1 - predict(model.rf.q, patientTest[,isort.rf], type="prob")[,positiveClass]
    class.rf = predict(model.rf.q, patientTest[,isort.rf], type="class")

    return(list(class=class.rf , proba=pred.rf, gene=NULL))
}


# ---------------------------------- #
# SVM classification
# ---------------------------------- #
library(kernlab)
QpredictorSVM <- function(patientTrain, patientTest, labelTrain, positiveClass="1", intnfolds=5, Clist=2^seq(-10,15), q=100, ...) {

    # Run Stability Selection on the training set
    imp.svm = runSS(mypredictor=predictorSVM, patientTrain=patientTrain, labelTrain=labelTrain, ...)
    # Consider importance
    # Select the Q most important feature
    isort.svm = sort(imp.svm$freq, decreasing=T, index.r=T)$ix[1:q]
    # Re-train with limited features
    list.svp = vector(mode="list",length=length(Clist))
    ii = 0
    for (cc in Clist) {
        ii = ii+1
        list.svp[[ii]] <- ksvm(patientTrain[,isort.svm],labelTrain,type="C-svc",kernel="rbf",C=cc,cross=intnfolds) #, prob.model=TRUE)
    }
    mycv = sapply(list.svp, function(x) cross(x))
    svp = list.svp[[which.min(mycv)]]
    #proba.svm = predict(svp, patientTest, type="proba")[,positiveClass]
    proba.svm = predict(svp, patientTest[,isort.svm], type="decision")[,1]
    class.svm = factor(predict(svp, patientTest[,isort.svm], type="response"))
    gene.svm = colnames(patientTrain[,isort.svm])[alphaindex(svp)[[1]]]

    return(list(class=class.svm , proba=proba.svm, gene=gene.svm))

}

# ---------------------------------- #
# Logistic classification
# ---------------------------------- #
library(glmnet)
library(parallel)
QpredictorLogistic <- function(patientTrain, patientTest, labelTrain, positiveClass="1", q=100, measure="deviance", intnfolds=5, alpha=1, ...) {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    # alpha=1/2 --> elastic net

    # Run Stability Selection on the training set
    imp.log = runSS(mypredictor=predictorLogistic, patientTrain=patientTrain, labelTrain=labelTrain, ...)
    # Consider importance
    # Select the Q most important feature
    isort.log = sort(imp.log$freq, decreasing=T, index.r=T)$ix[1:q]
    # Re-train with limited features
    cvfit = cv.glmnet(patientTrain[,isort.log], factor(labelTrain), family="binomial", type.measure=measure, alpha=alpha, nfolds=intnfolds, parallel=TRUE)
    # Test
    proba.log = predict(cvfit, newx = patientTest[,isort.log], s = "lambda.min", type="response")[,positiveClass]
    class.log = factor(predict(cvfit, newx = patientTest[,isort.log], s = "lambda.min", type="class")[,positiveClass])
    aa.log = predict(cvfit, newx = patientTest[,isort.log], s = "lambda.min", type="coef")[,1]
    gene.log = names(aa.log[aa.log>0][-1])

    return(list(class=class.log , proba=proba.log, gene=gene.log))
}
