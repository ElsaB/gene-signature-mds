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
predictorRF <- function(patientTrain, patientTest, labelTrain, positiveClass="1", ntree=500) {

    # Train a random forest model
    model.rf = randomForest(patientTrain, factor(labelTrain), ntree = ntree)
    # Make predictions and output probabilities
    pred.rf = 1 - predict(model.rf, patientTest, type="prob")[,positiveClass]
    class.rf = predict(model.rf, patientTest, type="class")

    return(list(class=class.rf , proba=pred.rf, gene=NULL))
}


# ---------------------------------- #
# SVM classification
# ---------------------------------- #
library(kernlab)
predictorSVM <- function(patientTrain, patientTest, labelTrain, positiveClass="1", intnfolds=5, Clist=2^seq(-10,15)) {

    list.svp = vector(mode="list",length=length(Clist))
    ii = 0
    for (cc in Clist) {
        ii = ii+1
        list.svp[[ii]] <- ksvm(patientTrain,labelTrain,type="C-svc",kernel="rbf",C=cc,cross=intnfolds) #, prob.model=TRUE)
    }
    mycv = sapply(list.svp, function(x) cross(x))
    svp = list.svp[[which.min(mycv)]]

    #proba.svm = predict(svp, patientTest, type="proba")[,positiveClass]
    proba.svm = predict(svp, patientTest, type="decision")[,1]
    class.svm = factor(predict(svp, patientTest, type="response"))
    gene.svm = colnames(patientTrain)[alphaindex(svp)[[1]]]

    return(list(class=class.svm , proba=proba.svm, gene=gene.svm))

}

# ---------------------------------- #
# Logistic classification
# ---------------------------------- #
library(glmnet)
library(parallel)
predictorLogistic <- function(patientTrain, patientTest, labelTrain, positiveClass="1", alpha=1, intnfolds=5, measure="deviance") {
    # alpha=1 --> l1 penalty
    # alpha=0 --> l2 penalty
    # alpha=1/2 --> elastic net

    # Train a logistic model with internal cross validation
    cvfit = cv.glmnet(patientTrain, factor(labelTrain), family = "binomial", type.measure = measure, alpha=alpha, nfolds=intnfolds, parallel=TRUE)
    # Make predictions and output probabilities
    proba.log = predict(cvfit, newx = patientTest, s = "lambda.min", type="response")[,positiveClass]
    class.log = factor(predict(cvfit, newx = patientTest, s = "lambda.min", type="class")[,positiveClass])
    aa.log = predict(cvfit, newx = patientTest, s = "lambda.min", type="coef")[,1]
    gene.log = names(aa.log[aa.log>0][-1])

    return(list(class=class.log , proba=proba.log, gene=gene.log))
}
