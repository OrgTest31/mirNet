#' @title Do differential gene/microRNA expression analysis
#'
#' @description \code{pairComp} uses functions in \code{limma} package to easily compute the moderated t-statistics and p-values from differential gene/microRNA expression tests comparing between different phenotypes even when sample size is small.
#'
#' @param data A matrix, the normalized gene/microRNA expression dataset, should be a numeric matrix, with rows referring to genes/microRNAs and columns to samples.
#' @param class A vector of sample phenotypes. Sample phenotype in a scientific research could be treatment/control, normal/cancer or smoker/non-smoker. Different phenotypes could each be encoded as 0/1 when inputting to \code{LimmaFn}, for example, Normal-0; Cancer-1.
#'
#' @return A table with rows for all genes (ranked by significance) and columns of log2 fold-change, average expression, moderated t-statistic, p-value, adjusted p-value (default to Benjamini–Hochberg procedure). The table is the output of \code{\link[limma]{topTable}} function.
#'
#' @details This function computes the moderated t-statistic for users using empirical Bayes method, it is especially useful when the sample size is too small to perform parametric tests.
#'
#' Given a normalized gene expression or DNA methylation data matrix and a vector indicating sample phenotype, \code{LimmaFn} first fits a linear model using \code{\link[limma]{lmFit}}, then it refits the model and do comparisons between any two different phenotypes with \code{\link[limma]{contrasts.fit}}, finally it estimates moderated t-statistics for each comparison from the fitted model using empirical Bayes method (\code{\link[limma]{eBayes}}) and output the result from the \code{\link[limma]{topTable}} function.
#'
#' Note that doing the \code{contrasts.fit} step will not make a difference if you do comparison between two different sample status (treatment/control). However, When there are more than two sample status in your data set, this step will do comparison between every two status. And resulted summary tables will be stored in a list.
#'
#' @seealso \code{\link[limma]{lmFit}} for fitting a linear model, \code{\link[limma]{contrasts.fit}} for refitting, \code{\link[limma]{eBayes}} for Bayes method, \code{\link[limma]{topTable}} for the output table.
#'
#' @importFrom stats cor
#' @importFrom stats lm
#' @importFrom stats model.matrix
#' @importFrom stats pbinom
#' @importFrom stats pnorm
#' @importFrom stats sd
#'
#' @export pairComp
#'
#' @examples
#' # prepare your normalized data matrix
#' data.m <- matrix(rnorm(120), nrow = 20, ncol = 6)
#' # prepare the phenotype info ("C"-control; "T"-treatment)
#' class.v <- c('C', 'C', 'C', 'T', 'T', 'T')
#' # run function
#' lim.o <- pairComp(data = data.m, class = class.v)

pairComp <- function(data, class, padj = 'fdr'){
    
    # design matrix
    classFac.v = as.factor(classes)
    classLev.v <- levels(classFac.v)
    design.m <- model.matrix(~ 0 + classFac.v)
    colnames(design.m) <- classLev.v

    # pair-wise contrast matrix
    N <- length(classLev.v)
    compsNum <- 0.5*N*(N - 1)

    colN.v <- vector()
    for(i in 1:(N - 1)) colN.v <- c(colN.v, paste0(classLev.v[i], '-', classLev.v[(i + 1):N]))
    contr.m <- matrix(0, nrow = N, ncol = compsNum, dimnames = list(classLev.v, colN.v))

    for (i in 1:(N - 1)){

        clab <- classLev.v[(i + 1):N] 
        colN.v0 <- paste0(classLev.v[i], '-', clab)

        contr.m[classLev.v[i], colN.v0] <- 1
        for (j in clab) contr.m[j, colN.v0[clab == j]] <- -1
    }
    
    # linear model fit and empirical Bayesian estimation 
    # of adjusted variance and differential expression
    fit.o <- limma::contrasts.fit(limma::lmFit(data, design.m), contr.m)
    ebayes.o <- limma::eBayes(fit.o)

    # select top genes/microRNAs for every comparison
    res.l <- vector(length = compsNum, mode = 'list')
    names(res.l) <- colnames(contr.m)
    for(i in 1:compsNum){
        res.l[[i]] <- limma::topTable(ebayes.o, coef = i, adjust.method = padj, number = nrow(data))
    }
    list(res = res.l, contrast = contr.m)
}

