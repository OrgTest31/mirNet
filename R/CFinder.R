#' @title Identify potential confounding factors or batch effects
#'
#' @description \code{CFinder} first does a principal component analysis via singular value decomposition, and compute correlation of a significant PCs (dimensionality estimated by \code{\link[isva]{EstDimRM}}) with sample class information.
#'
#' @param data A matrix, the normalized gene/microRNA expression dataset, should be a numeric matrix, with rows referring to genes/microRNAs and columns to samples.
#' @param pheno A vector of sample phenotypes. Sample phenotype in a scientific research could be categorical (treatment/control), or continuous (age). If you have multiple phenotypes, use a list with names denoting the corresponding phenotype.
#' @param type A string, 'categorical' if sample phenotype if categorical, 'continuous' if sample phenotype is continuous.
#'
#' @return A vector containing p-values for top k significant PCs.
#'
#' @details This function computes the moderated t-statistic for users using empirical Bayes method, it is especially useful when the sample size is too small to perform parametric tests.
#'
#' Given a normalized gene/microRNA expression data matrix and a vector indicating sample phenotype/class, \code{CFinder} first centers the input data matrix, then estimate the number of significant PCs using \code{\link[isva]{EstDimRM}}. Next, it does a singular value decomposition to the data matrix. Then it computes the correlation of significant PCs with potential confounding factors by using linear correlation (for continueous phenotype data like age) or Kruskal-Wallis rank sum test (for categorical phenotype data like tumor subtype). 
#'
#' @seealso \code{\link[isva]{EstDimRM}} for estimating the dimensionality of a dataset.
#'
#' @importFrom stats lm
#' @importFrom isva EstDimRMT
#' @importFrom stats sd
#'
#' @export CFinder
#'
#' @examples
#' # prepare your normalized data matrix
#' data.m <- matrix(rnorm(120), nrow = 20, ncol = 6)
#' # prepare the phenotype info ("C"-control; "T"-treatment)
#' class.v <- c('C', 'C', 'C', 'T', 'T', 'T')
#' # run function
#' p.v <- CFinder(data = data.m, pheno = class.v, type = 'categorical')

CFinder <- function(data, pheno, type){

    data.cent <- data - apply(data, 1, mean)
    svd.o <- svd(data.cent)

    # estimate data dimensionality
    k <- isva::EstDimRMT(data.cent)$dim
    
    # vector storing computed p-values
    p.v <- vector()

    if(type == 'continuous'){
        for(i in 1:k){
            lm.o <- lm(svd.o$v[,i] ~ pheno)
            p.v[i] <- summary(lm.o)$coeff[2, 4]
        }

    }

    else if(type == 'categorical'){
        for(i in 1:k){
            p.v[i] <- kruskal.test(svd.o$v[, i] ~ as.factor(pheno))$p.value
        }
    }

    else{
        stop('Wrong "type" value.')
    }

    p.v

}


