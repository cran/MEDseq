% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{MEDseq_fit}
\alias{MEDseq_fit}
\alias{summary.MEDseq}
\alias{print.MEDseq}
\title{MEDseq: Mixtures of Exponential-Distance Models with Covariates}
\usage{
MEDseq_fit(seqs, 
           G = 1L:9L, 
           modtype = c("CC", "UC", "CU", "UU", 
                       "CCN", "UCN", "CUN", "UUN"), 
           gating = NULL, 
           weights = NULL, 
           ctrl = MEDseq_control(...), 
           covars = NULL, 
           ...)

\method{summary}{MEDseq}(object,
        classification = TRUE,
        parameters = FALSE,
        network = FALSE,
        SPS = FALSE,
        ...)

\method{print}{MEDseq}(x,
      digits = 3L,
      ...)
}
\arguments{
\item{seqs}{A state-sequence object of class \code{"stslist"} as created by the \code{\link[TraMineR]{seqdef}} function in the \pkg{TraMineR} package (which is reexported by \pkg{MEDseq} for convenience). Note that the data set must have equal sequence lengths, the intervals are assumed to be evenly spaced, and missingness is not allowed.}

\item{G}{A positive integer vector specifying the numbers of mixture components (clusters) to fit. Defaults to \code{G=1:9}.}

\item{modtype}{A vector of character strings indicating the type of MEDseq models to be fitted, in terms of the constraints or lack thereof on the precision parameters. By default, all valid model types are fitted (except some only where \code{G > 1} or \code{G > 2}, see \code{Note}). 
The models are named \code{"CC"}, \code{"CU"}, \code{"UC"}, \code{"UU"}, \code{"CCN"}, \code{"CUN"}, \code{"UCN"}, and \code{"UUN"}. The first letter denotes whether the precision parameters are constrained/unconstrained across clusters. The second letter denotes whether the precision parameters are constrained/unconstrained across sequence positions (i.e. time points). The third letter denotes whether one of the components is constrained to have zero-precision/infinite variance. Such a noise component assumes sequences in that cluster follow a uniform distribution.}

\item{gating}{A \code{\link[stats]{formula}} for determining the model matrix for the multinomial logistic regression in the gating network when fixed covariates enter the mixing proportions. Defaults to \code{~1}, i.e. no covariates. This will be ignored where \code{G=1}. Continuous, categorical, and/or ordinal covariates are allowed. Logical covariates will be coerced to factors. Interactions, transformations, and higher order terms are permitted: the latter \strong{must} be specified explicitly using the \code{AsIs} operator (\code{\link{I}}). The specification of the LHS of the formula is ignored. Intercept terms are included by default.}

\item{weights}{Optional numeric vector containing observation-specific sampling weights, which are accounted for in the model fitting and other functions where applicable. \code{weights} are always internally normalised to sum to the sample size. See the \code{unique} argument to \code{\link{MEDseq_control}} to see how incorporating weights also yields computational benefits. Note that \code{weights} must \strong{always} be explicitly supplied here; it is not enough to use weights when constructing the state sequence object via \code{\link[TraMineR]{seqdef}} (reexported by \pkg{MEDseq} for convenience). If you \emph{are} using a weighted \code{"stslist"} state sequence object and do not specify \code{weights}, you will be prompted to explicitly specify \code{weights=attr(seqs, "weights")} for a weighted model or \code{weights=NULL} for an unweighted model.}

\item{ctrl}{A list of control parameters for the EM/CEM and other aspects of the algorithm. The defaults are set by a call to \code{\link{MEDseq_control}}.}

\item{covars}{An optional data frame (or a matrix with named columns) in which to look for the covariates in the \code{gating} network formula, if any. If not found in \code{covars}, any supplied \code{gating} covariates are taken from the environment from which \code{MEDseq_fit} is called. Try to ensure the names of variables in \code{covars} do not match any of those in \code{seqs}.}

\item{...}{Catches unused arguments (see \code{\link{MEDseq_control}}).}

\item{x, object, digits, classification, parameters, network, SPS}{Arguments required for the \code{print} and \code{summary} functions: \code{x} and \code{object} are objects of class \code{"MEDseq"} resulting from a call to \code{\link{MEDseq_fit}}, while \code{digits} gives the number of decimal places to round to for printing purposes (defaults to \code{3}). \code{classification}, \code{parameters}, and \code{network} are logicals which govern whether a table of the MAP classification of observations, the mixture component parameters, and the gating network coefficients are printed, respectively. \code{SPS} governs the printing of the relevant quantities in \code{"summaryMEDseq"} objects when any of \code{classification}, \code{parameters}, &/or \code{network} are \code{TRUE} (see \code{\link{MEDseq_clustnames}} and \code{\link[TraMineR]{seqformat}}).}
}
\value{
A list (of class \code{"MEDseq"}) with the following named entries (of which some may be missing, depending on the \code{criterion} employed), mostly corresponding to the chosen optimal model (as determined by the \code{criterion} within \code{\link{MEDseq_control}}):
\item{\code{call}}{The matched call.}
\item{\code{data}}{The input data, \code{seqs}.}
\item{\code{modtype}}{A character string denoting the MEDseq model type at which the optimal \code{criterion} occurs.}
\item{\code{G}}{The optimal number of mixture components according to \code{criterion}.}
\item{\code{params}}{A list with the following named components:
\describe{
\item{\code{theta}}{A matrix with \code{G} rows and T columns, where T is the number of sequence positions, giving the central sequences of each cluster. The mean of the noise component is not reported, as it does not contribute in any way to the likelihood. A dedicated \code{print} function is provided.}
\item{\code{lambda}}{A matrix of precision parameters. Will contain \code{1} row if the 1st letter of \code{modtype} is "C" and \code{G} columns otherwise. Will contain \code{1} column if the 2nd letter of \code{modtype} is "C" and T columns otherwise, where T is the number of sequence positions. Precision parameter values of zero are reported for the noise component, if any. Note that values of \code{Inf} are also possible, corresponding to zero-variance, which is most likely under the \code{"UU"} or \code{"UUN"} models. A dedicated \code{print} function is provided.}
\item{\code{tau}}{The mixing proportions: either a vector of length \code{G} or, if \code{gating} covariates were supplied, a matrix with an entry for each observation (rows) and component (columns).}}
}
\item{\code{gating}}{An object of class \code{"MEDgating"} (for which dedicated \code{print}, \code{summary}, and \code{\link[=predict.MEDgating]{predict}} methods exist) and either \code{"multinom"} or \code{"glm"} (only for single-component models) giving the \code{\link[nnet]{multinom}} regression coefficients of the \code{gating} network. If \code{gating} covariates were \emph{NOT} supplied (or the best model has just one component), this corresponds to a RHS of \code{~1}, otherwise the supplied \code{gating} formula. As such, a fitted \code{gating} network is always returned even in the absence of supplied covariates or clusters. If there is a noise component (and the option \code{noise.gate=TRUE} is invoked), its coefficients are those for the \emph{last} component. \strong{Users are cautioned against making inferences about statistical significance from summaries of the coefficients in the gating network. Users are instead advised to use the function \code{\link{MEDseq_stderr}}}.}
\item{\code{z}}{The final responsibility matrix whose \code{[i,k]}-th entry is the probability that observation \emph{i} belongs to the \emph{k}-th component. If there is a noise component, its values are found in the \emph{last} column.}
\item{\code{MAP}}{The vector of cluster labels for the chosen model corresponding to \code{z}, i.e. \code{max.col(z)}. Observations belonging to the noise component, if any, will belong to component \code{0}.}
\item{\code{BIC}}{A matrix of \emph{all} BIC values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
\item{\code{ICL}}{A matrix of \emph{all} ICL values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
\item{\code{AIC}}{A matrix of \emph{all} AIC values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
\item{\code{DBS}}{A matrix of \emph{all} (weighted) mean/median DBS values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note} and \code{\link{dbs}}.}
\item{\code{DBSvals}}{A list of lists giving the observation-specific DBS values for \emph{all} fitted models. The first level of the list corresponds to numbers of components, the second to the MEDseq model types.}
\item{\code{dbs}}{The (weighted) mean/median DBS value corresponding to the optimal model. May not necessarily be the optimal DBS.}
\item{\code{dbsvals}}{Observation-specific DBS values corresponding to the optimum model, which may not be optimal in terms of DBS.}
\item{\code{ASW}}{A matrix of \emph{all} (weighted) mean/median ASW values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
\item{\code{ASWvals}}{A list of lists giving the observation-specific ASW values for \emph{all} fitted models. The first level of the list corresponds to numbers of components, the second to the MEDseq model types.}
\item{\code{asw}}{The (weighted) mean/median ASW value corresponding to the optimal model. May not necessarily be the optimal ASW.}
\item{\code{aswvals}}{Observation-specific ASW values corresponding to the optimum model, which may not be optimal in terms of ASW.}
\item{\code{LOGLIK}}{A matrix of \emph{all} maximal log-likelihood values with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
\item{\code{DF}}{A matrix giving the numbers of estimated parameters (i.e. the number of 'used' degrees of freedom) for \emph{all} visited models, with \code{length{G}} rows and \code{length(modtype)} columns. Subtract these numbers from the sample size to get the degrees of freedom. See \code{Note}.}
\item{\code{ITERS}}{A matrix giving the total number of EM/CEM iterations for \emph{all} visited models, with \code{length{G}} rows and \code{length(modtype)} columns. See \code{Note}.}
\item{\code{CV}}{A matrix of \emph{all} cross-validated log-likelihood values with \code{length{G}} rows and \code{length(modtype)} columns, if available. See \code{Note} and the arguments \code{do.cv} and \code{nfolds} to \code{\link{MEDseq_control}}.}
\item{\code{NEC}}{A matrix of \emph{all} NEC values with \code{length{G}} rows and \code{length(modtype)} columns, if available. See \code{Note} and the argument \code{do.nec} to \code{\link{MEDseq_control}}.}
\item{\code{bic}}{The BIC value corresponding to the optimal model. May not necessarily be the optimal BIC.}
\item{\code{icl}}{The ICL value corresponding to the optimal model. May not necessarily be the optimal ICL.}
\item{\code{aic}}{The AIC value corresponding to the optimal model. May not necessarily be the optimal AIC.}
\item{\code{loglik}}{The vector of increasing log-likelihood values for every EM/CEM iteration under the optimal model. The last element of this vector is the maximum log-likelihood achieved by the parameters returned at convergence.}
\item{\code{df}}{The number of estimated parameters in the optimal model (i.e. the number of 'used' degrees of freedom). Subtract this number from the sample size to get the degrees of freedom.}
\item{\code{iters}}{The total number of EM/CEM iterations for the optimal model.}
\item{\code{cv}}{The cross-validated log-likelihood value corresponding to the optimal model, if available. May not necessarily be the optimal one.}
\item{\code{nec}}{The NEC value corresponding to the optimal model, if available. May not necessarily be the optimal NEC.}
\item{\code{ZS}}{A list of lists giving the \code{z} matrices for \emph{all} fitted models. The first level of the list corresponds to numbers of components, the second to the MEDseq model types.}
\item{\code{uncert}}{The uncertainty associated with the \code{classification}.}
\item{\code{covars}}{A data frame gathering the set of covariates used in the \code{gating} network, if any. Will contain zero columns in the absence of gating covariates. Supplied gating covariates will be excluded if the optimal model has only one component. May have fewer columns than covariates supplied via the \code{covars} argument also, as only the included covariates are gathered here.}
Dedicated \code{\link[=plot.MEDseq]{plot}}, \code{print}, and \code{summary} functions exist for objects of class \code{"MEDseq"}.
}
\description{
Fits MEDseq models: mixtures of Exponential-Distance models with gating covariates and sampling weights. Typically used for clustering categorical/longitudinal life-course sequences. Additional arguments are available via the function \code{\link{MEDseq_control}}.
}
\details{
The function effectively allows 8 different MEDseq precision parameter settings for models with or without gating network covariates. By constraining the mixing proportions to be equal (see \code{equalPro} in \code{\link{MEDseq_control}}) an extra special case is facilitated in the latter case. 

While model selection in terms of choosing the optimal number of components and the MEDseq model type is performed within \code{\link{MEDseq_fit}}, using one of the \code{criterion} options within \code{\link{MEDseq_control}}, choosing between multiple fits with different combinations of covariates or different initialisation settings can be done by supplying objects of class \code{"MEDseq"} to \code{\link{MEDseq_compare}}.
}
\note{
Where \code{BIC}, \code{ICL}, \code{AIC}, \code{DBS}, \code{ASW}, \code{LOGLIK}, \code{DF}, \code{ITERS}, \code{CV}, and \code{NEC} contain \code{NA} entries, this corresponds to a model which was not run; for instance a UU model is never run for single-component models as it is equivalent to CU, while a UCN model is never run for two-component models as it is equivalent to CCN. As such, one can consider the value as not really missing, but equivalent to the corresponding value. On the other hand, \code{-Inf} represents models which were terminated due to error, for which a log-likelihood could not be estimated. These objects all inherit the class \code{"MEDCriterion"} for which dedicated \code{print} and \code{summary} methods exist. For plotting, please see \code{\link[=plot.MEDseq]{plot}}.
}
\examples{
\dontshow{if (interactive()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
# Load the MVAD data
data(mvad)
mvad$Location <- factor(apply(mvad[,5:9], 1L, function(x) 
                 which(x == "yes")), labels = colnames(mvad[,5:9]))
mvad          <- list(covariates = mvad[c(3:4,10:14,87)],
                      sequences = mvad[,15:86], 
                      weights = mvad[,2])
mvad.cov      <- mvad$covariates

# Create a state sequence object with the first two (summer) time points removed
states        <- c("EM", "FE", "HE", "JL", "SC", "TR")
labels        <- c("Employment", "Further Education", "Higher Education", 
                   "Joblessness", "School", "Training")
mvad.seq      <- seqdef(mvad$sequences[-c(1,2)], states=states, labels=labels)

# Fit a range of exponential-distance models without clustering
mod0          <- MEDseq_fit(mvad.seq, G=1)

\donttest{
# Fit a range of unweighted mixture models without covariates
# Only consider models with a noise component
# Supply some MEDseq_control() arguments
# mod1        <- MEDseq_fit(mvad.seq, G=9:10, modtype=c("CCN", "CUN", "UCN", "UUN"),
#                           algo="CEM", init.z="kmodes", criterion="icl")

# Fit a model with weights and a gating covariate
# Have the probability of noise-component membership be constant
mod2          <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, 
                            gating=~ gcse5eq, covars=mvad.cov, noise.gate=FALSE)
                            
# Examine this model in greater detail
summary(mod2, classification=TRUE, parameters=TRUE)
summary(mod2$gating, SPS=TRUE)
print(mod2$params$theta, SPS=TRUE)
plot(mod2, "clusters")}
\dontshow{\}) # examplesIf}
}
\references{
Murphy, K., Murphy, T. B., Piccarreta, R., and Gormley, I. C. (2021). Clustering longitudinal life-course sequences using mixtures of exponential-distance models. \emph{Journal of the Royal Statistical Society: Series A (Statistics in Society)}, 184(4): 1414-1451. <\doi{10.1111/rssa.12712}>.
}
\seealso{
\code{\link[TraMineR]{seqdef}} (reexported by \pkg{MEDseq} for convenience), \code{\link{MEDseq_control}}, \code{\link{MEDseq_compare}}, \code{\link{plot.MEDseq}}, \code{\link{predict.MEDgating}}, \code{\link{MEDseq_stderr}}, \code{\link{I}}, \code{\link{MEDseq_clustnames}}, \code{\link[TraMineR]{seqformat}}
}
\author{
Keefe Murphy - <\email{keefe.murphy@mu.ie}>
}
\keyword{clustering}
\keyword{main}
