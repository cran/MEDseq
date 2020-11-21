## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=7, fig.align = 'center', fig.show='hold',
                      dev.args=list(type="cairo"), warning=FALSE, message=FALSE, 
                      progress=FALSE, collapse=TRUE, comments="#>")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('devtools')
#  devtools::install_github('Keefe-Murphy/MEDseq')

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages('MEDseq')

## -----------------------------------------------------------------------------
library(MEDseq)

## ---- echo=FALSE--------------------------------------------------------------
suppressMessages(library(TraMineR))

## -----------------------------------------------------------------------------
data(mvad, package="MEDseq")
mvad$Location <- factor(apply(mvad[,5L:9L], 1L, function(x) which(x == "yes")), 
                        labels = colnames(mvad[,5L:9L]))
mvad          <- list(covariates = mvad[c(3L:4L,10L:14L,87L)], 
                      sequences = mvad[,15L:86L],
                      weights = mvad[,2L])
mvad.cov      <- mvad$covariates
mvad.seq      <- seqdef(mvad$sequences[-c(1L,2L)],
                        states = c("EM", "FE", "HE", "JL", "SC", "TR"),
                        labels = c("Employment", "Further Education", "Higher Education", 
                                   "Joblessness", "School", "Training"))

## ---- eval=FALSE--------------------------------------------------------------
#  mod1 <- MEDseq_fit(mvad.seq, G=11, modtype="UUN", weights=mvad$weights, gating=~ gcse5eq,
#                     covars=mvad.cov, control=MEDseq_control(noise.gate=FALSE))

## ---- eval=FALSE--------------------------------------------------------------
#  # 10-component CUN model with no covariates.
#  # CUN models have a precision parameter for each sequence position (i.e. time point),
#  # though each time point's precision is common across clusters.
#  
#  mod2 <- MEDseq_fit(mvad.seq, G=10, modtype="CUN", weights=mvad$weights)
#  
#  # 12-component CC model with all covariates.
#  # CC models have a single precision parameter across all clusters and time points.
#  
#  mod3 <- MEDseq_fit(mvad.seq, G=12, modtype="CC", weights=mvad$weights,
#                     gating=~. - Grammar - Location, covars=mvad.cov)

## ---- include=FALSE-----------------------------------------------------------
load(file="mvad_mod1.rda")
load(file="mvad_mod2.rda")
load(file="mvad_mod3.rda")

## -----------------------------------------------------------------------------
(comp <- MEDseq_compare(mod1, mod2, mod3, criterion="bic"))

## -----------------------------------------------------------------------------
opt   <- comp$optimal
summary(opt, classification = TRUE, parameters = FALSE, network = FALSE)

## -----------------------------------------------------------------------------
print(opt$gating)

## ---- eval=FALSE--------------------------------------------------------------
#  plot(opt, type="clusters")

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("MVAD_Clusters.png")

## ---- eval=FALSE--------------------------------------------------------------
#  plot(opt, type="mean")

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("MVAD_Mean.png")

## ---- fig.height=8.5----------------------------------------------------------
plot(opt, type="dbsvals")

## -----------------------------------------------------------------------------
MEDseq_meantime(opt, MAP=TRUE, norm=TRUE)

## -----------------------------------------------------------------------------
data(biofam, package="MEDseq")
biofam     <- list(covariates = biofam[2L:9L], 
                   sequences = biofam[10L:25L] + 1L)
biofam.cov <- biofam$covariates[,colSums(is.na(biofam$covariates)) == 0]
biofam.seq <- seqdef(biofam$sequences,
                     states = c("P", "L", "M", "L+M", 
                                "C", "L+C", "L+M+C", "D"),
                     labels = c("Parent", "Left", "Married", 
                                "Left+Marr", "Child", "Left+Child", 
                                "Left+Marr+Child", "Divorced"))

## ---- eval=FALSE--------------------------------------------------------------
#  # The UUN model includes a noise component.
#  # Otherwise, the model has a precision parameter for each time point in each cluster.
#  
#  bio <- MEDseq_fit(biofam.seq, G=10, modtype="UUN", gating=~ birthyr,
#                    covars=biofam.cov, noise.gate=FALSE)

## ---- echo=FALSE--------------------------------------------------------------
bio <- MEDseq_fit(biofam.seq, G=10, modtype="UUN", gating=~ birthyr, 
                  covars=biofam.cov, noise.gate=FALSE, verbose=FALSE)

## ---- echo=FALSE--------------------------------------------------------------
bio$call <- bio$call[-length(bio$call)]

## -----------------------------------------------------------------------------
print(bio)

## ---- eval=FALSE--------------------------------------------------------------
#  plot(bio, type="clusters", seriated="both")

## ---- echo=FALSE--------------------------------------------------------------
knitr::include_graphics("BIO_Clusters.png")

## -----------------------------------------------------------------------------
plot(bio, type="precision", quant.scale=TRUE, seriated="clusters")

## ---- fig.height=8.5----------------------------------------------------------
plot(bio, type="aswvals")

## -----------------------------------------------------------------------------
seqHtplot(biofam.seq)

## ---- fig.height=8.5----------------------------------------------------------
plot(bio, type="Ht")

