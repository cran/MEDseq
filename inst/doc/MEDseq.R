## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(fig.width=7, fig.height=7, fig.align = 'center', fig.show='hold',
                      dev.args=list(type="cairo"), warning=FALSE, message=FALSE, 
                      progress=FALSE, collapse=TRUE, comments="#>")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages('devtools')
#  devtools::install_github('Keefe-Murphy/MEDseq')

## ---- eval=FALSE---------------------------------------------------------
#  install.packages('MEDseq')

## ------------------------------------------------------------------------
library(MEDseq)

## ---- echo=FALSE---------------------------------------------------------
suppressMessages(library(TraMineR))

## ------------------------------------------------------------------------
data(mvad, package="MEDseq")
mvad$Location <- factor(apply(mvad[,5L:9L], 1L, function(x) which(x == "yes")), 
                        labels = colnames(mvad[,5L:9L]))
mvad          <- list(covariates = mvad[c(3L:4L,10L:14L,87L)], 
                      sequences = mvad[,15L:86L],
                      weights = mvad[,2L])
mvad.cov      <- mvad$covariates
mvad.seq      <- seqdef(mvad$sequences[,-1L],
                        states = c("EM", "FE", "HE", "JL", "SC", "TR"),
                        labels = c("Employment", "Further Education", "Higher Education", 
                                   "Joblessness", "School", "Training"))

## ---- eval=FALSE---------------------------------------------------------
#  mod1 <- MEDseq_fit(mvad.seq, G=10, modtype="UCN", weights=mvad$weights,
#                     gating=~ fmpr + gcse5eq + livboth, covars=mvad.cov)

## ---- echo=FALSE---------------------------------------------------------
mod1 <- MEDseq_fit(mvad.seq, G=10, modtype="UCN", weights=mvad$weights, 
                   gating=~ fmpr + gcse5eq + livboth, covars=mvad.cov, verbose=FALSE)

## ------------------------------------------------------------------------
# 9-component CUN model with no covariates.
# The CUN model has a precision parameter for each sequence position,
# though each time-point's precision is common across clusters.

mod2 <- MEDseq_fit(mvad.seq, G=9,  modtype="CUN", weights=mvad$weights, verbose=FALSE)

# 11-component CC model with all coviarates.
# The CC model has only a single precision parameter across all clusters and time-points.

mod3 <- MEDseq_fit(mvad.seq, G=11, modtype="CC", weights=mvad$weights,
                   gating=~. - Grammar - Location, covars=mvad.cov, verbose=FALSE)

## ---- eval=FALSE---------------------------------------------------------
#  (comp <- MEDseq_compare(mod1, mod2, mod3, criterion="dbs"))
#  opt   <- comp$optimal
#  summary(opt)

## ---- echo=FALSE---------------------------------------------------------
(comp <- MEDseq_compare(mod1, mod2, mod3, criterion="dbs"))
opt   <- comp$optimal
suppressMessages(summary(opt))

## ------------------------------------------------------------------------
coef(opt$gating)

## ------------------------------------------------------------------------
plot(opt, type="clusters")

## ------------------------------------------------------------------------
plot(opt, type="mean")

## ---- fig.height=8.5-----------------------------------------------------
plot(opt, type="dbsvals")

## ------------------------------------------------------------------------
MEDseq_meantime(opt, MAP=TRUE, norm=TRUE)

## ------------------------------------------------------------------------
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

## ---- eval=FALSE---------------------------------------------------------
#  # The UUN model includes a noise component.
#  # Otherwise, the model has a precision parameter for each sequence position in each cluster.
#  
#  bio <- MEDseq_fit(biofam.seq, G=10, modtype="UUN", gating=~ birthyr,
#                    covars=biofam.cov, control=MEDseq_control(noise.gate=FALSE))

## ---- echo=FALSE---------------------------------------------------------
bio <- MEDseq_fit(biofam.seq, G=10, modtype="UUN", gating=~ birthyr, 
                  covars=biofam.cov, control=MEDseq_control(noise.gate=FALSE), verbose=FALSE)

## ---- echo=FALSE---------------------------------------------------------
bio$call <- bio$call[-length(bio$call)]
print(bio)

## ------------------------------------------------------------------------
plot(bio, type="clusters", seriate="both")

## ------------------------------------------------------------------------
plot(bio, type="precision", log.scale=TRUE)

## ---- fig.height=8.5-----------------------------------------------------
plot(bio, type="aswvals")

## ------------------------------------------------------------------------
seqHtplot(biofam.seq)

## ---- fig.height=8.5-----------------------------------------------------
plot(bio, type="Ht")

