library(xnet)
library(pROC)

# VALIDATION METHODS
# ------------------

validate.interactions <- function(model){
  Yloo <- loo(model, exclusion="interaction")
  Y <- model@y
  nonzerorows <- colSums(Y) > 0
  return(auc(c(Y[,nonzerorows])>0, c(Yloo[,nonzerorows])))
}

validate.hosts <- function(model){
  Yloo <- loo(model, exclusion="row")
  Y <- model@y
  m <- dim(Y)[2]
  aucs <- c()
  viruses <- colnames(Y)
  for(i in 1:m){
    if(var(Y[,i]) > 0){
      a <- auc(Y[,i]>0, Yloo[,i])
      aucs <- c(aucs, a)
      if(TRUE){
        print(paste(viruses[i], ": auc =", a))
      }
    }
  }
  return(mean(aucs))
}

# TRAITS
# ------

model.traits <- tskrr(as.matrix(Ytraits), K.traits, G, lambda = c(0.1, 1))
(auc.traits.interactions <- validate.interactions(model.traits))
(auc.traits.hosts <- validate.hosts(model.traits))

# PHYLOGENY
# --------

model.phylo <- tskrr(as.matrix(Ytree), K.tree, G, lambda = c(1, 0.1))
(auc.phylo.interactions <- validate.interactions(model.phylo))
(auc.phylo.hosts <- validate.hosts(model.phylo))


# COMBINATION TRAITS AND PHYLO 
# ----------------------------

model.both <- tskrr(as.matrix(Yboth), K.both, G, lambda = c(1, 0.1))
(auc.both.interactions <- validate.interactions(model.both))
(auc.both.hosts <- validate.hosts(model.both))

predictions.both <- predict(model.both, k=K.both.test)
write.csv(predictions.both, "05_results/scores_tskrr_both.csv")


# WITHOUT CITATIONS
# -----------------

model.both.nocites <- tskrr(as.matrix(Yboth), 0.5 * (K.traits.complete[hosts.both, hosts.both] + K.tree), G, lambda = c(1, 0.1))
(auc.both.interactions <- validate.interactions(model.both.nocites))
(auc.both.hosts <- validate.hosts(model.both.nocites))

predictions.both.nocites <- predict(model.both, k=K.traits.complete[, hosts.both])
write.csv(predictions.both.nocites, "05_results/scores_tskrr_both_nocites.csv")
