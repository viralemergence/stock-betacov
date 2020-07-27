

# effect of shuffling kernel matrices
shuffleGram <- function(K){
  indices <-  1:dim(K)[1]
  indices <-  sample(indices)
  return(K[indices,][,indices])
}

# effect of shuffling a trait
shufflecol <- function(X, col, makeGram=TRUE){
  indices <- 1:dim(X)[1]
  indices <- sample(indices)
  X[,col] <- X[indices,col]
  if(makeGram){
    D <- as.matrix(dist(X))^2
    return(exp(-D / median(D)))
  }else{
    return(X)
  }
}


n.perm <- 100
validation <- validate.hosts
lambda <- c(1, 0.1)
ref.perf <- validation(model.both)

n.traits <- dim(battraits.imp)[2]

# shuffling the matrices

imp.traits <- c()
imp.phylo <- c()
imp.both <- c()

for(rep in 1:n.perm){
  K.traits.perm <- 0 * K.traits + shuffleGram(K.traits)
  K.phylo.perm <- 0 * K.tree + shuffleGram(K.tree)
  K.both.perm <- 0 * K.both + shuffleGram(K.both)
  
  # traits
  model.perm <- tskrr(as.matrix(Yboth), 0.5 * (K.traits.perm + K.tree), G, lambda = lambda)
  imp.traits <- c(imp.traits, ref.perf - validation(model.perm))
  
  # phylo
  model.perm <- tskrr(as.matrix(Yboth), 0.5 * (K.traits + K.phylo.perm), G, lambda = lambda)
  imp.phylo <- c(imp.phylo, ref.perf - validation(model.perm))
  
  # both
  model.perm <- tskrr(as.matrix(Yboth), K.both.perm, G, lambda = lambda, testlabels = F)
  imp.both <- c(imp.both, ref.perf - validation(model.perm))
}

importance.scores <- data.frame(traits=imp.traits, phylogeny=imp.phylo, both=imp.both)

X <- battraits.imp[hosts.traits,]
#X[,"cites"] <- 0


for(i in 1:n.traits){
  imp.var <- c()
  for(rep in 1:n.perm){
    K.traits.perm <- shufflecol(X, i)
    model.perm <- tskrr(as.matrix(Yboth), 0.5 * (K.traits.perm + K.tree), G, lambda = lambda)
    imp.var <- c(imp.var, ref.perf - validation(model.perm))
  }
  importance.scores[colnames(X)[i]] <- imp.var
}


importance.scores.sorted <- importance.scores[,order(colMeans(-importance.scores))]

summary(importance.scores.sorted)

write.csv(importance.scores, "05_results/importance_scores.csv")

boxplot(importance.scores.sorted, las=2, col = c("cyan", "orange", "red", rep("green", dim(X)[2])), pch=19, ylab="drop in AUC",
        main="Variable importance\nby random permutation")

