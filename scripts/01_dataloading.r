
# Loading data for Stock Models ####

library(missForest)
library(cluster)
library(ape)

conflict_prefer("as.matrix", "base")
conflict_prefer("dist", "stats")

setwd(paste0(here::here(), "/Github/Repos/stock-betacov"))

# INTERACTIONS
# ------------

read.csv(paste0(here::here(), '/Github/Repos/virionette/03_interaction_data/virionette.csv')) -> 
  
  interactions

interactions %<>% 
  mutate_at("host_species", ~.x %>% str_trim %>% str_replace_all(" ", "_"))

Y <- interactions %>% 
  dplyr::select(host_species, virus_genus) %>% 
  table() %>% 
  matrix(
    
    nrow = interactions %>% dplyr::select(host_species, virus_genus) %>% table() %>% nrow,
    ncol = interactions %>% dplyr::select(host_species, virus_genus) %>% table() %>% ncol
    
  )

dimnames(Y) <- interactions %>% 
  dplyr::select(host_species, virus_genus) %>% 
  table() %>% dimnames

Y %<>% as.data.frame()

# Y <- read.csv("03_interaction_data/Y.csv", row.names = 73)
hosts <- row.names(Y)
viruses <- colnames(Y)

n <- dim(Y)[1]
m <- dim(Y)[2]


# KERNEL VIRUSES
# ---------------

# compute genus sim
famgenus <- data.frame(genus=c())

# interactions <- read.csv("03_interaction_data/virionette.csv")

read.csv(paste0(here::here(), '/Github/Repos/virionette/03_interaction_data/virionette.csv')) -> 
  
  interactions

interactions %>% mutate_at("virus_family", 
                           ~factor(.x, levels = unique(c("", levels(.x))))) -> 
  
  interactions

virus.geni <- c() 

# for(v in viruses){
#   virus.geni <- c(virus.geni, interactions[interactions$virus_genus==v,"virus_family"][1])
# }

virus.geni <- viruses %>% 
  map(~interactions %>% filter(virus_genus == .x) %>% pull(virus_family) %>% unique %>% 
        ifelse(is.na(.), "", .)) %>% 
  map_chr(1) %>% 
  ifelse(. == "", 1, .) %>% as.numeric

# map_dbl(~interactions %>% filter(virus_genus == .x) %>% pull(virus_family) %>% unique)

G.fam <- matrix(0, nrow=m, ncol=m)

for(i in 1:m){
  for(j in 1:m){
    G.fam[i,j] <- virus.geni[i] == virus.geni[j]
  }
}

# smoother kernel
G <- G.fam +  0.1 * diag(m) + 0.1
rownames(G) <- colnames(Y)
colnames(G) <- colnames(Y)


# TRAITS HOSTS
# ------------

brtvar <- c('X5.1_AdultBodyMass_g',
            'X8.1_AdultForearmLen_mm',
            'X26.1_GR_Area_km2',
            'X26.2_GR_MaxLat_dd',
            'X26.3_GR_MinLat_dd',
            'X26.4_GR_MidRangeLat_dd',
            'X26.5_GR_MaxLong_dd',
            'X26.6_GR_MinLong_dd',
            'X26.7_GR_MidRangeLong_dd',
            'X27.2_HuPopDen_Mean_n.km2',
            'X27.4_HuPopDen_Change',
            'X28.1_Precip_Mean_mm',
            'X28.2_Temp_Mean_01degC',
            'X30.1_AET_Mean_mm',
            'Diet.Inv',
            'Diet.Fruit',
            'Diet.Nect',
            'Activity.Crepuscular',
            'ForStrat.Value_A',
            'ForStrat.Value_Ar',
            'ForStrat.Value_G',
            'ForStrat.Value_S', 'cites')

if(file.exists("04_predictors/battraits-completed.csv")){
  
  battraits.imp <- read.csv("04_predictors/battraits-completed.csv", row.names = 1)
  
}else{
  
  "04_predictors/Han-BatTraits.csv" %>% 
    read.csv(row.names = 1) ->
    # read.csv() -> 
    
    battraits1
  
  ManipulateColumns <- battraits1 %>% 
    dplyr::select(4:66) %>%
    names %>% setdiff(c("IUCN", "ForStrat.Value"))
  
  read.csv(paste0(here::here(), '/Github/Repos/virionette/04_predictors/Han-BatTraits.csv')) -> 
    
    battraits
  
  # UNCOMMENT IF YOU NEED TO battraits-completed.csv, takes 10'
  
  battraits$Activity.Nocturnal <- as.factor(battraits$Activity.Nocturnal)
  battraits$Activity.Diurnal <- as.factor(battraits$Activity.Diurnal)
  battraits$Activity.Crepuscular <- as.factor(battraits$Activity.Crepuscular)
  
  battraits.imputed <- missForest(battraits1[,ManipulateColumns], ntree=200)
  
  vartoremove <- c("X21.1_PopulationDensity_n", "X12.2_Terrestriality", "Diet.Vunk", "Diet.Scav")
  battraits.imp <- battraits.imputed$ximp[,-c(19, 23, 50, 51)]
  
  write.csv(battraits.imp, "04_predictors/battraits-completed.csv")
  
}

read.csv(paste0(here::here(), 
                '/Github/Repos/virionette/04_predictors/Citations.csv'), 
         row.names = 2) -> 
  
  citations

# citations <- read.csv("04_predictors/Citations.csv", row.names = 2)[row.names(battraits.imp),]

citations <- citations[row.names(battraits.imp),]

# subset of the variables
battraits.imp <- battraits.imp[,as.factor(brtvar)]
battraits.imp$cites <- citations$cites

battraits.imp <- scale(battraits.imp)

D.traits <- as.matrix(dist(battraits.imp))^2

# version with mean cites
battraits.imp.meancites <- as.data.frame(battraits.imp[,])
battraits.imp.meancites[,"cites"] <- mean(battraits.imp.meancites[,"cites"])


#D.traits <- as.matrix(daisy(battraits.imp, metric="gower", stand=T))
gamma <- median(D.traits)
K.traits.complete <- exp(-D.traits / gamma) + 0.1 * diag(dim(D.traits)[1]) + 0.1

hosts.traits <- intersect(hosts, row.names(battraits.imp))
Ytraits <- Y[hosts.traits,]
n.traits <- dim(Ytraits)[1]

K.traits <- K.traits.complete[hosts.traits, hosts.traits]

# now with mean traits

D.traits <- as.matrix(dist(battraits.imp.meancites))^2
K.traits.complete <- exp(-D.traits / gamma) + 0.1 * diag(dim(D.traits)[1]) + 0.1

K.traits.test <- K.traits.complete[,hosts.traits]


# PHYLOGENY HOSTS
# ---------------

readRDS(paste0(here::here(), 
               '/Github/Repos/virionette/04_predictors/bat-supertree_clean.rds')) -> 
  
  tree

# tree <- readRDS("04_predictors/bat-supertree_clean.rds")
tree <- as.phylo(tree)

D.tree <- cophenetic(tree)
K.tree.complete <- exp(-D.tree / median(D.tree)) + 0.1 + 0.1 * diag(dim(D.tree)[1])

hosts.tree <- intersect(hosts, row.names(D.tree))

Ytree <- Y[hosts.tree,]
n.tree <- dim(Ytree)[1]

K.tree <- K.tree.complete[hosts.tree, hosts.tree]

K.tree.test <- K.tree.complete[,hosts.tree]


# COMBINATION TRAITS AND PHYLO 
# ----------------------------

hosts.both <- intersect(hosts.traits, hosts.tree)  # host for which we have labels
hosts.both.test <- intersect(row.names(D.traits), row.names(D.tree))  # all the hosts
Yboth <- Y[hosts.both,]
n.both <- dim(Yboth)[1]

K.both <- 0.5 * (K.traits[hosts.both, hosts.both] + K.tree[hosts.both, hosts.both])

K.both.test <- 0.5 * (K.traits.complete[hosts.both.test,hosts.both] + K.tree.complete[hosts.both.test,hosts.both])

