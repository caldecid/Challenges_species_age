
# Testing the function wit different extinction and missing species -------

##sourcing the libraries and the directories
source(file.path(getwd(), "/source.R"))


# no extinction scenario -------------------------------------------------

## 25% missing species
trees.25.missing <- sim.bd.taxa(134, numbsim = 1000,
                                lambda = 0.3, mu = 0)

## 50% missing species
trees.50.missing <- sim.bd.taxa(200, numbsim = 1000,
                                lambda = 0.3, mu = 0)

## 25% missing species
trees.25.extant <- lapply(trees.25.missing, prune.fossil.tips)

##root ages
root.25 <- sapply(trees.25.missing, tree.max)

### 50% missing species
trees.50.extant <- lapply(trees.50.missing, prune.fossil.tips)

## root ages
root.50 <- sapply(trees.50.missing, tree.max)

# bifurcating  ------------------------------------------------------------

## 25% missing species
tax.25.bif <- lapply(trees.25.missing, sim.taxonomy, beta = 1)

##get ages extant species
ages.25.bif <- Map(getAgesExtantSpecies, tax.25.bif, trees.25.missing,
                   Tol = 1e-6)

## 50% missing species
tax.50.bif <- lapply(trees.50.missing, sim.taxonomy, beta = 1)

## get ages extant species
ages.50.bif <- Map(getAgesExtantSpecies, tax.50.bif, trees.50.missing,
                   Tol = 1e-6)

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bif, trees.25.missing,
                          SamplingFrac = 0.25)


ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.25$root.age <- rep(root.25, each = 100)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.25 <- ages.bif.incomp.25 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.25$extinction <- rep("no_ext", nrow(ages.bif.incomp.25))

##species name
ages.bif.incomp.25$species <- paste0(ages.bif.incomp.25$label,".",
                                     ages.bif.incomp.25$tree)

##sampling fraction
ages.bif.incomp.25$fraction <- rep(0.25, nrow(ages.bif.incomp.25))


###speciation mode
ages.bif.incomp.25$speciation <- rep("bif", nrow(ages.bif.incomp.25))

#####sampling fraction = 0.50##########

##get ages extant species with an uniform incomple sampling (0.50)
ages.bif.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bif, trees.50.missing,
                          SamplingFrac = 0.50)


ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.50$root.age <- rep(root.50, each = 100)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.50 <- ages.bif.incomp.50 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.50$extinction <- rep("no_ext", nrow(ages.bif.incomp.50))

##species name
ages.bif.incomp.50$species <- paste0(ages.bif.incomp.50$label,".",
                                     ages.bif.incomp.50$tree)

##sampling fraction
ages.bif.incomp.50$fraction <- rep(0.50, nrow(ages.bif.incomp.50))

###speciation mode
ages.bif.incomp.50$speciation <- rep("bif", nrow(ages.bif.incomp.50))

##binding dataframes
ages_incomplete_no <- rbind(ages.bif.incomp.25, ages.bif.incomp.50)

##saving
write_csv(ages_incomplete_no, file = 
            "results/data/processed/incomplete_sampling/ages.incomplete.no.ex.csv")


# Probabilistic function --------------------------------------------------

###for 25% incomplete sampling
ages.incomplete_no.f25 <- ages_incomplete_no %>% filter(fraction == "0.25")

##correcting using the exact rho

f25_no.list <- list()

for(i in 1:nrow(ages.incomplete_no.f25)){
  f25_no.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                       mu = 0,
                                       node_age = ages.incomplete_no.f25$Incomp.age[i],
                                       rho = 0.75)
}

##mean age
f25_no.mean.age <- sapply(f25_no.list, function(x) sum(x$time * x$prob))

##median age
f25_no.median.age <- sapply(f25_no.list, function(x)
                              weighted.median(x = x$time,
                                              w = x$prob))




##binding columns
ages.f25_no <- cbind(ages.incomplete_no.f25,
                     f25_no.mean.age, f25_no.median.age) %>% 
  rename(mean.age = f25_no.mean.age,
         median.age = f25_no.median.age)



# 50% missing species -----------------------------------------------------

ages.incomplete_no.f50 <- ages_incomplete_no %>% filter(fraction == "0.5")

##correcting using the exact rho

f50_no.list <- list()

for(i in 1:nrow(ages.incomplete_no.f50)){
  f50_no.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                          mu = 0,
                                          node_age = ages.incomplete_no.f50$Incomp.age[i],
                                          rho = 0.50)
}

##mean age
f50_no.mean.age <- sapply(f50_no.list, function(x) sum(x$time * x$prob))

##median age
f50_no.median.age <- sapply(f50_no.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))



##binding columns
ages.f50_no <- cbind(ages.incomplete_no.f50,
                     f50_no.mean.age, f50_no.median.age) %>% 
  rename(mean.age = f50_no.mean.age,
         median.age = f50_no.median.age)


ages.random_no_ext <- rbind(ages.f25_no, ages.f50_no)

write_csv(ages.random_no_ext,
      file = "results/data/processed/incomplete_sampling/ages.random_no_ext.csv")

ages.random_no_ext <- read_csv("results/data/processed/incomplete_sampling/ages.random_no_ext.csv")

###calculating MAPE

ages.mape.no_ext <- ages.random_no_ext %>% group_by(fraction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Incomp.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)


ages.mape.no_ext$extinction <- "no_ext"



##binding mape dataframes
ages.mape.sampling_no <- ages.mape.no_ext %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.mape.sampling_no$estimate <- factor(ages.mape.sampling_no$estimate,
                                          levels = c("mape.incomplete",
                                                     "mape.mean",
                                                     "mape.median"),
                                          ordered = TRUE)

# Intermediate extinction -------------------------------------------------

# Testing the new function with incomplete sampling ---------------------------
ages.incomplete <- read_csv("results/data/processed/incomplete_sampling/ages.incomplete_int.csv")

##Filtering bifurcating speciation simulation
ages.incomplete.bif <- ages.incomplete %>% filter(speciation == "bif")

ages.incomplete.bif$fraction <- as.factor(ages.incomplete.bif$fraction)

###for 25% incomplete sampling
ages.incomplete.25 <- ages.incomplete.bif %>% filter(fraction == "0.25")

##correcting using the exact rho

f25.list <- list()

for(i in 1:nrow(ages.incomplete.25)){
  f25.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                       mu = 0.15,
                                       node_age = ages.incomplete.25$Incomp.age[i],
                                       rho = 0.75)
}

##mean age
f25.mean.age <- sapply(f25.list, function(x) sum(x$time * x$prob))

##median age
f25.median.age <- sapply(f25.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))

##binding columns
ages.f25 <- cbind(ages.incomplete.25, f25.mean.age, f25.median.age) %>% 
  rename(mean.age = f25.mean.age,
         median.age = f25.median.age)



###for 50% incomplete sampling
ages.incomplete.50 <- ages.incomplete.bif %>% filter(fraction == "0.5")

##correcting
f50.list <- list()

for(i in 1:nrow(ages.incomplete.50)){
  f50.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                       mu = 0.15,
                                       node_age = ages.incomplete.50$Incomp.age[i],
                                       rho = 0.50)
}

##mean age
f50.mean.age <- sapply(f50.list, function(x) sum(x$time * x$prob))

##median age
f50.median.age <- sapply(f50.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))


##binding columns
ages.f50 <- cbind(ages.incomplete.50, f50.mean.age, f50.median.age) %>% 
  rename(mean.age = f50.mean.age,
         median.age = f50.median.age)


##binding dataframes
ages.random_int_ext <- rbind(ages.f25, ages.f50) 

ages.random_int_ext$extinction

##saving
write_csv(ages.random_int_ext,
         file = "results/data/processed/incomplete_sampling/ages.random_int_ext.csv")

ages.random_int_ext <- read_csv(file = "results/data/processed/incomplete_sampling/ages.random_int_ext.csv")

ages.mape.int <- ages.random_int_ext %>% group_by(fraction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Incomp.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)

ages.mape.int$extinction <- "intermediate"

##pivot longer
ages.mape.sampling_int <- ages.mape.int %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.mape.sampling_int$estimate <- factor(ages.mape.sampling_int$estimate,
                                           levels = c("mape.incomplete",
                                                      "mape.mean",
                                                      "mape.median"),
                                           ordered = TRUE)



# High extinction ---------------------------------------------------------

## 25% missing species
trees.25.missing <- sim.bd.taxa(134, numbsim = 1000,
                                lambda = 0.3, mu = 0.25)

## 50% missing species
trees.50.missing <- sim.bd.taxa(200, numbsim = 1000,
                                lambda = 0.3, mu = 0.25)

## 25% missing species
trees.25.extant <- lapply(trees.25.missing, prune.fossil.tips)

##root ages
root.25 <- sapply(trees.25.missing, tree.max)

### 50% missing species
trees.50.extant <- lapply(trees.50.missing, prune.fossil.tips)

## root ages
root.50 <- sapply(trees.50.missing, tree.max)

# bifurcating  ------------------------------------------------------------

## 25% missing species
tax.25.bif <- lapply(trees.25.missing, sim.taxonomy, beta = 1)

##get ages extant species
ages.25.bif <- Map(getAgesExtantSpecies, tax.25.bif, trees.25.missing,
                   Tol = 1e-6)

## 50% missing species
tax.50.bif <- lapply(trees.50.missing, sim.taxonomy, beta = 1)

## get ages extant species
ages.50.bif <- Map(getAgesExtantSpecies, tax.50.bif, trees.50.missing,
                   Tol = 1e-6)

###########uniform incomple sampling of extant species##########################

##get ages extant species with an uniform incomple sampling (0.25)
ages.bif.incomp.25 <- Map(getAgesExtantIncompleteSampling, 
                          ages.25.bif, trees.25.missing,
                          SamplingFrac = 0.25)


ages.bif.incomp.25 <- do.call("rbind", ages.bif.incomp.25) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.25$root.age <- rep(root.25, each = 100)

##tree numbers
ages.bif.incomp.25$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.25 <- ages.bif.incomp.25 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.25$extinction <- rep("high", nrow(ages.bif.incomp.25))

##species name
ages.bif.incomp.25$species <- paste0(ages.bif.incomp.25$label,".",
                                     ages.bif.incomp.25$tree)

##sampling fraction
ages.bif.incomp.25$fraction <- rep(0.25, nrow(ages.bif.incomp.25))


###speciation mode
ages.bif.incomp.25$speciation <- rep("bif", nrow(ages.bif.incomp.25))

#####sampling fraction = 0.50##########

##get ages extant species with an uniform incomple sampling (0.50)
ages.bif.incomp.50 <- Map(getAgesExtantIncompleteSampling, 
                          ages.50.bif, trees.50.missing,
                          SamplingFrac = 0.50)


ages.bif.incomp.50 <- do.call("rbind", ages.bif.incomp.50) %>% 
  rename(True.age = Age,
         Estimated.age = tip_length_reconstr,
         Incomp.age = IncompletePhyloAge) %>%
  select(label, True.age, Estimated.age, Incomp.age)

##root age
ages.bif.incomp.50$root.age <- rep(root.50, each = 100)

##tree numbers
ages.bif.incomp.50$tree <- rep(paste0("tree.", 1:1000), each = 100)


##scale ages
ages.bif.incomp.50 <- ages.bif.incomp.50 %>%
  mutate(rTrue.age = True.age/root.age,
         rPhylo.age = Estimated.age/root.age,
         rIncomp.age = Incomp.age/root.age)

##extinction
ages.bif.incomp.50$extinction <- rep("high", nrow(ages.bif.incomp.50))

##species name
ages.bif.incomp.50$species <- paste0(ages.bif.incomp.50$label,".",
                                     ages.bif.incomp.50$tree)

##sampling fraction
ages.bif.incomp.50$fraction <- rep(0.50, nrow(ages.bif.incomp.50))

###speciation mode
ages.bif.incomp.50$speciation <- rep("bif", nrow(ages.bif.incomp.50))

##binding dataframes
ages_incomplete_high <- rbind(ages.bif.incomp.25, ages.bif.incomp.50)

##saving
write_csv(ages_incomplete_high, file = 
            "results/data/processed/incomplete_sampling/ages.incomplete.high.ex.csv")


# Probabilistic function --------------------------------------------------

###for 25% incomplete sampling
ages.incomplete_high.f25 <- ages_incomplete_high %>% filter(fraction == "0.25")

##correcting using the exact rho

f25_high.list <- list()

for(i in 1:nrow(ages.incomplete_high.f25)){
  f25_high.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                          mu = 0.25,
                                          node_age = ages.incomplete_high.f25$Incomp.age[i],
                                          rho = 0.75)
}

##mean age
f25_high.mean.age <- sapply(f25_high.list, function(x) sum(x$time * x$prob))

##median age
f25_high.median.age <- sapply(f25_high.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))



##binding columns
ages.f25_high <- cbind(ages.incomplete_high.f25,
                     f25_high.mean.age, f25_high.median.age) %>% 
  rename(mean.age = f25_high.mean.age,
         median.age = f25_high.median.age)


# 50% missing species -----------------------------------------------------

ages.incomplete_high.f50 <- ages_incomplete_high %>% filter(fraction == "0.5")

##correcting using the exact rho

f50_high.list <- list()

for(i in 1:nrow(ages.incomplete_high.f50)){
  f50_high.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                          mu = 0.25,
                                          node_age = ages.incomplete_high.f50$Incomp.age[i],
                                          rho = 0.50)
}

##mean age
f50_high.mean.age <- sapply(f50_high.list, function(x) sum(x$time * x$prob))

##median age
f50_high.median.age <- sapply(f50_high.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))



##binding columns
ages.f50_high <- cbind(ages.incomplete_high.f50,
                     f50_high.mean.age) %>% 
  rename(mean.age = f50_high.mean.age,
         median.age = f50_high.median.age)


ages.random_high_ext <- rbind(ages.f25_high, ages.f50_high)

write_csv(ages.random_high_ext,
          file = "results/data/processed/incomplete_sampling/ages.random_high_ext.csv")

ages.random_high_ext <- read_csv(file = "results/data/processed/incomplete_sampling/ages.random_high_ext.csv")

###calculating MAPE

ages.mape.high <- ages.random_high_ext %>% group_by(fraction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Incomp.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)

ages.mape.high$extinction <- "high"


##binding mape dataframes
ages.mape.sampling_high <- ages.mape.high %>% 
  pivot_longer(cols = starts_with("mape."),
               values_to = "mape", names_to = "estimate")

ages.mape.sampling_high$estimate <- factor(ages.mape.sampling_high$estimate,
                                         levels = c("mape.incomplete",
                                                    "mape.mean",
                                                    "mape.median"),
                                         ordered = TRUE)

ages.mape.sampling_high$fraction <- as.factor(ages.mape.sampling_high$fraction)





# Fully sampled evaluation ------------------------------------------------

##reading from conservation status
ages.full_samp <- read_csv(file = "results/data/processed/incomplete_sampling/ages.full_samp.csv")


# no extinction -----------------------------------------------------------

ages.full_no <- ages.full_samp %>% filter(extinction == "no")

ages.full_no$extinction <- "no_ext"


# probability function ----------------------------------------------------

full_no.list <- list()

for(i in 1:nrow(ages.full_no)){
  full_no.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                            mu = 0,
                          node_age = ages.full_no$Estimated.age[i])
}

##mean age
full_no.mean.age <- sapply(full_no.list, function(x) sum(x$time * x$prob))

##median age
full_no.median.age <- sapply(full_no.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))

##dataframe
ages.full_no_prob <- cbind(ages.full_no,
                           full_no.mean.age, full_no.median.age) %>% 
                     rename(mean.age = full_no.mean.age,
                            median.age = full_no.median.age)


# intermediate extinction -------------------------------------------------

ages.full_int <- ages.full_samp %>% filter(extinction == "intermediate")


# probability function ----------------------------------------------------

full_int.list <- list()

for(i in 1:nrow(ages.full_int)){
  full_int.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                           mu = 0.15,
                                           node_age = ages.full_int$Estimated.age[i])
}

##mean age
full_int.mean.age <- sapply(full_int.list, function(x) sum(x$time * x$prob))

##median age
full_int.median.age <- sapply(full_int.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))

##dataframe
ages.full_int_prob <- cbind(ages.full_int,
                           full_int.mean.age, full_int.median.age) %>% 
  rename(mean.age = full_int.mean.age,
         median.age = full_int.median.age)


# High extinction ---------------------------------------------------------

ages.full_high <- ages.full_samp %>% filter(extinction == "high")


# probability function ----------------------------------------------------

full_high.list <- list()

for(i in 1:nrow(ages.full_high)){
  full_high.list[[i]] <- get_sp_age_prob_new(lam = 0.3,
                                            mu = 0.25,
                                            node_age = ages.full_high$Estimated.age[i])
}

##mean age
full_high.mean.age <- sapply(full_high.list, function(x) sum(x$time * x$prob))

##median age
full_high.median.age <- sapply(full_high.list, function(x)
  weighted.median(x = x$time,
                  w = x$prob))

##dataframe
ages.full_high_prob <- cbind(ages.full_high,
                            full_high.mean.age, full_high.median.age) %>% 
  rename(mean.age = full_high.mean.age,
         median.age = full_high.median.age)

##binding dataframes
ages.full_total_prob <- rbind(ages.full_no_prob,
                              ages.full_int_prob,
                              ages.full_high_prob)

write_csv(ages.full_total_prob,
      file = "results/data/processed/incomplete_sampling/ages.full_total_prob.csv")

ages.full_total_prob <- 
  read_csv(file = "results/data/processed/incomplete_sampling/ages.full_total_prob.csv")
####MAPE

ages.mape.full <- ages.full_total_prob %>% group_by(extinction, tree) %>% 
  summarise(mape.incomplete = 
              mean(abs(True.age - Estimated.age)/True.age)*100,
            mape.mean = 
              mean(abs(True.age - mean.age)/True.age)*100,
            mape.median = 
              mean(abs(True.age - median.age)/True.age)*100)
            
ages.mape.full$fraction <- "full"


##binding mape dataframes
ages.mape.full.sampling <- ages.mape.full %>% 
              pivot_longer(cols = starts_with("mape."),
                           values_to = "mape", names_to = "estimate")

ages.mape.full.sampling$estimate <- factor(ages.mape.full.sampling$estimate,
                                           levels = c("mape.incomplete",
                                                      "mape.mean",
                                                      "mape.median"),
                                           ordered = TRUE)

ages.mape.full.sampling$extinction <- factor(ages.mape.full.sampling$extinction,
                                             levels = c("no_ext",
                                                        "intermediate",
                                                        "high"))
ages.mape.full.sampling$fraction <- "full"




###binding all age mape
ages.mape.total <- rbind(ages.mape.full, ages.mape.no_ext,
                         ages.mape.int, ages.mape.high)

ages.mape.total$extinction <- factor(ages.mape.total$extinction,
                                     levels = c("no_ext",
                                                "intermediate",
                                                "high"),
                                     ordered = TRUE)

ages.mape.total$fraction <- factor(ages.mape.total$fraction,
                                     levels = c("full",
                                                "0.25",
                                                "0.5"),
                                     ordered = TRUE)

######binding all ages.mape.sampling df
ages.mape.sampling.tot <- rbind(ages.mape.sampling_no, ages.mape.sampling_int,
                                ages.mape.sampling_high, ages.mape.full.sampling)

##arranging factor
ages.mape.sampling.tot$extinction <- factor(ages.mape.sampling.tot$extinction,
                                            levels = c("no_ext",
                                                       "intermediate",
                                                       "high"))
ages.mape.sampling.tot$fraction <- factor(ages.mape.sampling.tot$fraction,
                                          levels = c("full",
                                                     "0.25",
                                                     "0.5"))

# plots -------------------------------------------------------------------
## themes
mynamestheme <- theme(strip.text = element_text(family = "serif", size = (9)),
                      plot.title = element_text(family = "serif", size = (12),
                                                face = "bold", hjust = 0.5),
                      axis.title = element_text(family = "serif", size = (11),
                                                face = "bold"),
                      axis.text = element_text(family = "serif", size = (9)),
                      legend.title = element_text(family = "serif", size = (9),
                                                  face = "bold"),
                      legend.text = element_text(family = "serif", size = (9)),
                      legend.background = element_rect(fill="grey",
                                                       size=.5, linetype="dotted"),
                      legend.position = "bottom")

##########MAPE of phylogenetic age

png("text/figures/MAPE.sampling.extinction.png", 
    width = 17, height = 12, units = "cm", 
    pointsize = 8, res = 300)


ggplot(ages.mape.total, aes(x = fraction, y = mape.incomplete,
                            fill = fraction))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~extinction, labeller = as_labeller(c("no_ext" = "No extinction",
                                                   "intermediate" = "Intermediate extinction",
                                                   "high" = "High extinction")))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    name="Missing species",
                    breaks=c("full", "0.25", "0.5"),
                    labels=c("0%", "25%",
                             "50%"))+
  ylim(0,1050)+
  ylab("MAPE")+
  xlab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank(),
        legend.background = element_rect(fill=NA,
                                         size=.5, linetype="dotted"),
        legend.position = c(0.1, 0.8))

dev.off()


##Fully sampled

full.plot <- ages.mape.sampling.tot %>% filter(fraction == "full") %>% 
         ggplot(aes(y = mape, x = estimate,
                                          fill = estimate))+
          geom_boxplot(outlier.shape = NA)+
        #geom_jitter(size = 0.1, alpha = 0.3)+
        facet_grid(fraction~extinction,
                   labeller = as_labeller(c("no_ext" = "No extinction",
                                            "intermediate" = "Intermediate extinction",
                                            "high" = "High extinction",
                                            "full" = "Fully sampled")))+
        scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                          name="Estimation",
                          breaks=c("mape.incomplete", "mape.mean", "mape.median"),
                          labels=c("Phylogenetic", "Mean",
                                   "Median"))+
        ylim(0, 130)+
        ylab(NULL)+
        xlab(NULL)+
        theme_bw()+
        mynamestheme+
        theme(axis.text.x = element_blank(),
              legend.position = "none")


#### 25% missing species

f25.plot <- ages.mape.sampling.tot %>% filter(fraction == "0.25") %>% 
  ggplot(aes(y = mape, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.25" = "25% missing species")))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    name="Estimation",
                    breaks=c("mape.incomplete", "mape.mean", "mape.median"),
                    labels=c("Phylogenetic", "Mean",
                             "Median"))+
  ylim(0, 350)+
  ylab(NULL)+
  xlab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank(),
        legend.position = "none")

#### 50% missing species

f50.plot <- ages.mape.sampling.tot %>% filter(fraction == "0.5") %>% 
  ggplot(aes(y = mape, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.5" = "50% missing species")))+
  scale_fill_manual(values = c("#1b9e77", "#d95f02", "#7570b3"),
                    name="Estimation",
                    breaks=c("mape.incomplete", "mape.mean", "mape.median"),
                    labels=c("Phylogenetic", "Mean",
                             "Median"))+
  ylim(0, 1050)+
  ylab(NULL)+
  xlab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank())



##plot grid
png("text/figures/MAPE.ext.samp.correct.png", 
    width = 17, height = 20, units = "cm", 
    pointsize = 8, res = 300)


##left
left <- text_grob("MAPE",
                    family = "serif", size = 13,  face = "bold",
                  rot = 90)

##grid for plotting both figures
grid.arrange(full.plot, f25.plot, f50.plot,
             nrow = 3, left = left)

dev.off()


###DELTA FIGURE

ages.delta <- ages.mape.total %>% group_by(fraction, extinction) %>%  
  summarise(delta.mean = mape.mean - mape.incomplete,
            delta.median =  mape.median - mape.incomplete) %>% 
      pivot_longer(cols = starts_with("delta."),
               values_to = "delta", names_to = "estimate")

##plots

##Fully sampled
full.delta <- ages.delta %>% filter(fraction == "full") %>% 
      ggplot(aes(y = delta, x = estimate,
                 fill = estimate))+
        geom_boxplot(outlier.shape = NA)+
        facet_grid(fraction~extinction,
                   labeller = as_labeller(c("no_ext" = "No extinction",
                                            "intermediate" = "Intermediate extinction",
                                            "high" = "High extinction",
                                            "full" = "Fully sampled"
                                           )),
                   scale = "free_y")+
        scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                          name="Estimation",
                          breaks=c("delta.mean", "delta.median"),
                          labels=c("Mean",
                                   "Median"))+
        geom_hline(yintercept = 0, linetype = "dashed",
                   size = 1, colour = "red")+
        #ylab(expression(bold(Delta* "MAPE")))+
        ylim(-100, 10)+
        xlab(NULL)+
        ylab(NULL)+
        theme_bw()+
        mynamestheme+
        theme(axis.text.x = element_blank(),
              legend.position = "none")


##25% missing species
f25.delta <- ages.delta %>% filter(fraction == "0.25") %>% 
  ggplot(aes(y = delta, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.25" = "25% missing species"
             )))+
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                    name="Estimation",
                    breaks=c("delta.mean", "delta.median"),
                    labels=c("Mean",
                             "Median"))+
  geom_hline(yintercept = 0, linetype = "dashed",
             size = 1, colour = "red")+
  #ylab(expression(bold(Delta* "MAPE")))+
  ylim(-350, 20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank(),
        legend.position = "none")

##50% missing species
f50.delta <- ages.delta %>% filter(fraction == "0.5") %>% 
  ggplot(aes(y = delta, x = estimate,
             fill = estimate))+
  geom_boxplot(outlier.shape = NA)+
  facet_grid(fraction~extinction,
             labeller = as_labeller(c("no_ext" = "No extinction",
                                      "intermediate" = "Intermediate extinction",
                                      "high" = "High extinction",
                                      "0.5" = "50% missing species"
             )))+
  scale_fill_manual(values = c("#d8b365", "#5ab4ac"),
                    name="Estimation",
                    breaks=c("delta.mean", "delta.median"),
                    labels=c("Mean",
                             "Median"))+
  geom_hline(yintercept = 0, linetype = "dashed",
             size = 1, colour = "red")+
  #ylab(expression(bold(Delta* "MAPE")))+
  ylim(-1000, 20)+
  xlab(NULL)+
  ylab(NULL)+
  theme_bw()+
  mynamestheme+
  theme(axis.text.x = element_blank())

##plot grid

png("text/figures/delta.MAPE.samp.correct.png", 
    width = 17, height = 20, units = "cm", 
    pointsize = 8, res = 300)


##left
left <- text_grob(expression(bold(Delta* "MAPE")),
                  family = "serif", size = 13,  face = "bold",
                  rot = 90)

##grid for plotting both figures
grid.arrange(full.delta, f25.delta, f50.delta,
             nrow = 3, left = left)

dev.off()
