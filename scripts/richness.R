# > RICHNESS ####
# Gather observations into narrow dataset
obs <- data.table(gather(sample.bdms, Taxon, Obs, -SiteId, -SampId))

rank <- lapply(unique(obs$Taxon), function(j){
  taxon.rank(j)
})
rank <- rbindlist(rank)

# Obtain all observations along with full taxonomy and taxonomic level of identification
obs <- merge(obs, rank, by = "Taxon", all.x = TRUE)

# Table: Community richness by taxonomic resolution
# Overall richness? EPT species richness? Family richness?
# Taxa richness
rank %>%
  summarise(n = n_distinct(Taxon)) 
# EPT species richness
rank %>%
  filter(Rank == 'Species' , Order == 'Ephemeroptera' | Order == 'Plecoptera' | Order == 'Trichoptera') %>%
  summarise(n = n_distinct(Species)) 
# Family richness (including lower levels)
rank %>%
  filter(Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
  summarise(n = n_distinct(Family)) 
# Family richness (excluding lower levels)
rank %>%
  filter(Rank == "Family") %>%
  summarise(n = n_distinct(Family)) 
# Number of taxa ID'ed at each taxonomic rank
rank %>%
  group_by(Rank) %>%
  summarise(n = n_distinct(Taxon))

# Plot - P1 map richness ####
# Observed site richness
obs.site.richness.species <- obs %>%
  filter(Obs == 1, Rank == "Species") %>%
  group_by(SiteId) %>%
  mutate(ObsSpecies = n_distinct(Species)) %>%
  ungroup() %>%
  select(SiteId, ObsSpecies) %>%
  distinct()

obs.site.richness <- obs %>%
  filter(Obs == 1, Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
  group_by(SiteId) %>%
  mutate(ObsFamilies = n_distinct(Family)) %>%
  ungroup() %>% # remove grouping information from data frame
  select(SiteId, ObsFamilies) %>%
  distinct() %>%
  left_join(obs.site.richness.species, by = "SiteId") %>%
  left_join(inputs$xy, by = "SiteId")
rm(obs.site.richness.species)

ch <- fortify(inputs$ch)
g <- ggplot()
g <- g + geom_polygon(data = ch, aes(x = long, y = lat, group = group), fill=NA, color="black")
g <- g + geom_point(data = obs.site.richness, aes(X, Y, size = ObsFamilies), alpha = 0.35)
g <- g + labs(title = "",
              x = "",
              y = "")
g <- g + scale_y_continuous(breaks=NULL)
g <- g + scale_x_continuous(breaks=NULL)
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(size = "Number of families")
g <- g + theme_minimal(base_size = 15)
g <- g + theme(plot.title = element_text(hjust = 0.5))
pdf('P1 - map family richness.pdf', paper = 'special', width = 11, onefile = TRUE)
print(g)
dev.off()

g <- ggplot()
g <- g + geom_polygon(data = ch, aes(x = long, y = lat, group = group), fill=NA, color="black")
g <- g + geom_point(data = obs.site.richness, aes(X, Y, size = ObsSpecies), alpha = 0.35)
g <- g + labs(title = "",
              x = "",
              y = "")
g <- g + scale_y_continuous(breaks=NULL)
g <- g + scale_x_continuous(breaks=NULL)
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(size = "Number of species")
g <- g + theme_minimal(base_size = 15)
g <- g + theme(plot.title = element_text(hjust = 0.5))
pdf('P1 - map species richness.pdf', paper = 'special', width = 11, onefile = TRUE)
print(g)
dev.off()

# > Observed richness ####
# Observed species richness by sample
# Note n_distinct() approach doesn't return samples where richness is zero;
# use !NA and sum() since identical taxa don't occur more than once at a sample
obs.richness.species <- obs %>%
  filter(!is.na(Obs), Rank == "Species") %>%
  group_by(SampId) %>%
  # Species richness calculated by summing !NA species observations per sample
  summarise(ObsSpecies = sum(Obs))

# Observed family richness by sample
obs.richness <- obs %>%
  filter(Obs==1, Rank == "Family" | Rank == "Genus" | Rank == "Species") %>%
  group_by(SampId) %>%
  summarise(ObsFamilies = n_distinct(Family)) %>%
  full_join(obs.richness.species, SampId, by = "SampId")
rm(obs.richness.species)

# # Observed rare family richness (including genera and species)
# n.rare.families <- obs %>%
#   filter(Obs == 1, Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
#   group_by(SampId, Family) %>%
#   # For each sample and for each family, how many distinct families? (should return values of 1)
#   summarise(n = n_distinct(Family)) %>%
#   group_by(SampId, Family) %>%
#   # Ensure n_distinct() is set to 1
#   summarise(n.family = ifelse(n > 0, 1, 0)) %>%
#   # Get the number of occurrences per family
#   group_by(Family) %>% # Overall occurrences by family
#   summarise(n = sum(n.family))
# 
# n.families <- n.rare.families$n
# names(n.families) <- n.rare.families$Family
# n.families <- sort(n.families, decreasing = TRUE)
# rm(n.rare.families)

# # Observed rare species richness
# n.species <- occur.freq(sample.bdms)
# obs.rare.species <- obs  %>%
#   mutate(n = n.species[Taxon]) %>%
#   filter(!is.na(Obs), Rank == "Species", n <= 56) %>%
#   group_by(SampId) %>%
#   summarise(ObsRareSpecies = sum(Obs))

# # Observed rare family richness
# obs.rare <- obs  %>%
#   filter(Obs==1, Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
#   # Calculate the occurrence frequency based on family
#   mutate(n = n.families[Family]) %>%
#   filter(n <= 56) %>%
#   group_by(SampId) %>%
#   # n_distinct doesn't return samples with zero richness; replace NAs with zero after full_join
#   summarise(ObsRareFamilies = n_distinct(Family)) %>%
#   full_join(obs.rare.species, by = "SampId") %>%
#   mutate(ObsRareFamilies = ifelse(is.na(ObsRareFamilies), 0, ObsRareFamilies))
# 
# obs.richness <- full_join(obs.richness, obs.rare, by = "SampId")
# rm(obs.rare.species, obs.rare)

# > Predicted richness ####
# Calculate "predicted" family/species richness by sample
# Combine jSDM probabilities (w/o RE)
jsdm.re.prob <- jsdm.re$probability %>%
  select(SampId, Taxon, Pred) %>%
  mutate(RE = TRUE)

jsdm.prob <- jsdm$probability %>% 
  select(SampId, Taxon, Pred) %>%
  mutate(RE = FALSE) %>%
  bind_rows(jsdm.re.prob) %>%
  as.tibble()
rm(jsdm.re.prob)


# Predicted family richness (RE: FALSE)
jsdm.family.prob <- obs %>%
  # Filter the obs/taxonomy and join the joint model probabilities 
  filter(!is.na(Obs), Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
  left_join(filter(jsdm.prob, RE == FALSE), by = c("Taxon", "SampId")) %>%
  # Calculate the predicted richness by family at each sample
  group_by(Family, SampId) %>%
  summarise(Pred = 1-prod(1-Pred)) %>%
  ungroup() %>%
  select(SampId, Family, Pred)

jsdm.richness.family <- jsdm.family.prob %>%
  group_by(SampId) %>%
  summarise(PredFamilies = sum(Pred)) %>%
  select(SampId, PredFamilies) %>%
  mutate(RE = FALSE)

# Predicted family richness (RE: TRUE)
jsdm.family.prob <- obs %>%
  # Filter the obs/taxonomy and join the joint model probabilities 
  filter(!is.na(Obs), Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
  left_join(filter(jsdm.prob, RE == TRUE), by = c("Taxon", "SampId")) %>%
  # Calculate the predicted richness by family at each sample
  group_by(Family, SampId) %>%
  summarise(Pred = 1-prod(1-Pred)) %>%
  ungroup() %>%
  select(SampId, Family, Pred)

jsdm.richness.family <- jsdm.family.prob %>%
  group_by(SampId) %>%
  summarise(PredFamilies = sum(Pred)) %>%
  select(SampId, PredFamilies) %>%
  mutate(RE = TRUE) %>%
  bind_rows(jsdm.richness.family)
rm(jsdm.family.prob)

# Predicted total species richness by sample (RE: FALSE)
jsdm.richness.species <- obs %>%
  filter(!is.na(Obs), Rank == "Species") %>%
  left_join(filter(jsdm.prob, RE == FALSE), by = c("Taxon", "SampId")) %>%
  group_by(SampId) %>%
  summarise(PredSpecies = sum(Pred)) %>%
  select(SampId, PredSpecies) %>%
  mutate(RE = FALSE)

# Predicted total species richness by sample (RE: TRUE)
jsdm.richness.species <- obs %>%
  filter(!is.na(Obs), Rank == "Species") %>%
  left_join(filter(jsdm.prob, RE == TRUE), by = c("Taxon", "SampId")) %>%
  group_by(SampId) %>%
  summarise(PredSpecies = sum(Pred)) %>%
  select(SampId, PredSpecies) %>%
  mutate(RE = TRUE) %>%
  bind_rows(jsdm.richness.species)

jsdm.richness <- left_join(jsdm.richness.family, jsdm.richness.species, by = c("SampId", "RE"))
rm(jsdm.richness.family, jsdm.richness.species)

# # Combine predicted richness plots
# # > Rare families ####
# # Predicted rare family richness by sample (RE: FALSE)
# jsdm.family.prob <- obs %>%
#   # Filter observations by taxonomy/rarity and join the joint model probabilities 
#   filter(!is.na(Obs), Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
#   mutate(n = n.families[Family]) %>%
#   filter(n <= 56) %>%
#   left_join(filter(jsdm.prob, RE == FALSE), by = c("Taxon", "SampId")) %>%
#   # Calculate the predicted richness by family at each sample
#   group_by(Family, SampId) %>%
#   summarise(PredRareFamilies = 1-prod(1-Pred)) %>%
#   ungroup() %>%
#   select(SampId, Family, PredRareFamilies)
# 
# jsdm.richness.family <- jsdm.family.prob %>%
#   group_by(SampId) %>%
#   summarise(PredRareFamilies = sum(PredRareFamilies)) %>%
#   select(SampId, PredRareFamilies) %>%
#   mutate(RE = FALSE)
# 
# # Predicted rare family richness by sample (RE: TRUE)
# jsdm.family.prob <- obs %>%
#   # Filter observations by taxonomy/rarity and join the joint model probabilities 
#   filter(!is.na(Obs), Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
#   mutate(n = n.families[Family]) %>%
#   filter(n <= 56) %>%
#   left_join(filter(jsdm.prob, RE == TRUE), by = c("Taxon", "SampId")) %>%
#   # Calculate the predicted richness by family at each sample
#   group_by(Family, SampId) %>%
#   summarise(PredRareFamilies = 1-prod(1-Pred)) %>%
#   ungroup() %>%
#   select(SampId, Family, PredRareFamilies)
# 
# jsdm.richness.family <- jsdm.family.prob %>%
#   group_by(SampId) %>%
#   summarise(PredRareFamilies = sum(PredRareFamilies)) %>%
#   select(SampId, PredRareFamilies) %>%
#   mutate(RE = TRUE) %>%
#   bind_rows(jsdm.richness.family)
# 
# jsdm.richness <- left_join(jsdm.richness, jsdm.richness.family, by = c("SampId", "RE"))
# rm(jsdm.richness.family)
# 
# # # > Rare species ####
# # Predicted total species richness by sample (RE: FALSE)
# jsdm.richness.species <- obs %>%
#   filter(!is.na(Obs), Rank == "Species") %>%
#   mutate(n = n.species[Taxon]) %>%
#   filter(n <= 56) %>%
#   left_join(filter(jsdm.prob, RE == FALSE), by = c("Taxon", "SampId")) %>%
#   group_by(SampId) %>%
#   summarise(PredRareSpecies = sum(Pred)) %>%
#   select(SampId, PredRareSpecies) %>%
#   mutate(RE = FALSE)
# 
# # Predicted total species richness by sample (RE: TRUE)
# jsdm.richness.species <- obs %>%
#   filter(!is.na(Obs), Rank == "Species") %>%
#   mutate(n = n.species[Taxon]) %>%
#   filter(n <= 56) %>%
#   left_join(filter(jsdm.prob, RE == TRUE), by = c("Taxon", "SampId")) %>%
#   group_by(SampId) %>%
#   summarise(PredRareSpecies = sum(Pred)) %>%
#   select(SampId, PredRareSpecies) %>%
#   mutate(RE = TRUE) %>%
#   bind_rows(jsdm.richness.species)
# 
# jsdm.richness <- left_join(jsdm.richness, jsdm.richness.species, by = c("SampId", "RE"))
# rm(jsdm.richness.species)

# Combine observed and predicted richness ####
rdt <- left_join(jsdm.richness, obs.richness, by = "SampId")

# Split/gather observed and predicted richness data
rdt.pred <- rdt %>%
  select(SampId, RE, contains("Pred")) %>%
  gather(Variable, Value, -SampId, -RE)

rdt$Rank <- ""
rdt.pred$Rank[rdt.pred$Variable=="PredFamilies"] <- "Families"
rdt.pred$Rank[rdt.pred$Variable=="PredSpecies"] <- "Species"
rdt.pred$Variable <- "Predicted"

rdt.obs <- rdt %>%
  select(SampId, RE, contains("Obs")) %>%
  gather(Variable, Value, -SampId, -RE)

rdt.obs$Rank <- ""
rdt.obs$Rank[rdt.obs$Variable=="ObsFamilies"] <- "Families"
rdt.obs$Rank[rdt.obs$Variable=="ObsSpecies"] <- "Species"

rdt.obs$Variable <- "Observed"

rdt.plot <- bind_rows(rdt.pred, rdt.obs)
# rm(rdt.pred, rdt.obs)

rdt.plot$RE.label <- ifelse(rdt.plot$RE==TRUE, "With random effect", "Without random effect")
rdt.plot$RE.label <- as.factor(rdt.plot$RE.label)
rdt.plot$Rank[rdt.plot$Rank=="Families"] <- "All families"
rdt.plot$Rank[rdt.plot$Rank=="Species"] <- "EPT species"
rdt.plot$Rank <- as.factor(rdt.plot$Rank)
# rdt.plot$Rank <- factor(rdt.plot$Rank, levels = c("Families" = "All families", "Species" = "EPT species"))

rdt.plot <- spread(rdt.plot, Variable, Value)
setDT(rdt.plot)


# Visually check that subset has the 4 distinct cases in plot
rdt.rmse <- rdt.plot[1:4, ]
rmse.val <- c()
test <- filter(rdt.plot, RE==FALSE, Rank=="All families")
rmse.val[1] <- rmse(test$Observed-test$Predicted)
test <- filter(rdt.plot, RE==FALSE, Rank=="EPT species")
rmse.val[2] <- rmse(test$Observed-test$Predicted)
test <- filter(rdt.plot, RE==TRUE, Rank=="All families")
rmse.val[3] <- rmse(test$Observed-test$Predicted)
test <- filter(rdt.plot, RE==TRUE, Rank=="EPT species")
rmse.val[4] <- rmse(test$Observed-test$Predicted)
rdt.rmse$rmse <- rmse.val
rdt.rmse$Label <- paste("RMSE =", round(rdt.rmse$rmse, 2))
rm(test, rmse.val) # remove this shameful hack from my sight...

rdt.rmse$Observed <- 10
rdt.rmse$Predicted <- 30

# Plot all richness data ####
g <- ggplot(data = rdt.plot, aes(x = Observed, y = Predicted))
g <- g + geom_point(alpha = 0.35)
g <- g + facet_grid(RE.label ~ Rank, scales="free_x")
# g <- g + facet_grid(Rank ~ RE.label)
g <- g + geom_smooth(method = 'lm', se=FALSE)
g <- g + geom_text(data=rdt.rmse, aes(x = Observed, y = Predicted, label = Label), size = 6, parse = FALSE, inherit.aes=FALSE)
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
g <- g + theme_bw(base_size = 18)
g <- g + theme(axis.text=element_text(size=18))
g <- g + labs(x = "Observed sample richness",
              y = "Predicted sample richness")
print(g)

# Plot % rare vs richness ####
prare.ps <- rdt %>%
  filter(RE==T) %>%
  mutate(Fraction = PredRareSpecies/PredSpecies,
         Richness = PredSpecies,
         Type = "Predicted",
         Rank = "Species") %>%
  select(SampId, Fraction, Richness, Type, Rank)

prare.pf <- rdt %>%
  filter(RE==T) %>%
  mutate(Fraction = PredRareFamilies/PredFamilies,
         Richness = PredFamilies,
         Type = "Predicted",
         Rank = "Family") %>%
  select(SampId, Fraction, Richness, Type, Rank)

prare.os <- rdt %>%
  filter(RE==T) %>%
  mutate(Fraction = ObsRareSpecies/ObsSpecies,
         Richness = ObsSpecies,
         Type = "Observed",
         Rank = "Species") %>%
  select(SampId, Fraction, Richness, Type, Rank)

prare.of <- rdt %>%
  filter(RE==T) %>%
  mutate(Fraction = ObsRareFamilies/ObsFamilies,
         Richness = ObsFamilies,
         Type = "Observed",
         Rank = "Family") %>%
  select(SampId, Fraction, Richness, Type, Rank)

prare <- bind_rows(prare.ps, prare.pf, prare.os, prare.of)
rm(prare.ps, prare.pf, prare.os, prare.of)

setDT(prare)
lm.data <- data.table(Type = rep(c("Observed", "Predicted"), 2),
                      Rank = c("Family", "Family", "Species", "Species"),
                      Eqn = c(write_lm("Richness", "Fraction", prare[Type=="Observed" & Rank =="Family",]),
                              write_lm("Richness", "Fraction", prare[Type=="Predicted" & Rank =="Family",]),
                              write_lm("Richness", "Fraction", prare[Type=="Observed" & Rank =="Species",]),
                              write_lm("Richness", "Fraction", prare[Type=="Predicted" & Rank =="Species",])),
                      x = rep(22, 4),
                      y = rep(0.8, 4))

g <- ggplot(prare, aes(x = Richness, y = Fraction))
g <- g + geom_point(alpha = 0.35)
g <- g + facet_grid(Rank ~ Type)
g <- g + geom_smooth(method='lm')
g <- g + geom_text(data=lm.data, aes(x = x, y = y, label = Eqn), parse = TRUE, inherit.aes=FALSE)
g <- g + theme_bw(base_size = 17)
g <- g + labs(title="Proportion of rare taxa richness by sample in joint model with random effects",
              x = "Richness per sample",
              y = "Rare taxa richness as fraction of total richness (rare/all)")
g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
print(g)

rdt.plot %>%
  group_by(RE, Rank) %>%
  summarise(min = min(Predicted),
            mean = mean(Predicted),
            median = median(Predicted),
            max = max(Predicted)) %>%
  arrange(Rank)

# iSDM richness ####
# Pred family richness ####
isdm.family.prob <- obs %>%
  # Filter the obs/taxonomy and join the joint model probabilities 
  filter(!is.na(Obs), Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
  left_join(bdm.isdm.prob, by = c("Taxon", "SampId")) %>%
  # Calculate the predicted richness by family at each sample
  group_by(Family, SampId) %>%
  summarise(Pred = 1-prod(1-Pred)) %>%
  ungroup() %>%
  select(SampId, Family, Pred)

isdm.richness.family <- isdm.family.prob %>%
  group_by(SampId) %>%
  summarise(PredFamilies = sum(Pred)) %>%
  select(SampId, PredFamilies)

# Pred species richness ####
isdm.richness.species <- obs %>%
  filter(!is.na(Obs), Rank == "Species") %>%
  left_join(bdm.isdm.prob, by = c("Taxon", "SampId")) %>%
  group_by(SampId) %>%
  summarise(PredSpecies = sum(Pred)) %>%
  select(SampId, PredSpecies)

isdm.richness <- left_join(isdm.richness.family, isdm.richness.species, by = c("SampId"))
rm(isdm.richness.family, isdm.richness.species)

# > Rare families ####
# Pred rare family richness
isdm.family.prob <- obs %>%
  # Filter observations by taxonomy/rarity and join the joint model probabilities 
  filter(!is.na(Obs), Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
  mutate(n = n.families[Family]) %>%
  filter(n <= 56) %>%
  left_join(bdm.isdm.prob, by = c("Taxon", "SampId")) %>%
  # Calculate the predicted richness by family at each sample
  group_by(Family, SampId) %>%
  summarise(PredRareFamilies = 1-prod(1-Pred)) %>%
  ungroup() %>%
  select(SampId, Family, PredRareFamilies)

isdm.richness.family <- isdm.family.prob %>%
  group_by(SampId) %>%
  summarise(PredRareFamilies = sum(PredRareFamilies)) %>%
  select(SampId, PredRareFamilies)

isdm.richness <- left_join(isdm.richness, isdm.richness.family, by = c("SampId"))
rm(isdm.richness.family)



# > Rare species ####
# Predicted total species richness by sample (RE: FALSE)
isdm.richness.species <- obs %>%
  filter(!is.na(Obs), Rank == "Species") %>%
  mutate(n = n[Taxon]) %>%
  filter(n <= 56) %>%
  left_join(bdm.isdm.prob, by = c("Taxon", "SampId")) %>%
  group_by(SampId) %>%
  summarise(PredRareSpecies = sum(Pred)) %>%
  select(SampId, PredRareSpecies)

# Predicted total species richness by sample (RE: TRUE)
isdm.richness.species <- obs %>%
  filter(!is.na(Obs), Rank == "Species") %>%
  mutate(n = n[Taxon]) %>%
  filter(n <= 56) %>%
  left_join(bdm.isdm.prob, by = c("Taxon", "SampId")) %>%
  group_by(SampId) %>%
  summarise(PredRareSpecies = sum(Pred)) %>%
  select(SampId, PredRareSpecies) %>%
  bind_rows(isdm.richness.species)

isdm.richness <- left_join(isdm.richness, isdm.richness.species, by = c("SampId"))
rm(isdm.richness.species)

isdm.richness <- distinct(isdm.richness)
isdm.richness <- left_join(isdm.richness, obs.richness, by = "SampId")

plot(isdm.richness$ObsFamilies, isdm.richness$PredFamilies)
abline(a = 0, b = 1)

# Test data for plotting
test <- jsdm.richness %>%
  filter(RE==FALSE) %>%
  select(-RE)

test$Model <- "jSDM"
isdm.richness$Model <- "iSDM"

test <- bind_rows(test, isdm.richness)
test <- left_join(test, obs.richness, by="SampId")

# PLOTS ####
# All families ####
setDT(test)
lm.data <- data.table(Model = c("iSDM", "jSDM"),
                      Eqn = c(write_lm("ObsFamilies", "PredFamilies", test[Model=="iSDM",]),
                              write_lm("ObsFamilies", "PredFamilies", test[Model=="jSDM",])),
                      x = c(7,7),
                      y = c(28, 31))

g <- ggplot(test, aes(x=ObsFamilies, y=PredFamilies, color=Model))
g <- g + geom_point(alpha=0.4)
g <- g + geom_smooth(method="lm")
# g <- g + geom_text(data=lm.data, aes(x = x, y = y, label = Eqn, color=Model), parse = TRUE, inherit.aes=FALSE)
print(g)

# All species ####
lm.data <- data.table(Model = c("iSDM", "jSDM"),
                      Eqn = c(write_lm("ObsSpecies", "PredSpecies", test[Model=="iSDM",]),
                              write_lm("ObsSpecies", "PredSpecies", test[Model=="jSDM",])),
                      x = c(7,7),
                      y = c(28, 31))

g <- ggplot(test, aes(x=ObsSpecies, y=PredSpecies, color=Model))
g <- g + geom_point(alpha=0.4)
g <- g + geom_smooth(method="lm")
g <- g + geom_text(data=lm.data, aes(x = x, y = y, label = Eqn, color=Model), parse = TRUE, inherit.aes=FALSE)
print(g)

# Rare families ####
lm.data <- data.table(Model = c("iSDM", "jSDM"),
                      Eqn = c(write_lm("ObsRareFamilies", "PredRareFamilies", test[Model=="iSDM",]),
                              write_lm("ObsRareFamilies", "PredRareFamilies", test[Model=="jSDM",])),
                      x = c(7,7),
                      y = c(28, 31))

g <- ggplot(test, aes(x=ObsRareFamilies, y=PredRareFamilies, color=Model))
g <- g + geom_point(alpha=0.4)
g <- g + geom_smooth(method="lm")
g <- g + geom_text(data=lm.data, aes(x = x, y = y, label = Eqn, color=Model), parse = TRUE, inherit.aes=FALSE)
print(g)

# Rare species ####
lm.data <- data.table(Model = c("iSDM", "jSDM"),
                      Eqn = c(write_lm("ObsRareSpecies", "PredRareSpecies", test[Model=="iSDM",]),
                              write_lm("ObsRareSpecies", "PredRareSpecies", test[Model=="jSDM",])),
                      x = c(7,7),
                      y = c(28, 31))

g <- ggplot(test, aes(x=ObsRareSpecies, y=PredRareSpecies, color=Model))
g <- g + geom_point(alpha=0.4)
g <- g + geom_smooth(method="lm")
g <- g + geom_text(data=lm.data, aes(x = x, y = y, label = Eqn, color=Model), parse = TRUE, inherit.aes=FALSE)
print(g)

# Combined all families/species plot
test.family <- test[, c("SampId", "PredFamilies", "ObsFamilies", "Model"), with=F]
colnames(test.family) <- c("SampId", "Pred", "Obs", "Model")
test.family$Rank <- "Family"

test.species <- test[, c("SampId", "PredSpecies", "ObsSpecies", "Model"), with=F]
colnames(test.species) <- c("SampId", "Pred", "Obs", "Model")
test.species$Rank <- "Species"

plot.data <- bind_rows(test.family, test.species)
rm(test.family, test.species)

lm.family <- data.table(Model = c("iSDM", "jSDM"),
                        Rank = c("Family", "Family"),
                        Eqn = c(write_lm("Obs", "Pred", plot.data[Model=="iSDM" & Rank=="Family",]),
                                write_lm("Obs", "Pred", plot.data[Model=="jSDM" & Rank=="Family",])),
                        x = c(9,9),
                        y = c(28, 31))

lm.species <- data.table(Model = c("iSDM", "jSDM"),
                         Rank = c("Species", "Species"),
                         Eqn = c(write_lm("Obs", "Pred", plot.data[Model=="iSDM" & Rank=="Species",]),
                                 write_lm("Obs", "Pred", plot.data[Model=="jSDM" & Rank=="Species",])),
                         x = c(9,9),
                         y = c(20, 23))

lm.data <- bind_rows(lm.family, lm.species)
rm(lm.family, lm.species)

g <- ggplot(plot.data, aes(x=Obs, y=Pred, color=Model))
g <- g + geom_point(alpha=0.3)
g <- g + facet_wrap(~ Rank, scales="free")
g <- g + geom_smooth(method="lm")
g <- g + geom_text(data=lm.data, aes(x = x, y = y, label = Eqn, color=Model), size = 6, parse = TRUE, inherit.aes=FALSE, show.legend = FALSE)
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)

g <- g + labs(x = "Observed richness", y = "Predicted richness")
g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
g <- g + theme(axis.text=element_text(size=18))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + theme_bw(base_size = 18)
print(g)

# # Plot random effect vs influence factors
# re <- data.table(SiteId = names(bdm.ch.re$eps), RandomEffect = bdm.ch.re$eps)
# re.data <- left_join(predictors, re, by = "SiteId")
# re.data <- select(re.data, SampId, RandomEffect, Temp, Temp2, FV, F10m, IAR, Urban, LUD)
# 
# # pairs(select(re.data, -SampId), pch = 20)
# pdf('pairwise_scatterplot_RE.pdf')
# pairs.panels(select(re.data, -SampId), cex.labels=1.5, hist.col="grey")
# dev.off()

# Plot: % rare vs RE
# name:value vector of site:sample pairs
ss <- sample.bdms$SiteId
names(ss) <- sample.bdms$SampId

prare$Eps <- bdm.ch.re$eps.maxpost[ss[prare$SampId]]

dt <- filter(prare, Type =="Predicted", Rank=="Species")
g <- ggplot(data = dt, aes(x = Eps, y = Fraction))
g <- g + geom_point(alpha = 0.3)
# g <- g + geom_density2d(color = "red", size = 1.3)
g <- g + geom_smooth(method="lm")
g <- g + theme_bw()
print(g)

# Richness model ####
predictors.p1 <- prepare.inputs(c("Temp","Temp2", "FV", "F10m", "IAR", "Urban", "LUD"), y = sample.bdms, center = TRUE)

rm.data <- rdt %>%
  filter(RE == FALSE) %>%
  select(-PredFamilies, -PredSpecies) %>% # Remove joint model results
  left_join(select(predictors.p1, -SiteId), by="SampId") %>%
  as.tibble()

# Fit the negative binomial models: family and species
rmf <- glm.nb(ObsFamilies ~ Temp + Temp2 + FV + F10m + IAR + Urban + LUD,data=rm.data,link=log)
rms <- glm.nb(ObsSpecies ~ Temp + Temp2 + FV + F10m + IAR + Urban + LUD,data=rm.data,link=log)

rm.data$PredFamilies <- exp(predict(rmf))
rm.data$PredSpecies <- exp(predict(rms))

rm <- rm.data %>%
  select(SampId, RE, PredFamilies, ObsFamilies, PredSpecies, ObsSpecies) %>%
  gather(Variable, Value, -SampId, -RE)

rm$Rank[rm$Variable=="PredFamilies"] <- "All families"
rm$Rank[rm$Variable=="ObsFamilies"] <- "All families"
rm$Rank[rm$Variable=="PredSpecies"] <- "EPT species"
rm$Rank[rm$Variable=="ObsSpecies"] <- "EPT species"


rm$Variable[rm$Variable=="PredFamilies"] <- "Predicted"
rm$Variable[rm$Variable=="ObsFamilies"] <- "Observed"
rm$Variable[rm$Variable=="PredSpecies"] <- "Predicted"
rm$Variable[rm$Variable=="ObsSpecies"] <- "Observed"

rm <- spread(rm, Variable, Value)
rm$RE.label <- "RM"

rm.rmse <- data.table(Rank = c("All families", "EPT species"),
                      Predicted = c(25, 25),
                      Observed = c(10,10))
rm.rmse$rmse <- c(rmse(filter(rm, Rank == "All families")$Observed - filter(rm, Rank == "All families")$Predicted),
                  rmse(filter(rm, Rank == "EPT species")$Observed - filter(rm, Rank == "EPT species")$Predicted))
rm.rmse$Label <- paste("RMSE =", round(rm.rmse$rmse, 2))

g <- ggplot(data = rm, aes(x = Observed, y = Predicted))
g <- g + geom_point(alpha = 0.35)
# g <- g + facet_grid(. ~ Rank, scales="free_x")
g <- g + facet_wrap(~ Rank, scales="free_x")
g <- g + geom_smooth(method = 'lm', se=FALSE)
g <- g + geom_text(data=rm.rmse, aes(x = Observed, y = Predicted, label = Label), size = 6, parse = FALSE, inherit.aes=FALSE)
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
g <- g + theme_bw(base_size = 18)
g <- g + theme(axis.text=element_text(size=18))
g <- g + labs(x = "Observed sample richness",
              y = "Predicted sample richness")
print(g)

# max.families <- max(test$ObsFamilies)
# plot(test$ObsFamilies,exp(predict(nbf)),xlim=c(0,max.families),ylim=c(0,max.families),main="Families",
#      xlab="Observed Richness",ylab="Predicted Richness")
# lines(c(0,max.families),c(0,max.families))
# abline(lm.nbf)

# Families
# Negative binomial GLM
rmse(rm.data$ObsFamilies-exp(predict(rmf)))
# Joint model
rmse(rm.data$ObsFamilies-rm.data$PredFamilies)

# max.species <- max(test$ObsSpecies)
# plot(test$ObsSpecies,exp(predict(nbs)),xlim=c(0,max.species),ylim=c(0,max.species),main="EPT Species",
#      xlab="Observed Richness",ylab="Predicted Richness")
# lines(c(0,max.species),c(0,max.species))
# abline(lm.nbs)

# Species
# Negative binomial GLM
rmse(test$ObsSpecies-exp(predict(nbs)))
# Joint model
rmse(test$ObsSpecies-test$PredSpecies)

# Write the table to disk
dt <- data.table(Rank = c("Family", "Family", "Species", "Species"), 
                 Model = c("nbGLM", "jSDM", "nbGLM", "jSDM"),
                 RMSE = c(rmse(test$ObsFamilies-exp(predict(nbf))), rmse(test$ObsFamilies-test$PredFamilies), rmse(test$ObsSpecies-exp(predict(nbs))), rmse(test$ObsSpecies-test$PredSpecies)))

dt$nRMSE <- c(rep(mean(test$ObsFamilies),2), rep(mean(test$ObsSpecies),2))
dt$nRMSE <- dt$RMSE/dt$nRMSE

write.csv(dt, 'RMSE_richness.csv', row.names = F)
