# deploy.jsdm() ####
K <- c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "WV", "Temp", "Temp2")
deploy.jsdm(K, sample.bdms, "bdms", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "bdms_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "bdms_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "bdms_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdmf, "bdmf", center=T, cv=0)
deploy.jsdm(K, sample.bdmf, "bdmf_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdmf, "bdmf_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdmf, "bdmf_train3", center=T, cv=3)

deploy.jsdm(K, sample.invf, "invf", center=T, cv=0)
deploy.jsdm(K, sample.invf, "invf_train1", center=T, cv=1)
deploy.jsdm(K, sample.invf, "invf_train2", center=T, cv=2)
deploy.jsdm(K, sample.invf, "invf_train3", center=T, cv=3)

deploy.jsdm(K, sample.invfp, "invfp", center=T, cv=0)
deploy.jsdm(K, sample.invfp, "invfp_train1", center=T, cv=1)
deploy.jsdm(K, sample.invfp, "invfp_train2", center=T, cv=2)
deploy.jsdm(K, sample.invfp, "invfp_train3", center=T, cv=3)

# occur.freq(sample) ####
n.bdmf <- occur.freq(sample.bdmf)
n.bdms <- occur.freq(sample.bdms)
n.invf <- occur.freq(sample.invf)
n.invfp <- occur.freq(sample.invfp)

# get.taxonomy(sample) ####
taxonomy.bdmf <- get.taxonomy(sample.bdmf)
taxonomy.bdms <- get.taxonomy(sample.bdms)
taxonomy.invf <- get.taxonomy(sample.invf)

taxonomy.ept <- taxonomy.bdms[Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera"),]

# extract.jsdm(results) ####
# Results: paper 1 BDM species
p1.results.re <- extract.jsdm('outputs/jsdm_p1', 'RE', epsilon=TRUE)
p1.results <- extract.jsdm('outputs/jsdm_p1', 'noRE', epsilon=FALSE)

# Results: paper 2 BDM family/species
# jsdm.results <- extract.jsdm(trials)
jsdm.bdmf <- extract.jsdm("outputs/jsdm_p2", "bdmf", epsilon=FALSE)
jsdm.bdms <- extract.jsdm("outputs/jsdm_p2", "bdms", epsilon=FALSE)

# Calibration: combined families
jsdm.invf <- extract.jsdm("outputs/jsdm_p2", "invf")
jsdm.invfp <- extract.jsdm("outputs/jsdm_p2", "invf_plat")

# cv.jsdm(jsdm.results) ####
jsdm.cv.bdm <- cv.jsdm(c("bdmf", "bdms"))
jsdm.cv.inv <- cv.jsdm(c("invf", "invf_plat"))

# map.jsdm, plot.prob ####
map.jsdm(jsdm.bdmf, 'P2 maps bdmf')
map.jsdm(jsdm.bdms, 'P2 maps bdms')
map.jsdm(jsdm.invf, 'P2 maps invf')
map.jsdm(jsdm.invfp, 'P2 maps invf plateau')

plot.prob(jsdm.bdmf, 'P2 prob bdmf')
plot.prob(jsdm.invf, 'P2 prob invf')
plot.prob(jsdm.bdms, 'P2 prob bdms')
plot.prob(jsdm.invfp, 'P2 prob invf plateau')

# plot.comm(results) ####
plot.comm(jsdm.bdmf, 'P2 parameters bdmf')
plot.comm(jsdm.bdms, 'P2 parameters bdms')
plot.comm(jsdm.invf, 'P2 parameters invf')
plot.comm(jsdm.invfp, 'P2 parameters invf plateau')

# parameter distributions
jsdm.beta <- bind_rows(jsdm.bdmf$parameters, jsdm.invf$parameters, jsdm.invfp$parameters)

# do the same taxa show different responses in different spatial distributions?
plot.data <- jsdm.beta[Taxon %in% Reduce(intersect, list(names(n.bdmf),names(n.invf),names(n.invfp))),]
# slopes.plot$Label <- factor(slopes.plot$Variable, levels = v)
# levels(slopes.plot$Label) <- labeller(levels(slopes.plot$Label))

plot.data$Label <- factor(plot.data$Variable, levels = K)
levels(plot.data$Label) <- labeller(levels(plot.data$Label))
plot.data$Trial[plot.data$Trial=="bdmf"] <- "BDM families"
plot.data$Trial[plot.data$Trial=="invf"] <- "Combined families CH"
plot.data$Trial[plot.data$Trial=="invf_plat"] <- "Combined families plateau"

g <- ggplot(data=plot.data)
g <- g + geom_density(aes(x = Parameter, fill = Trial), alpha = 0.4)
g <- g + facet_wrap(~ Label, scales='free', labeller=label_parsed)
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + theme_bw(base_size=16)
g <- g + labs(title="Max. posterior taxon-specific distributions of 92 common families",
              x="Influence Factor", y="Density")
print(g)

# Plot: scatterplots ####
pdf('P2 scatterplot by trial.pdf', onefile = TRUE)
pairs.panels(select(prepare.inputs(c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp", "Temp2"), sample.bdmf, center=FALSE), -SiteId, -SampId), density = TRUE, scale=FALSE, hist.col="grey", cex.cor=1.5, cex.labels=1.5, main = "BDM families")
pairs.panels(select(prepare.inputs(c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp", "Temp2"), sample.bdms, center=FALSE), -SiteId, -SampId), density = TRUE, scale=FALSE, hist.col="grey", cex.cor=1.5, cex.labels=1.5, main = "BDM species")
pairs.panels(select(prepare.inputs(c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp", "Temp2"), sample.invf, center=FALSE), -SiteId, -SampId), density = TRUE, scale=FALSE, hist.col="grey", cex.cor=1.5, cex.labels=1.5, main = "combined families CH")
pairs.panels(select(prepare.inputs(c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp", "Temp2"), sample.invf.plat, center=FALSE), -SiteId, -SampId), density = TRUE, scale=FALSE, hist.col="grey", cex.cor=1.5, cex.labels=1.5, main = "combined families plateau")
dev.off()


# Plot: taxonomic res ####
# Filter BDMs by taxonomy.ept orders
taxonomy.ept <- taxonomy.bdms[Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera"),]

# Select taxa in BDMf and BDMs that belong to taxonomy.ept orders
n.bdmf.ept <- n.bdmf[names(n.bdmf) %in% taxonomy.ept$Family]
n.bdms.ept <- n.bdms[names(n.bdms) %in% taxonomy.ept$Taxon]

# P2 deviance by resolution ####
# Get BDMf taxonomy.ept deviance
dev.bdmf <- filter(jsdm.bdmf$deviance, Taxon %in% names(n.bdmf.ept))
dev.bdmf$Family <- dev.bdmf$Taxon
dev.bdmf$Rank <- "Family"

# Obtain data for BDMs
dev.bdms <- taxonomy.ept %>%
  left_join(jsdm.bdms$deviance, by = "Taxon") %>%
  filter(Rank != "Family") %>%
  select_(.dots = c(colnames(dev.bdmf), "Family", "Rank"))

plot.data <- bind_rows(dev.bdmf, dev.bdms)
rm(dev.bdmf, dev.bdms)

# Assign occurrence frequency based on whether taxon is in BDMf or not
plot.data$n <- ifelse(plot.data$Trial == 'bdmf', n.bdmf.ept[plot.data$Taxon], n.bdms.ept[plot.data$Taxon])
# Label and order based on family occurrence frequency in BDM families
plot.data$Label <- factor(paste(plot.data$Family, ' - ', n.bdmf.ept[plot.data$Family]), levels = paste(names(sort(n.bdmf.ept, decreasing = TRUE)), ' - ', sort(n.bdmf.ept, decreasing = TRUE)))

# P2 Deviance per family and by rank
g <- ggplot(plot.data, aes(x = Label, y = std.res.dev, color = Rank, size = n, label = Taxon))
g <- g + geom_point(alpha=0.4)
# g <- g + scale_size_continuous(range = c(2,7))
g <- g + scale_radius(range = c(2,8))
g <- g + theme_bw(base_size = 15)
g <- g + theme(axis.text.x = element_text(angle=45,hjust=1))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title="Quality of fit for EPT taxa in BDM",
              y="Standardized deviance",
              x="Family - occurrence frequency",
              size="Occurrence\nfrequency",
              color="Rank")
# g <- ggplotly(g)
print(g)

# plot_ly(data = dev.bdm,
#         x = ~Family, y = ~std.res.dev, color = ~factor(Rank), size = ~n, label=~Taxon,
#         text = ~paste("Family: ", Family, '<br>Std. deviance:', std.res.dev),
#         type = 'scatter')  


# P2 D2 per family and by rank
g <- ggplot(plot.data, aes(x = Label, y = D2, color = Rank, size = n, label = Taxon))
g <- g + geom_point(alpha=0.4)
g <- g + scale_radius(range = c(2,8))
g <- g + theme_bw(base_size = 15)
g <- g + theme(axis.text.x = element_text(angle=45,hjust=1))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title="Explanatory power for EPT taxa in BDM",
              y=expression(paste("D"^2)),
              x="Family - occurrence frequency",
              size="Occurrence\nfrequency",
              color="Rank")
# g <- ggplotly(g)
print(g)


# Plot D2 vs deviance in BDM
plot.data <- dev[Taxon %in% names(n.bdmf.ept) | Taxon %in% names(n.bdms.ept),]

g <- ggplot(plot.data)
g <- g + geom_point(aes(x = std.res.dev, y = D2, color = Trial), size =2)
g <- g + scale_size_continuous(range = c(2,7))
g <- g + theme_bw(base_size = 15)
g <- g + theme(axis.text.x = element_text(angle=45,hjust=1))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title="Explanatory power vs quality of fit for EPT",
              x="Standardized deviance",
              y="D2")
g <- g + scale_color_brewer(palette = "Set1")
print(g)

# P2 d/D2 by trial ####
t <- bind_rows(taxonomy.bdms, taxonomy.bdmf, taxonomy.invf, taxonomy.invfp)
t <- unique(t)

plot.data <- bind_rows(jsdm.bdmf$deviance, jsdm.bdms$deviance, jsdm.invf$deviance, jsdm.invfp$deviance)
plot.data$Trial <- ifelse(plot.data$Trial=="bdmf", "BDMf", plot.data$Trial)
plot.data$Trial <- ifelse(plot.data$Trial=="bdms", "BDMs", plot.data$Trial)
plot.data$Trial <- ifelse(plot.data$Trial=="invf", "CFCH", plot.data$Trial)
plot.data$Trial <- ifelse(plot.data$Trial=="invf_plat", "CFp", plot.data$Trial)
plot.data$Trial <- factor(plot.data$Trial, levels=c("BDMs", "BDMf", "CFCH", "CFp"))

g <- ggplot(data = plot.data)
g <- g + geom_violin(aes(x = Trial, y = D2, fill = Trial))
g <- g + geom_boxplot(aes(x = Trial, y = D2), width=0.20)
g <- g + labs(y = expression(paste("D"^2)),
              fill = "Dataset")
g <- g + geom_hline(yintercept = 0)
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + theme_bw(base_size = 28)
g <- g + theme(plot.title = element_blank(),
               axis.title.x = element_blank())
print(g)

g <- ggplot(data = plot.data)
g <- g + geom_violin(aes(x = Trial, y = std.res.dev, fill = Trial), show.legend = FALSE)
g <- g + geom_boxplot(aes(x = Trial, y = std.res.dev), width=0.20, show.legend = FALSE)
g <- g + labs(y = "Standardized deviance")
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + theme_bw(base_size = 28)
g <- g + theme(plot.title = element_blank(),
               axis.title.x = element_blank())
print(g)



# # P2 D2 vs dev by trial ###
dev.family <- bind_rows(jsdm.bdmf$deviance, jsdm.invf$deviance, jsdm.invfp$deviance)
dev.family$rel.freq[dev.family$Trial=='bdmf'] <- dev.family$n[dev.family$Trial=='bdmf']/nrow(sample.bdmf)
dev.family$rel.freq[dev.family$Trial=='invf'] <- dev.family$n[dev.family$Trial=='invf']/nrow(sample.invf)
dev.family$rel.freq[dev.family$Trial=='invf_plat'] <- dev.family$n[dev.family$Trial=='invf_plat']/nrow(sample.invfp)
dev.family <- dev.family[Taxon %in% names(n.invf),]

g <- ggplot(dev.family)
g <- g + geom_path(lineend="round", linejoin="round", aes(x = std.res.dev, y = D2, group = Taxon, size=rel.freq), alpha=0.3, show.legend=FALSE)
g <- g + geom_point(aes(x = std.res.dev, y = D2, size = rel.freq, color=Trial), alpha=0.5)
# g <- g + geom_density(aes(x = D2, fill=Trial), alpha=0.5)
# g <- g + geom_polygon(aes(x = std.res.dev, y = D2, group=Taxon, size=rel.freq), alpha=0.3)
# g <- g + facet_wrap(~ Trial)
g <- g + scale_radius(range = c(2,7), breaks=seq(0,1,0.2), limits=c(0,1))
g <- g + theme_bw(base_size = 15)
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title="Explanatory power vs quality of fit",
              y=expression(paste("D"^2)),
              # y="D2",
              x="Standardized deviance",
              size="Relative occurrence\nfrequency")
# g <- ggplotly(g)
print(g)

# P2 pred by trial ####
dev.cv <- bind_rows(jsdm.cv.bdm$deviance, jsdm.cv.inv$deviance)
plot.data <- dev.cv %>%
  group_by(Taxon, Type, Trial) %>%
  summarise(mean.std.dev = mean(std.deviance)) %>%
  # spread(Type, mean.std.dev) %>%
  ungroup()

plot.data$n <- NA
plot.data$n <- ifelse(plot.data$Trial=='bdms', n.bdms[plot.data$Taxon], plot.data$n)
plot.data$n <- ifelse(plot.data$Trial=='bdmf', n.bdmf[plot.data$Taxon], plot.data$n)
plot.data$n <- ifelse(plot.data$Trial=='invf', n.invf[plot.data$Taxon], plot.data$n)
plot.data$n <- ifelse(plot.data$Trial=='invf_plat', n.invfp[plot.data$Taxon], plot.data$n)

# g <- ggplot(plot.data, aes(x = Training, y = Testing, color = Trial, size = n))
# g <- g + geom_point(alpha = 0.3)
# g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
# g <- g + theme(strip.background=element_rect(fill="black"), strip.text=element_text(color="white", face="bold"), axis.text=element_text(size=18))
# 
# g <- g + theme_bw(base_size = 18)
# g <- g + guides(color = guide_legend(override.aes = list(size=6)))
# g <- g + labs(x = "Mean standardized deviance for calibration",
#               y = "Mean standardized deviance for prediction",
#               colour = "Trial",
#               size = "Occurrence\nfrequency")
# g <- g + scale_fill_brewer(palette = "Set1")
# g <- g + scale_size_continuous(range = c(2, 8))
# print(g)

g <- ggplot(data = plot.data)
g <- g + geom_boxplot(aes(x = Trial, y = mean.std.dev, fill = Trial))
g <- g + facet_grid(. ~ Type)
g <- g + labs(title="Cross-validation performance by trial",
              subtitle="3-fold cross-validation",
              y = "Mean standardized deviance",
              x = "Trial")
g <- g + theme_bw(base_size = 15)
g <- g + scale_fill_brewer(palette = "Set1")
print(g)

# P2 pred by resolution ####
dev.bdmf <- jsdm.cv.bdm$deviance %>%
  group_by(Type, Trial, Taxon) %>%
  summarise(mean.std.deviance = mean(std.deviance)) %>%
  ungroup() %>%
  filter(Trial == 'bdmf', Taxon %in% names(n.bdmf.ept)) %>%
  mutate(Family = Taxon, Rank = "Family")

dev.bdms <- jsdm.cv.bdm$deviance %>%
  group_by(Type, Trial, Taxon) %>%
  summarise(mean.std.deviance = mean(std.deviance)) %>%
  ungroup() %>%
  filter(Trial == 'bdms', Taxon %in% names(n.bdms.ept)) %>%
  left_join(taxonomy.ept, by = "Taxon")

plot.data <- dev.bdms[, colnames(dev.bdmf)] %>%
  bind_rows(dev.bdmf)
rm(dev.bdmf, dev.bdms)

# For families in BDMf and BDMs, keep only BDMs performance
test <- plot.data %>% 
  filter(Rank=="Family") %>% 
  group_by(Family) %>% 
  summarise(n=n()) %>% 
  filter(n > 2)
test <- filter(plot.data, Taxon %in% test$Family & Trial == 'bdms')
plot.data <- filter(plot.data, !(Taxon %in% unique(test$Family)))
plot.data <- bind_rows(plot.data, test)
plot.data <- filter(plot.data, Family != "Lepidostomatidae")
rm(test)

# Assign occurrence frequency based on whether taxon is in BDMf or not
plot.data$n <- ifelse(plot.data$Trial == 'bdmf', n.bdmf.ept[plot.data$Taxon], n.bdms.ept[plot.data$Taxon])
# Label and order based on family occurrence frequency in BDM families
plot.data$TaxonLabel <- factor(paste(plot.data$Family, ' - ', n.bdmf.ept[plot.data$Family]), levels = paste(names(sort(n.bdmf.ept, decreasing = TRUE)), ' - ', sort(n.bdmf.ept, decreasing = TRUE)))
plot.data$Type <- factor(plot.data$Type, levels=c("Training", "Testing"))

# Facet by Training and Testing
g <- ggplot(data = plot.data)
g <- g + geom_point(aes(x = TaxonLabel, y = mean.std.deviance, color = Rank, size = n), alpha = 0.4)
g <- g + facet_grid(Type ~ .)
g <- g + scale_radius(range = c(2,8))
g <- g + theme_bw(base_size = 15)
g <- g + theme(axis.text.x = element_text(angle=45,hjust=1))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title="Predictive performance of EPT taxa in BDM",
              y="Mean standardized deviance",
              x="Family - occurrence frequency",
              size="Occurrence\nfrequency",
              color="Rank")
print(g)

# Don't facet by training and testing
plot.data <- filter(plot.data, Type=="Testing")

g <- ggplot(data = plot.data)
g <- g + geom_point(aes(x = TaxonLabel, y = mean.std.deviance, color = Rank, size = n), alpha = 0.4)
g <- g + scale_radius(range = c(2,8))
g <- g + theme_bw(base_size = 15)
g <- g + theme(plot.title = element_blank(),
               axis.text.x = element_text(angle=45,hjust=1))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title="Predictive performance of EPT taxa in BDM",
              y="Mean standardized deviance",
              x="Family - occurrence frequency",
              size="Occurrence\nfrequency",
              color="Rank")
print(g)

g <- ggplotly(g)
print(g)

# Invf cross-validation ####
dev.inv <- jsdm.cv.inv$deviance %>%
  group_by(Type, Trial, Taxon) %>%
  summarise(mean.std.dev = mean(std.deviance)) %>%
  ungroup()

g <- ggplot(data = dev.inv)
g <- g + geom_boxplot(aes(x = Trial, y = mean.std.dev, fill = Trial))
g <- g + facet_grid(. ~ Type)
g <- g + labs(title="Cross-validation performance for combined family data",
              subtitle="3-fold cross-validation",
              y = "Mean standardized deviance",
              x = "Dataset")
g <- g + scale_fill_brewer(palette = "Set1")
print(g)

# CH vs plateau predictive deviance ####

# Get predictive performance for Invf CH (filter by plateau)
# Does all CH help predictive performance for combined families?

# Calculate deviance for 'invf' only in the plateau
invf.dev <- jsdm.cv.inv$probability %>%
  filter(SiteId %in% na.omit(env$SiteId[env$BIOGEO=="Mittelland"]), Trial=='invf') %>%
  group_by(Taxon, Fold) %>% # get deviance per taxon in each fold
  summarise(residual.deviance = sum(-2*log(ifelse(Obs == 1, Pred, 1-Pred)), na.rm = TRUE),
         n.samples = length(na.omit(Obs))) %>%
  mutate(std.deviance = residual.deviance/n.samples, Trial = "invf") %>%
  ungroup()

plot.data <- jsdm.cv.inv$deviance %>%
  filter(Trial =='invf_plat') %>%
  select(Taxon, Fold, residual.deviance, n.samples, std.deviance, Trial) %>%
  bind_rows(invf.dev) %>%
  group_by(Taxon, Trial) %>%
  summarise(mean.std.dev = mean(std.deviance)) %>%
  ungroup()
  # mutate(Trial = ifelse(Trial=='invf', 'families CH \n(plateau only)', 'families plateau \n (plateau only)'))


# test <- plot.data %>%
#   spread(Trial, mean.std.dev) %>%
#   mutate(difference = invf-invf_plat) %>%
#   arrange(desc(difference)) %>%
#   select(Taxon, difference)

plot.data <- plot.data %>%
  left_join(test, by = "Taxon") %>%
  arrange(desc(difference))

boxplot(mean.std.dev ~ Trial , data = plot.data, 
        main = "Predictive performance by trial",
        xlab = "Trial (plateau sites only)",
        ylab = "Mean standardized deviance")

# plot.data$Label <- factor(plot.data$Taxon, levels(na.omit(test$Taxon)))

g <- ggplot(plot.data)
g <- g + geom_point(aes(x = reorder(Taxon, -difference), y = mean.std.dev, color = Trial), size = 3)
g <- g + theme_bw(base_size = 15)
# g <- g + theme(axis.text.x = element_text(angle=45,hjust=1))
g <- g + coord_flip()
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + scale_color_brewer(palette = "Set1")
print(g)


# linear.predictor(results) ####
slopes.bdmf <- linear.predictor(jsdm.bdmf)
slopes.bdms <- linear.predictor(jsdm.bdms)
slopes.invf <- linear.predictor(jsdm.invf)
slopes.invfp <- linear.predictor(jsdm.invfp)

slopes.bdmf$Trial <- "BDMf"
slopes.bdms$Trial <- "BDMs"
slopes.invf$Trial <- "CFCH"
slopes.invfp$Trial <- "CFp"

slopes.all <- bind_rows(slopes.bdmf, slopes.bdms, slopes.invf, slopes.invfp)

# Drop taxonomy.ept taxa?
# dt <- dt[Taxon %in% names(n.bdms.taxonomy.ept),]
# dt <- dt[Taxon %in% names(n.bdmf.taxonomy.ept),]

# Get the max and min z-values by variable and by trial
plot.data <- slopes.all %>%
  group_by(Taxon, Trial, Variable) %>%
  summarise(z.min = min(z, na.rm=T), z.max = max(z, na.rm=T)) %>%
  mutate(z.range = z.max-z.min) %>%
  ungroup()

plot.data$Variable <- factor(plot.data$Variable, levels=c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp"))
plot.data$Trial <- factor(plot.data$Trial, levels=c("BDMs", "BDMf", "CFCH", "CFp"))
plot.data$Group <- ifelse(plot.data$Trial=="BDMf" | plot.data$Trial=="BDMs", "BDM", "CF")
g <- ggplot(data = plot.data) 
g <- g + geom_boxplot(aes(x=Group, y=z.range, fill=Trial)) # position = position_dodge2(preserve = "total")
g <- g + facet_grid(. ~ Variable)
g <- g + coord_cartesian(ylim=c(0,19)) # 
g <- g + theme_bw(base_size = 24)
g <- g + theme(axis.title.x=element_blank(),
                axis.text.x=element_text(angle=45, hjust=1)
                # axis.ticks.x=element_blank()
               )
g <- g + scale_fill_brewer(palette = "Set1")

g <- g + labs(y = expression(paste("z"["range"])),
              fill = "Dataset")
print(g)

# Plot the distribution of influence factors over the programs
x <- prepare.inputs(c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp", "Temp2"), sample.invf, center=FALSE)
x$Trial <- "Families CH"
x.bdm <- prepare.inputs(c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp", "Temp2"), sample.bdms, center=FALSE)
x.bdm$Trial <- "BDM"
x.invfp <- prepare.inputs(c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp", "Temp2"), sample.invf.plat, center=FALSE)
x.invfp$Trial <- "Families plateau"
x <- bind_rows(x, x.bdm, x.invfp)
x <- gather(x, Variable, Value, -SiteId, -SampId, -Trial)
x <- filter(x, Variable != "Temp2")
rm(x.bdm, x.invfp)
setDT(x)

g <- ggplot(data = x)
g <- g + geom_violin(aes(x = Trial, y = Value, fill=Trial))
g <- g + geom_boxplot(aes(x = Trial, y = Value), width = 0.15)
g <- g + facet_wrap(~ Variable, scales = "free_x", strip.position="bottom")
g <- g + theme_bw()
g <- g + theme(strip.background = element_blank(), strip.placement = "outside")
g <- g + coord_flip()
g <- g + labs(title = "Distribution of influence factors by trial",
              x = "Trial",
              y = "")
print(g)

# Spatial clusters of impacts ####
dt <- slopes.bdms %>%
  group_by(SampId, Variable) %>%
  summarise(z.range.comm = max(z)-min(z)) %>%
  filter(Variable=='Temp')

dt <- left_join(dt, select(sample.bdms, SiteId, SampId), by = "SampId")
dt <- left_join(dt, inputs$xy, by="SiteId")

g <- ggplot(inputs$ch, aes(POINT_X, POINT_Y))
g <- g + geom_path(lineend = "round")
g <- g + geom_point(data = dt, aes(X, Y, size = z.range.comm), alpha = 0.35)
g <- g + scale_size_continuous(limits = c(min(dt$z.range.comm), max(dt$z.range.comm)), range = c(2, 7))

g <- g + scale_colour_brewer(palette = "Set1")
g <- g + scale_y_continuous(breaks=NULL)
g <- g + scale_x_continuous(breaks=NULL)
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + theme_minimal()
g <- g + labs(main = "",
              x = "",
              y = "")
print(g)

# Site impacts ####
# Get the range of linear predictor per taxon and variable,
# then get the mean z-range per variable for community,
# then get the weight for each influence factor
w.bdms <- slopes.all %>%
  filter(Trial=='BDMs') %>%
  group_by(Taxon, Variable) %>%
  summarise(z.min = min(z, na.rm=T), z.max = max(z, na.rm=T)) %>%
  mutate(z.range = z.max-z.min) %>%
  ungroup() %>%
  group_by(Variable) %>%
  summarise(z.mean.range = mean(z.range)) %>%
  ungroup() %>%
  arrange(desc(z.mean.range)) %>%
  mutate(Weight=z.mean.range/sum(z.mean.range))

w.v.bdms <- w.bdms$Weight
names(w.v.bdms) <- w.bdms$Variable

# Tidy influence factors, normalize them, and add the weight per variable
predictors <- prepare.inputs(K, sample.bdms, center=F) %>%
  gather(Variable, Value, -SiteId, -SampId) %>%
  filter(Variable != "Temp2") %>%
  group_by(Variable) %>%
  mutate(Value = (Value/max(Value))*100) %>%
  ungroup() %>%
  mutate(Weight=w.v.bdms[Variable]) %>%
  na.omit() %>%
  setDT()

predictors <- predictors %>%
  mutate(wv = Value*Weight) %>%
  group_by(SampId) %>%
  summarise(Impact = sum(Value*Weight)) %>%
  left_join(select(sample.bdms, SiteId, SampId), by="SampId") %>%
  select(SiteId, SampId, Impact) %>%
  left_join(inputs$xy, by="SiteId")

write.csv(predictors, "predictors_impact.csv", row.names=F)

# Plot inputs ####
inputs.datasets <- data.table()
datasets <- c("sample.bdmf", "sample.invf", "sample.invfp")
for (d in 1:length(datasets)){
  dataset <- datasets[d]
  K <- c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "Temp")
  dt <- prepare.inputs(k, get(dataset), center = F)
  dt <- gather(dt, Variable, Value, -SiteId, -SampId)
  dt$Trial <- dataset
  inputs.datasets <- bind_rows(inputs.datasets, dt)
}


inputs.datasets$Trial <- ifelse(inputs.datasets$Trial=="sample.bdmf", "BDM", inputs.datasets$Trial)
inputs.datasets$Trial <- ifelse(inputs.datasets$Trial=="sample.invf", "CFCH", inputs.datasets$Trial)
inputs.datasets$Trial <- ifelse(inputs.datasets$Trial=="sample.invfp", "CFp", inputs.datasets$Trial)

inputs.datasets <- filter(inputs.datasets, Variable != "Temp2")
# inputs.datasets$Label <- labeller(inputs.datasets$Variable)

inputs.datasets$Label <- factor(inputs.datasets$Variable, levels = K[K !="Temp2"])
levels(inputs.datasets$Label) <- labeller(levels(inputs.datasets$Label))

# Boxplots
g <- ggplot(inputs.datasets)
g <- g + geom_boxplot(aes(x=Trial, y=Value, fill=Trial), position = position_dodge2(preserve = "total"))
g <- g + facet_wrap(~Label, scales="free", labeller=label_parsed, strip.position="bottom")

g <- g + theme_bw(base_size=20)
g <- g + labs(fill = "Dataset")
g <- g + theme(strip.background = element_blank(), strip.placement = "outside",
               plot.title = element_blank(),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
g <- g + scale_x_discrete(limits=c("CFp", "CFCH", "BDM"))
# g <- g + scale_fill_brewer(palette="Set2")
g <- g + scale_fill_manual(values=c("BDM"="#377eb8", "CFCH"="#4daf4a", "CFp"="#984ea3"))
g <- g + coord_flip()
g

# Density plots
g <- ggplot(inputs.datasets)
g <- g + geom_density(aes(x=Value, fill=Trial), alpha=0.5)
g <- g + facet_wrap(~Label, scales="free", labeller=label_parsed, strip.position="bottom")
g <- g + theme_bw(base_size=20)
g <- g + labs(fill = "Dataset", y = "Density")
g <- g + theme(strip.background = element_blank(), strip.placement = "outside",
               plot.title = element_blank(),
               axis.title.x = element_blank())
g
# Plotting maps of each family with its genera and species is problematic (to code) because there is a wide variation in the number of interesting taxa per family AND the interesting responses per taxon...
# library(gridExtra)
# library(grid)
# library(lattice)

# Titles and/or annotations
# gs <- lapply(1:3, function(ii) 
#   grobTree(rectGrob(gp=gpar(fill=ii, alpha=0.5)), textGrob(ii)))
# 
# lay <- rbind(c(NA,NA,2,2),
#              c(1,1,2,2),
#              c(1,1,3,3),
#              c(NA,NA,3,3))
# 
# select_grobs <- function(lay) {
#   id <- unique(c(t(lay))) 
#   id[!is.na(id)]
# }
# 
# grid.arrange(grobs=gs[select_grobs(lay)], layout_matrix=lay)

# g.baetidae <- map.jsdm.taxon(jsdm.bdmf, "Baetidae", legend=T)
# g.alpinus <- map.jsdm.taxon(jsdm.bdms, "Baetis_alpinus", legend=F)
# g.rhodani <- map.jsdm.taxon(jsdm.bdms, "Baetis_rhodani", legend=F)
# g.muticus <- map.jsdm.taxon(jsdm.bdms, "Baetis_muticus", legend=F)
# 
# gs <- lapply(list(g.baetidae, g.alpinus, g.rhodani, g.muticus), function(g){
#   ggplotGrob(g)
# })
# 
# 
# lay <- rbind(c(1,1,1,2,2),
#              c(1,1,1,2,2),
#              c(3,3,NA,4,4),
#              c(3,3,NA,4,4))
# select_grobs <- function(lay) {
#   id <- unique(c(t(lay)))
#   id[!is.na(id)]
# }

# grid.arrange(grobs=gs[select_grobs(lay)], layout_matrix=lay)
# grid.arrange(g.baetidae, g.alpinus, g.rhodani, g.muticus, nrow=2)

# Grid-arranged maps ####
g.baetidae <- map.jsdm.taxon(jsdm.bdmf, "Baetidae", legend=T)
g.alpinus <- map.jsdm.taxon(jsdm.bdms, "Baetis_alpinus", legend=T)
g.rhodani <- map.jsdm.taxon(jsdm.bdms, "Baetis_rhodani", legend=T)
g.muticus <- map.jsdm.taxon(jsdm.bdms, "Baetis_muticus", legend=T)
ggarrange(g.baetidae, g.alpinus, g.rhodani, g.muticus, ncol=2, nrow=2, common.legend = TRUE, legend="right")

# Calculate how often taxon responses are significant per variable
# A quick hack...
dt <- data.table()
setDT(plot.data)
for (k in 1:length(unique(plot.data$Variable))){
  variable <- unique(plot.data$Variable)[k]
  for (j in 1:length(unique(plot.data$Taxon))){
    taxon <- unique(plot.data$Taxon)[j]
    test <- plot.data[Taxon==taxon & Variable==variable, ]
    
    test <- t(as.matrix(test[, c("beta.mean", "quantile5", "quantile95"), with=F]))
    test <- ifelse(test["beta.mean",][1] < test[c("quantile95"),][2] & test["beta.mean",][1] > test[c("quantile5"),][2], "similar", "different")
    d <- data.table(Taxon=taxon, Variable=variable, Similarity=test)
    dt <- bind_rows(dt, d)
    rm(taxon, test, d)
  }
  rm(variable)
}

test <- dt %>%
  group_by(Variable, Similarity) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

test$Variable <- factor(test$Variable, levels=c(unique(test$Variable)))
g <- ggplot(test)
g <- g + geom_bar(aes(x=Variable, y=count, fill=Similarity), stat="identity")
print(g)

dt.heatmap <- dt
dt.heatmap$n <- n.invf[dt.heatmap$Taxon]
dt.heatmap$Labels <- factor(paste(dt.heatmap$Taxon, ' - ', n.invf[dt.heatmap$Taxon]), levels = paste(names(sort(n.invf)), ' - ', sort(n.invf)))

g <- ggplot(dt.heatmap)
g <- g + geom_tile(aes(x=Variable, y=Labels, fill=Similarity), alpha=0.4)
g <- g + scale_fill_manual(values=c("red", "blue"))
g <- g + scale_x_discrete(position = "top") 
g <- g + labs(title="Significant differences in response between combined families CH vs plateau",
              subtitle="A parameter is statistically similar if the posterior mean is within the 95% CI of its counterpart")

# Parameters ~ taxon. res. ####
beta.samples.bdmf <- extract.beta(jsdm.bdmf)
beta.samples.bdms <- extract.beta(jsdm.bdms)

# BDMf: specify Rank, Family for each taxon, limit to EPT taxa
beta.samples.bdmf <- beta.samples.bdmf %>%
  # group_by(Variable, Taxon) %>%
  # summarise(mean.post = mean(Value), min.post = min(Value), max.post = max(Value)) %>%
  mutate(Rank = "Family", Family = Taxon) %>%
  filter(Taxon %in% names(n.bdmf.ept))

# beta.samples.bdmf$Rank <- "Family"
# beta.samples.bdmf$Family <- beta.samples.bdmf$Taxon
# beta.samples.bdmf <- beta.samples.bdmf[Taxon %in% names(n.bdmf.ept),]

# BDMs: get Rank, Family  for each taxon, limit to EPT taxa
beta.samples.bdms <- beta.samples.bdms %>%
  # group_by(Variable, Taxon) %>%
  # summarise(mean.post = mean(Value), min.post = min(Value), max.post = max(Value)) %>%
  filter(Taxon %in% unique(taxonomy.ept$Taxon)) %>%
  left_join(select(taxonomy.ept, Taxon, Rank, Family), by="Taxon") %>%
  mutate(n=n.bdms.ept[Taxon]) %>%
  filter(n > 56) %>%
  select(-n)

plot.data <- bind_rows(beta.samples.bdmf, beta.samples.bdms)

# Calculate summary statistics for geom_pointrange(x,y,ymin,ymax)
# plot.data <- plot.data %>%
#   group_by(Trial, Variable, Taxon) %>%
#   summarise(mean.post = mean(Value), min.post = min(Value), max.post = max(Value))

plot.data$n.family <- n.bdmf.ept[plot.data$Family]

plot.data$Label <- factor(paste(plot.data$Family, ' - ', plot.data$n.family), levels = paste(names(sort(n.bdmf.ept, decreasing = TRUE)), ' - ', sort(n.bdmf.ept, decreasing = TRUE)))
# plot.data$Label <- factor(paste(plot.data$Family, ' - ', n.bdmf.ept[plot.data$Family]), levels = paste(names(sort(n.bdmf.ept, decreasing = TRUE)), ' - ', sort(n.bdmf.ept, decreasing = TRUE)))


response.bdm.tidy <- gather(response.bdm, Variable, Response, -Taxon)
response.bdm.tidy$Response[response.bdm.tidy$Response==0] <- "Not significant"
response.bdm.tidy$Response[response.bdm.tidy$Response==1] <- "Significant positive"
response.bdm.tidy$Response[response.bdm.tidy$Response==-1] <- "Significant negative"
# test <- filter(plot.data, Variable=='LUD' & Family=="Baetidae")
# test <- filter(plot.data, n.family > 56 & Variable != "Temp2")
test <- filter(plot.data, Variable != "Temp2" & Family == "Baetidae")



test <- left_join(test, response.bdm.tidy, by=c("Variable", "Taxon"))

g <- ggplot(test)
g <- g + geom_violin(aes(x = Taxon, y = Value, fill = Response, group=Taxon), alpha=0.6, trim=TRUE)
# g <- g + geom_pointrange(aes(x=Label, y=mean.post, ymin=min.post, ymax=max.post, color=Rank), alpha=0.5, position="jitter") # position=position_dodge(width=1)
g <- g + geom_hline(yintercept = 0, size=1)
g <- g + facet_grid(Variable ~ ., scales="free_y")
# g <- g + facet_grid(Variable ~ Label, scales="free_y")
# g <- g + facet_grid(. ~ Variable)
g <- g + theme_bw(base_size = 20)
g <- g + theme(plot.title=element_blank(),
               axis.text.x = element_text(angle=25,hjust=1),
               axis.title.x = element_blank(),
               axis.text.y = element_text(size=15))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(y="Parameter",
              color="Response")
g <- g + scale_fill_manual(values=c("Significant positive"="blue", "Not significant"="grey", "Significant negative"="red"))
print(g)

test <- response.bdm.tidy %>% group_by(Variable, Response) %>% summarise(n = n())

# Plot: responses ~ variable ####
g <- ggplot(test)
g <- g + geom_bar(aes(x=Variable,y=n,group=Response, fill=Response), alpha =0.4, stat="identity")
# g <- g + geom_text(stat='count',aes(x=Variable,label=..count.., group=Response),vjust=-1)
g <- g + theme_bw(base_size = 20)
g <- g + theme(plot.title=element_blank(),
               axis.text.x = element_text(angle=45,hjust=1),
               axis.title.x = element_blank(),
               axis.text.y = element_text(size=15))
g <- g + labs(y="Number of responses",
              fill="Response")
g <- g + scale_fill_manual(values=c("Significant positive"="blue", "Not significant"="grey", "Significant negative"="red"))
print(g)

# Plot: occurrence frequency ####
t.bdms <- taxonomy.bdms
t.bdms$Trial <- "BDMs"
t.bdms$n <- n.bdms[t.bdms$Taxon]

t.bdmf <- taxonomy.bdmf
t.bdmf$Trial <- "BDMf"
t.bdmf$n <- n.bdmf[t.bdmf$Taxon]
t.bdm <- bind_rows(t.bdmf, t.bdms)
rm(t.bdmf, t.bdms)

t.bdm <- t.bdm[t.bdm$Rank %in% c("Family", "Genus", "Species"),]
t.bdm <- arrange(t.bdm, Family, desc(n))

t.bdm <- t.bdm %>%
  group_by(Family) %>%
  mutate(n.taxa = n()) %>%
  filter(n.taxa > 2) %>%
  select(-n.taxa) %>%
  mutate(ID = row_number()+1) %>%
  mutate(alpha = n.bdmf[Family]/581)


g <- ggplot(t.bdm, aes(x=ID, y=n))
g <- g + geom_line(aes(group=Family, alpha=alpha), color="grey", show.legend=FALSE)
# g <- g + geom_smooth(aes(group=Family, alpha=alpha), method="lm", formula=(y~sqrt(x)), se=FALSE, color="grey", show.legend=FALSE)
# g <- g + geom_smooth(aes(group=Family, alpha=alpha), method="loess", se=F, color="grey", size=1.25, show.legend=F)
g <- g + geom_point(aes(color=Rank))
g <- g + theme_bw(base_size = 20)
g <- g + theme(plot.title=element_blank())
g <- g + labs(x="Taxon (ordered by occurrence frequency)",
              y="Occurrence frequency")
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
print(g)

# plot.data$Label <- factor(paste(plot.data$Family, ' - ', plot.data$n.family), levels = paste(names(sort(n.bdmf.ept, decreasing = TRUE)), ' - ', sort(n.bdmf.ept, decreasing = TRUE)))

# Given a vector of species in BDM, do a fancy automated plot:
# - get the family and species parameters, responses, rank, and occurrence frequency
# - get the full sample of the marginal posterior taxon-specific parameters
# - do plot parameters with specific order (family first, then species by frequency), formatted labels, and filled significance
# e.g., taxa <- c("Baetis_alpinus", "Baetis_rhodani", "Baetis_muticus", "Baetis_lutheri")
# scale argument controls degree of vertical overlap of density distributions
library(ggridges)

# Loop and plot all families
# Filter EPT families by number of taxa
families <- taxonomy.ept %>% 
  group_by(Family) %>% 
  summarise(species=n()) %>% 
  filter(Family %in% names(n.bdmf.ept)) %>%
  mutate(n.family = n.bdmf[Family]) %>%
  arrange(-n.family)

pdf("P2 parameters by family.pdf", height=9, width=18, onefile = TRUE)
for (j in 1:length(families$Family)){
  family <- families$Family[j]
  taxonomy <- taxonomy.ept[Family == family & Rank != "Family", ]
  
  if (!(family %in% c("Goeridae", "Hydroptilidae", "Lepidostomatidae"))){
    v <- taxonomy$Taxon
    # taxa <- taxa[taxa != "Family"]
    g <- plot.beta.taxa(v, scale = 2, size = 15, legend.position = "right")
    cat("Plot: ", family,"\n")
    print(g)
  }
}
dev.off()


plot.beta.taxa <- function(taxa, scale, size, legend.position){
  n <- c(n.bdmf, n.bdms)
  # reorder the taxa and ID the family
  j <- names(sort(n[taxa], decreasing=TRUE))
  taxonomy <- taxonomy.ept[Taxon %in% j, ]
  family <- unique(taxonomy$Family)
  
  # Extract parameters
  beta.samples.bdmf <- extract.beta(jsdm.bdmf)
  beta.samples.bdms <- extract.beta(jsdm.bdms)
  
  beta.samples.bdmf <- beta.samples.bdmf %>%
    filter(Taxon == family) %>%
    mutate(Rank = "Family", Family = Taxon)
  
  # BDMs: get Rank, Family  for each taxon, limit to EPT taxa
  beta.samples.bdms <- beta.samples.bdms %>%
    filter(Taxon %in% j) %>%
    left_join(select(taxonomy.ept, Taxon, Rank, Family), by="Taxon")
  
  beta.samples <- bind_rows(beta.samples.bdmf, beta.samples.bdms)
  
  # Extract and join significant responses
  response.bdmf <- extract.resp(jsdm.bdmf)
  response.bdmf <- filter(response.bdmf, Taxon != "Goeridae")
  response.bdms <- extract.resp(jsdm.bdms)
  
  response <- bind_rows(response.bdmf, response.bdms)
  response <- unique(response) # eliminate duplicate responses in BDMf and BDMs
  
  responses <- gather(response, Variable, Response, -Taxon)
  responses$Response[responses$Response==0] <- "Not significant"
  responses$Response[responses$Response==1] <- "Significant positive"
  responses$Response[responses$Response==-1] <- "Significant negative"
  
  plot.data <- beta.samples
  plot.data <- left_join(plot.data, responses, by=c("Variable", "Taxon"))
  
  # Labels for taxon, variable, response
  taxonomy$n <- n[taxonomy$Taxon]
  taxonomy <- arrange(taxonomy, n)
  
  taxonomy.taxa <- c(taxonomy$Rank, "Family")
  names(taxonomy.taxa) <- c(taxonomy$Taxon, family)
  taxonomy.taxa <- substring(taxonomy.taxa, 1, 1)
  
  taxon.labels <- paste(names(taxonomy.taxa), " (",taxonomy.taxa," ", n[names(taxonomy.taxa)],")", sep="")
  taxon.labels <- unique(taxon.labels)
  names(taxon.labels) <- names(taxonomy.taxa)
  
  plot.data$TaxonLabel <- taxon.labels[plot.data$Taxon]
  plot.data$TaxonLabel <- factor(plot.data$TaxonLabel, levels=taxon.labels)
  
  # plot.data$ResponseLabel <- factor(plot.data$Response, levels = c("Significant positive", "Not significant", "Significant negative"))
  # rm(beta.samples.bdmf, beta.samples.bdms, response.bdmf, response.bdms, response.bdm, responses)
  plot.data$VariableLabel <- factor(plot.data$Variable, levels = unique(plot.data$Variable))
  levels(plot.data$VariableLabel) <- labeller(levels(plot.data$VariableLabel))
  
  # Plot the data
  g <- ggplot(plot.data)
  g <- g + stat_density_ridges(aes(x = Value, y = TaxonLabel, fill = Response, group=Taxon), scale = scale, alpha = 0.5)
  # quantile_lines = TRUE, quantiles = c(0.05, 0.95)
  
  g <- g + geom_vline(xintercept = 0, alpha=0.7)
  g <- g + facet_grid(. ~ VariableLabel, scales="free", shrink=TRUE, labeller = label_parsed)
  # g <- g + facet_wrap(~ VariableLabel, scales="free_x", ncol=3, labeller = label_parsed)
  
  g <- g + theme_bw(base_size = size)
  g <- g + theme(plot.title=element_blank(),
                 axis.text.x = element_text(angle=35,hjust=1),
                 axis.text.y = element_text(hjust=1),
                 axis.title.y = element_blank(),
                 axis.title.x = element_text(size=24),
                 strip.text = element_text(size = 14),
                 legend.position=legend.position)
  g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
  g <- g + labs(x=expression(paste("Marginal posterior ",beta["jk"]^"taxa")),
                fill="")
  g <- g + scale_fill_manual(values=c("Significant positive"="blue", "Not significant"="grey", "Significant negative"="red"))
  g <- g + scale_y_discrete(expand = c(0.01, 0))
  g
}

pdf('P2 parameters Baetidae.pdf', height=7, width=15)
plot.beta.taxon(c("Baetis_alpinus", "Baetis_lutheri", "Baetis_muticus", "Baetis_rhodani"), scale=2, size=15, legend.position = "top")
dev.off()

# plot.beta.taxon(c("Habroleptoides_confusa", "Habrophlebia_lauta"))
plot.beta.taxon(c("Protonemura_lateralis", "Nemoura_mortoni", "Amphinemura", "Nemoura_minima", "Protonemura_brevistyla", "Protonemura_intricata"))


# IBCH families ####
library(ggridges)
taxa.ibch <- read.csv('IBCH_Taxa.dat', head=F, stringsAsFactors = F, sep=',')
taxa.ibch <- taxa.ibch$V1

# Extract parameters
beta.samples.bdmf <- extract.beta(jsdm.bdmf)
beta.samples.bdmf$Trial <- "BDMf"

beta.samples.invf <- extract.beta(jsdm.invf)
beta.samples.invf$Trial <- "CFCH"

beta.samples.ibch <- bind_rows(beta.samples.bdmf, beta.samples.invf)
rm(beta.samples.bdmf, beta.samples.invf)

beta.samples.ibch <- beta.samples.ibch[Taxon %in% taxa.ibch, ]

# Extract significant responses
response.bdmf <- extract.resp(jsdm.bdmf)
response.bdmf$Trial <- "BDMf"
response.invf <- extract.resp(jsdm.invf)
response.invf$Trial <- "CFCH"

response <- bind_rows(response.bdmf, response.invf)

responses <- gather(response, Variable, Response, -Taxon, -Trial)
responses$Response[responses$Response==0] <- "Not significant"
responses$Response[responses$Response==1] <- "Significant positive"
responses$Response[responses$Response==-1] <- "Significant negative"

plot.data <- beta.samples.ibch
plot.data <- left_join(plot.data, responses, by=c("Variable", "Taxon", "Trial"))
# rm(response.bdmf, response.invf, responses)
setDT(plot.data)

# Format labels
plot.data$TaxonLabel <- factor(plot.data$Taxon, levels=rev(taxa.ibch[taxa.ibch %in% unique(plot.data$Taxon)]))
plot.data$VariableLabel <- factor(plot.data$Variable, levels = K)
levels(plot.data$VariableLabel) <- labeller(levels(plot.data$VariableLabel))
plot.data$Trial <- as.factor(plot.data$Trial)


g <- ggplot(plot.data)
g <- g + stat_density_ridges(aes(x = Value, y = TaxonLabel, linetype=Trial, fill = Response), alpha=0.5, rel_min_height=0.01)
g <- g + geom_vline(xintercept = 0, alpha=0.7)
g <- g + facet_grid(. ~ VariableLabel, scales="free", shrink=TRUE, labeller = label_parsed)

g <- g + theme_bw(base_size = 20)
g <- g + theme(plot.title=element_blank(),
               axis.text.x = element_text(angle=35,hjust=1),
               axis.text.y = element_text(hjust=1),
               axis.title.y = element_blank(),
               axis.title.x = element_text(size=24),
               strip.text = element_text(size = 14),
               legend.position=legend.position)
g <- g + labs(x=expression(paste("Posterior ",beta["jk"]^"taxa")),
              fill="")
g <- g + scale_y_discrete(expand = c(0.01, 0)) # pad the y-axis limits
g <- g + scale_linetype_manual(values=c("BDMf"="dotted", "CFCH"="solid"))
g <- g + scale_fill_manual(values=c("Significant positive"="blue", "Not significant"="grey", "Significant negative"="red"))
# g <- g + scale_fill_manual(values=c("BDMf"="red", "CFCH"="blue"))
# g <- g + scale_size_manual(values=c("BDMf"=1, "CFCH"=2))
g


# Compare responses ####
families <- taxonomy.bdms %>% 
  group_by(Family) %>% 
  summarise(species=n()) %>%
  filter(!is.na(Family) & species > 1) %>%
  arrange(-species)

# Extract parameters
beta.samples.bdmf <- extract.beta(jsdm.bdmf)
setkey(beta.samples.bdmf, Taxon, Variable)
beta.samples.bdms <- extract.beta(jsdm.bdms)
setkey(beta.samples.bdms, Taxon, Variable)

# Extract significant responses
response.bdmf <- extract.resp(jsdm.bdmf)
response.bdms <- extract.resp(jsdm.bdms)

responses.diff <- data.table()
list.F <- lapply(families$Family, function(f){
  # Get the family response and taxonomy
  response <- response.bdmf[Taxon==f,]
  response <- select(response, -Taxon)
  response.family <- as.integer(response[1,])
  names(response.family) <- colnames(response)
  
  # Get the family taxonomy
  taxonomy <- taxonomy.bdms[Family==f,]
  
  # Loop through taxa in the family
  list.J <- lapply(taxonomy$Taxon, function(j){
    
    # Get the taxon response
    response <- response.bdms[Taxon==j,]
    response <- select(response, -Taxon)
    response.taxon <- as.integer(response[1,])
    names(response.taxon) <- colnames(response)
    
    # Loop through influence factors
    list.K <- lapply(K, function(k){
      
      # Get posterior distributions of taxon-specific parameters for family and taxon
      beta.sample.family <- beta.samples.bdmf[Taxon==f & Variable==k,]
      beta.sample.taxon <- beta.samples.bdms[Taxon==j & Variable==k,]
      
      # Perform Kolmogorov-Smirnov test
      # ks <- ks.test(x=beta.sample.family$Value, y=beta.sample.taxon$Value, alternative = "two.sided")
      # Perform Mann-Whitney U test (Wilcoxon rank sum test)
      # dt <- data.table(Family=f, Taxon=j, Variable=k, 
      #                  FamilyResponse=response.family[k], 
      #                  TaxonResponse=response.taxon[k], 
      #                  D=ks$statistic, 
      #                  p.value=ks$p.value)
      
      mwu <- wilcox.test(x=beta.sample.family$Value, y=beta.sample.taxon$Value)
      dt <- data.table(Family=f, Taxon=j, Variable=k,
                       FamilyResponse=response.family[k],
                       TaxonResponse=response.taxon[k],
                       W = mwu$statistic)
    })
    dt.K <- rbindlist(list.K)
  })
  dt.J <- rbindlist(list.J)
  cat("Statistics calculated", match(f, families$Family), "/", length(families$Family),"\n")
  return(dt.J)
})
response.diff <- rbindlist(list.F)
rm(list.F)

# Test p-value vs sample size
n <- 1:1000
n <- n*10
ps <- data.table()
for (i in 1:length(n)){
  sample.size <- n[i]
  bs.family.sample <- sample(beta.sample.family$Value, sample.size, replace=TRUE)
  bs.taxon.sample <- sample(beta.sample.taxon$Value, sample.size, replace=TRUE)
  
  mwu <- wilcox.test(x=bs.family.sample, y=bs.taxon.sample)
  dt <- data.table(p.value = mwu$p.value, statistic = mwu$statistic, sample.size = sample.size)
  ps <- bind_rows(ps, dt)
  cat(sample.size, "\n")
}


# What proportion of families contain species with different significant responses?
test <- response.ks %>%
  group_by(Family, Variable) %>%
  # mutate(n.taxa = uniqueN(Taxon)) %>%
  # Filter significant species responses that are different from their family
  filter(TaxonResponse != FamilyResponse, TaxonResponse != 0) %>%
  summarise(n.diff = n()) %>%
  left_join(families, by=c("Family"))


# NETWORK plots ####
# Experimental code
# library(htmlwidgets)
# JS <- htmlwidgets::JS

# # Example data.tree
# acme <- Node$new("Acme Inc.")
# accounting <- acme$AddChild("Accounting")
# software <- accounting$AddChild("New Software")
# standards <- accounting$AddChild("New Accounting Standards")
# research <- acme$AddChild("Research")
# newProductLine <- research$AddChild("New Product Line")
# newLabs <- research$AddChild("New Labs")
# it <- acme$AddChild("IT")
# outsource <- it$AddChild("Outsource")
# agile <- it$AddChild("Go agile")
# goToR <- it$AddChild("Switch to R")
# 
# print(acme)

test <- bind_rows(response.bdmf[, c("Taxon", "Temp")], response.bdms[, c("Taxon", "Temp")])
test <- unique(test)
testr <- test$Temp; names(testr) <- test$Taxon

test <- left_join(taxonomy.bdms, response.bdms[, c("Taxon", "Temp")], by="Taxon")
# test$response[test$Temp==1] <- "blue"
# test$response[test$Temp==0] <- "grey"
# test$response[test$Temp==-1] <- "red"
test$pathString <- paste5("Animalia", test$Phylum,test$Class,test$Order,test$Family, test$Genus, test$Taxon, sep = "/", na.rm = TRUE)
tree <- as.Node(test)
tree.df <- ToDataFrameNetwork(tree, "name")
test.response <- testr[tree.df$name]

test.response[test.response==1] <- "skyblue"
test.response[test.response==0] <- "grey"
test.response[test.response==-1] <- "red"
test.response[is.na(test.response)] <- "black"

rnw <- ToListExplicit(tree, unname = TRUE) # Convert to data.tree to list-of-list structure

# radialNetwork(rnw)


# Custom colors for the nodes
# https://stackoverflow.com/questions/42568997/r-networkd3-color-node-stroke-for-radialnetwork/42574609#42574609
jsarray <- paste0('["', paste(c("black", test.response), collapse = '", "'), '"]')
nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))

radialNetwork(rnw, 
              linkColour = "#ccc",
              nodeColour = nodeStrokeJS, # #fff
              nodeStroke = nodeStrokeJS,
              textColour = "black",
              opacity = 0.9,
              fontSize = 12)

# diagonalNetwork(rnw,
#                 linkColour = "#ccc",
#                 nodeColour = "#fff",
#                 nodeStroke = nodeStrokeJS,
#                 textColour = "black")

# radialNetwork(acmeNetwork)
# d3ForceNetwork(Links = MisLinks, Nodes = MisNodes,
#                Source = "source", Target = "target",
#                Value = "value", NodeID = "name",
#                Group = "group", width = 550, height = 400,
#                opacity = 0.9, zoom = TRUE)


# Color by BDMf v BDMs ####
dt <- get.taxonomy(sample.bdms)
dt$Dataset[!(dt$Taxon %in% taxonomy.bdmf$Taxon)] <- "BDM Species"

dt$Dataset[dt$Taxon %in% taxonomy.bdmf$Taxon] <- "Both datasets"

test <- filter(taxonomy.bdmf, !(Taxon %in% dt$Taxon))
test$Dataset <- "BDM Family"
dt <- bind_rows(dt, test)
rm(test)

dt$Species[!is.na(dt$Species)] <- labeller.species(dt$Species[!is.na(dt$Species)])
dt$pathString <- paste5("Animalia", dt$Phylum,dt$Class,dt$Order,dt$Family,dt$Genus,dt$Species, sep = "/", na.rm = TRUE)

dt$Taxon[!is.na(dt$Species)] <- labeller.species(dt$Taxon[!is.na(dt$Species)])
taxa <- dt$Dataset; names(taxa) <- dt$Taxon
taxa[taxa=="BDM Species"] <- "blue"
taxa[taxa=="Both datasets"] <- "orange"
taxa[taxa=="BDM Family"] <- "red"

tree <- as.Node(dt)
tree.df <- ToDataFrameNetwork(tree, "name")

tree.df$color <- taxa[tree.df$name]
tree.df$color[is.na(tree.df$color)] <- "grey"

rnw <- ToListExplicit(tree, unname = TRUE)

jsarray <- paste0('["', paste(c("grey", tree.df$color), collapse = '", "'), '"]')
nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))

radialNetwork(rnw, 
              linkColour = "#ccc",
              nodeColour = nodeStrokeJS, # "#fff"
              nodeStroke = nodeStrokeJS,
              textColour = "black",
              opacity = 0.9,
              fontSize = 10,
              height = 1080,
              width = 1080,
              margin=0.1)

# Combine BDMf and BDMs taxonomy and responses to influence factor 'k'
plot.network <- function(k){
  # Combine taxonomy/responses of BDMf and BDMs datasets 
  taxonomy.combined <- bind_rows(taxonomy.bdmf[!(Taxon %in% taxonomy.bdms$Taxon),], taxonomy.bdms)
  response.combined <- bind_rows(response.bdmf[!(Taxon %in% taxonomy.bdms$Taxon),], response.bdms)
  dt <- left_join(taxonomy.combined, response.combined, by="Taxon")
  rm(taxonomy.combined, response.combined)

  # Label formatting for all species
  dt$Species[!is.na(dt$Species)] <- labeller.species(dt$Species[!is.na(dt$Species)])
  dt$Taxon[!is.na(dt$Species)] <- labeller.species(dt$Taxon[!is.na(dt$Species)])

  dt$pathString <- paste5("Animalia", dt$Phylum,dt$Class,dt$Order,dt$Family, dt$Genus, dt$Species, sep = "/", na.rm = TRUE)
  tree <- as.Node(dt, pathDelimiter="/")
  
  # Add responses to the network data frame
  tree.df <- ToDataFrameNetwork(tree, "name")
  
  response <- dt[,c("Taxon", k)]
  colnames(response) <- c("name", "response")
  tree.df <- left_join(tree.df, response, by="name")
  rm(response)
  
  tree.df$color[tree.df$response==1] <- "blue"
  tree.df$color[tree.df$response==0] <- "grey"
  tree.df$color[tree.df$response==-1] <- "red"
  tree.df$color[is.na(tree.df$response)] <- "black"

  # First element of vector colors the root node, which is not contained in data tree conversion to data frame
  jsarray <- paste0('["', paste(c("black", tree.df$color), collapse = '", "'), '"]')
  nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
  
  rnw <- ToListExplicit(tree, unname = TRUE)
  n <- radialNetwork(rnw,
                     # width=900,
                     # height=900,
                     margin=c("top"=0.1, "right"=0.1, "bottom"=0.1, "left"=0.1),
                     linkColour = "#ccc",
                     nodeColour = "#fff",
                     nodeStroke = nodeStrokeJS,
                     textColour = "black",
                     opacity = 0.9,
                     fontSize = 10)
  
  saveNetwork(n, file = paste("radialnw_BDM_",k,".html", sep=''))
  cat("Plot:", k,"\n")
}

for (k in 1:length(K)){
  plot.network(K[k])
}


# Construct a network for BDM species
# Additional features to implement:
# - integrate BDM family responses from response.bdmf
# - size the points according to occurrence frequency (Nodesize	argument:
# character string specifying a column in the Nodes data frame to vary the node radius)
# - functional implementation
# - Shiny implementation with 
test <- left_join(taxonomy.bdms, response.bdms, by = "Taxon")
# test <- filter(test, Order == "Ephemeroptera")
test$pathString <- paste5("Animalia", test$Phylum,test$Class,test$Order,test$Family, test$Genus, test$Taxon, sep = "/", na.rm = TRUE)
tree <- as.Node(test)

# Convert the data.tree to a data.frame
tree.df <- ToDataFrameNetwork(tree, "name")
tree.df$source <- sapply(strsplit(tree.df$from, "/"), function(n){
  tail(n, n=1) # Get the last element from each vector in list
})

tree.df$target <- sapply(strsplit(tree.df$to, "/"), function(n){
  tail(n, n=1) # Get the last element from each vector in list
})


network <- data.table(source = tree.df$source, target = tree.df$target)

nodes <- data.table(name = unique(c(network$source, network$target)))
nodes$id <- 0:(nrow(nodes) - 1)

# create a data frame of the edges that uses id 0:9 instead of their names
edges <- network %>%
  left_join(nodes, by = c("source" = "name")) %>%
  select(-source) %>%
  rename(source = id) %>%
  left_join(nodes, by = c("target" = "name")) %>%
  select(-target) %>%
  rename(target = id)

# nodes$group <- 1
edges$width <- 1

# Color by response (0 = negative, 1 = not significant, 2 = positive, 3 = not modelled)
nodes <- nodes %>%
  left_join(select(response.bdms, Taxon, Temp), by = c("name" = "Taxon"))
nodes$group <- nodes$Temp; nodes$Temp <- NULL

nodes <- nodes %>%
  left_join(select(response.bdmf, Taxon, Temp), by = c("name" = "Taxon"))

nodes$group <- ifelse(is.na(nodes$Temp), nodes$group, nodes$Temp)
nodes$group <- nodes$group + 1
nodes$group[is.na(nodes$group)] <- 3 # replace these with BDM family responses

nodes$group[nodes$group==0] <- "significant negative"
nodes$group[nodes$group==1] <- "not significant"
nodes$group[nodes$group==2] <- "significant positive"
nodes$group[nodes$group==3] <- "taxonomy node"


ColourScale <- 'd3.scaleOrdinal()
            .domain(["significant negative", "not significant", "significant positive", "taxonomy node"])
           .range(["#D93434", "#B0B0B0", "#4545C7", "#A960D1"]);'

nodes$size <- (n.bdmf[nodes$name]/max(n.bdmf)*10)^2
nodes$size <- (n.bdms[nodes$name]/max(n.bdmf)*10)^2
nodes$size[is.na(nodes$size)] <- 1

forceNetwork(Links = edges, Nodes = nodes, 
             Source = "source",
             Target = "target",
             NodeID ="name",
             Group = "group",
             Value = "width",
             Nodesize = "size",
             opacity = 0.9,
             fontSize = 22,
             zoom = TRUE,
             legend = TRUE,
             charge = -60, # positive or negative node "attraction" to each other
             # linkDistance = networkD3::JS("function(d) { return 5*d.value; }"),
             colourScale = JS(ColourScale))

# Construct links
# tree.df$source <- unlist(strsplit(tree.df$from, "/"))[1]
# tree.df$target <- sapply(strsplit(tree.df$to, "/"), function(n){
#   tail(n, n=1) # Get the last element from each vector in list
# })

# Construct simple network
simpleNetwork(tree.df)


# # Construct a more complex network from BDM
# # links <- select(tree.df, from, to)
# 
# # links$source <- as.numeric(as.factor(links$from))
# # links$target <- as.numeric(as.factor(links$to))
# 
# # Convert from-to linkages from factor (altogether) to numeric
# taxa <- c(tree.df$from, tree.df$to)
# x <- as.numeric(as.factor(c(tree.df$from, tree.df$to)))-1
# links <- data.table(source = x[1:(length(x)/2)], target = x[((length(x)/2)+1):length(x)])
# links$source <- links$source - 1
# links$target <- links$target - 1
# 
# nodes <- data.table(name = unique(taxa), group = 1)
# # nodes <- data.frame(name = as.numeric(as.factor(test$Taxon)), group = 1)
# # nodes$name <- nodes$name-1
# 
# forceNetwork(Links=links,Nodes = nodes,
#              Source = 'source', Target = 'target', 
#              NodeID = 'name', Group = 'group', zoom = TRUE)

# # Advanced network example
# # Load data
# data(MisLinks)
# data(MisNodes)
# 
# # Plot
# forceNetwork(Links = MisLinks, Nodes = MisNodes,
#              Source = "source", Target = "target",
#              Value = "value", NodeID = "name",
#              Group = "group", opacity = 0.8)

# # Example from: https://stackoverflow.com/questions/45179424/r-networkd3-package-forcenetwork-coloring
# # Load package
# library(networkD3)
# library(dplyr) # to make the joins easier
# 
# # Create fake data
# source <- c("A", "A", "A", "A",
#             "B", "B", "C", "C", "D")
# target <- c("B", "C", "D", "J",
#             "E", "F", "G", "H", "I")
# networkData <- data.frame(source, target, stringsAsFactors = FALSE)
# 
# nodes <- data.frame(name = unique(c(source, target)), stringsAsFactors = FALSE)
# nodes$id <- 0:(nrow(nodes) - 1)
# 
# 
# # create a data frame of the edges that uses id 0:9 instead of their names
# edges <- networkData %>%
#   left_join(nodes, by = c("source" = "name")) %>%
#   select(-source) %>%
#   rename(source = id) %>%
#   left_join(nodes, by = c("target" = "name")) %>%
#   select(-target) %>%
#   rename(target = id)
# 
# edges$width <- 1
# 
# # make a grouping variable that will match to colours
# nodes$group <- ifelse(nodes$name %in% source, "lions", "tigers")
# 
# 
# ColourScale <- 'd3.scaleOrdinal()
# .domain(["lions", "tigers"])
# .range(["#FF6900", "#694489"]);'
# 
# forceNetwork(Links = edges, Nodes = nodes, 
#              Source = "source",
#              Target = "target",
#              NodeID ="name",
#              Group = "group",
#              Value = "width",
#              opacity = 0.9,
#              zoom = TRUE,
#              colourScale = JS(ColourScale))

# P2 dotplot BDMf vs Invf ####
beta.taxa.bdmf <- extract.beta(jsdm.bdmf)
beta.taxa.bdmf <- left_join(beta.taxa.bdmf, taxonomy.bdmf, by="Taxon")
beta.taxa.invf <- extract.beta(jsdm.invf)
beta.taxa.invf <- left_join(beta.taxa.invf, taxonomy.invf, by="Taxon")

beta.taxa.bdmf <- beta.taxa.bdmf %>%
  group_by(Trial, Taxon, Variable) %>%
  summarise(beta.mean = mean(Value), quantile5 = quantile(Value, probs=0.05), quantile95 = quantile(Value, probs=0.95)) %>%
  ungroup() %>%
  left_join(taxonomy.bdmf, by="Taxon") %>%
  mutate(n = n.bdmf[Taxon]/nrow(sample.bdmf))

beta.taxa.invf <- beta.taxa.invf %>%
  group_by(Trial, Taxon, Variable) %>%
  summarise(beta.mean = mean(Value), quantile5 = quantile(Value, probs=0.05), quantile95 = quantile(Value, probs=0.95)) %>%
  ungroup() %>%
  left_join(taxonomy.invf, by="Taxon") %>%
  mutate(n = n.invf[Taxon]/nrow(sample.invf))

beta.taxa.families <- bind_rows(beta.taxa.bdmf, beta.taxa.invf)
rm(beta.taxa.bdmf, beta.taxa.invf)

plot.data <- filter(beta.taxa.families, Rank == 'Family', Taxon %in% intersect(names(n.bdmf), names(n.invf)))
plot.data$Labels <- factor(paste(plot.data$Taxon, ' - ', n.invf[plot.data$Taxon]), levels = paste(names(sort(n.invf)), ' - ', sort(n.invf)))

g <- ggplot(plot.data)
g <- g + geom_pointrange(aes(x=Labels, y=beta.mean, ymin=quantile5, ymax=quantile95, color=Trial), alpha=0.4)
g <- g + geom_hline(yintercept=0)
g <- g + facet_grid(.~Variable, scales="free")
g <- g + theme_bw(base_size=12)
g <- g + theme(axis.text.x = element_text(angle=45,hjust=1))
g <- g + coord_flip()
g <- g + guides(colour = guide_legend(override.aes = list(size=3))) # guides(colour=FALSE)
g <- g + scale_color_brewer(palette = "Set1")
g <- g + labs(title="Mean of marginal posterior taxon-specific parameter samples (BDM vs combined families)")
pdf('P2 response dotplot invf-bdmf.pdf', height=12, width=14)
print(g)
dev.off()

# P2 dotplot CH vs plateau ####
beta.taxa.invf <- extract.beta(jsdm.invf)
beta.taxa.invf <- left_join(beta.taxa.invf, taxonomy.invf, by="Taxon")

beta.taxa.invfp <- extract.beta(jsdm.invfp)
beta.taxa.invfp <- left_join(beta.taxa.invfp, taxonomy.invfp, by="Taxon")

beta.taxa.invf <- beta.taxa.invf %>%
  group_by(Trial, Taxon, Variable) %>%
  summarise(beta.mean = mean(Value), quantile5 = quantile(Value, probs=0.05), quantile95 = quantile(Value, probs=0.95)) %>%
  ungroup() %>%
  left_join(taxonomy.invf, by="Taxon") %>%
  mutate(n = n.invf[Taxon]/nrow(sample.invf))

beta.taxa.invfp <- beta.taxa.invfp %>%
  group_by(Trial, Taxon, Variable) %>%
  summarise(beta.mean = mean(Value), quantile5 = quantile(Value, probs=0.05), quantile95 = quantile(Value, probs=0.95)) %>%
  ungroup() %>%
  left_join(taxonomy.invfp, by="Taxon") %>%
  mutate(n = n.invfp[Taxon]/nrow(sample.invfp))

beta.taxa.families.invfp <- bind_rows(beta.taxa.invf, beta.taxa.invfp)
rm(beta.taxa.invf, beta.taxa.invfp)

plot.data <- filter(beta.taxa.families.invfp, Rank == 'Family', Taxon %in% intersect(names(n.invf), names(n.invf)))
plot.data$Labels <- factor(paste(plot.data$Taxon, ' - ', n.invf[plot.data$Taxon]), levels = paste(names(sort(n.invf)), ' - ', sort(n.invf)))

g <- ggplot(plot.data)
g <- g + geom_pointrange(aes(x=Labels, y=beta.mean, ymin=quantile5, ymax=quantile95, color=Trial), alpha=0.4)
g <- g + geom_hline(yintercept=0)
g <- g + facet_grid(.~Variable, scales="free")
g <- g + theme_bw(base_size=12)
g <- g + theme(axis.text.x = element_text(angle=45,hjust=1))
g <- g + coord_flip()
g <- g + guides(colour = guide_legend(override.aes = list(size=3))) # guides(colour=FALSE)
g <- g + scale_color_brewer(palette = "Set1")
g <- g + labs(title="Mean of marginal posterior taxon-specific parameter samples (combined families CH vs plateau)")
pdf('P2 response dotplot ch-plateau.pdf', height=12, width=14)
print(g)
dev.off()

# plot.data <- dev
# plot.data <- filter(plot.data, Trial != "BDMs")
# plot.data$n <- n.invf[plot.data$Taxon]
# plot.data <- na.omit(plot.data)
# plot.data$TaxonLabel <- factor(paste(plot.data$Taxon, ' - ', n.invf[plot.data$Taxon]), levels = paste(names(sort(n.invf, decreasing = TRUE)), ' - ', sort(n.invf, decreasing = TRUE)))
# 
# test <- plot.data %>%
#   select(Taxon, Trial, D2, n, TaxonLabel) %>%
#   spread(Trial, D2) %>%
#   na.omit %>%
#   mutate(BDMf.CFCH = BDMf - CFCH,
#          BDMf.CFp = BDMf - CFp,
#          CFCH.CFp = CFCH - CFp) %>%
#   select(Taxon, BDMf.CFCH, BDMf.CFp, CFCH.CFp) %>%
#   gather(Variable, Value, -Taxon)
# 
# test$TaxonLabel <- factor(paste(test$Taxon, ' - ', n.invf[test$Taxon]), levels = paste(names(sort(n.invf, decreasing = TRUE)), ' - ', sort(n.invf, decreasing = TRUE)))
# 
# g <- ggplot(test)
# g <- g + geom_col(aes(x=TaxonLabel, y=Value))
# g <- g + facet_grid(Variable ~ .)
# print(g)
# 
# 
# 
# # Families of interest
# families <- c("Limnephilidae", "Heptageniidae", "Nemouridae", "Baetidae", "Rhyacophilidae", "Leuctridae", "Leptophlebiidae")
# dt <- bdm.fit[bdm.fit$Family %in% families, ]
# dt <- arrange(dt, Family, -n)
# 
# dt %>% group_by(Family) %>% summarise(n=n()) %>% arrange(-n)
# setDT(dt)
# 
# # Construct the inputs
# inputs.bdm <- prepare.inputs(K, sample.bdms, center = FALSE)
# inputs.bdm <- gather(inputs.bdm, Variable, Value, -SiteId, -SampId)
# inputs.bdm <- filter(inputs.bdm, Variable != "Temp2")
# 
# # labeller() returns expressions for nice plot labels for each influence factor
# inputs.bdm$Label <- factor(inputs.bdm$Variable, levels = K[K !="Temp2"])
# levels(inputs.bdm$Label) <- labeller(levels(inputs.bdm$Label))
# 
# # Construct the probabilities from BDM families and species
# prob.bdmf <- jsdm.bdmf$probability[Taxon %in% families,]
# prob.bdms <- jsdm.bdms$probability[Taxon %in% dt$Taxon,]
# prob.bdm <- bind_rows(prob.bdmf, prob.bdms)
# 
# # Join the inputs to the probabilities
# prob.bdm <- left_join(prob.bdm, inputs.bdm, by=c("SiteId", "SampId"))
# prob.bdm$Obs <- as.factor(prob.bdm$Obs)
# setDT(prob.bdm)
# rm(prob.bdmf, prob.bdms, inputs.bdm)
# 
# resp.bdmf <- response.bdmf[Taxon %in% families, ]
# resp.bdms <- response.bdms[Taxon %in% dt$Taxon, ]
# resp.bdm <- bind_rows(resp.bdmf, resp.bdms)
# resp.bdm <- gather(resp.bdm, Variable, Response, -Taxon)
# setDT(resp.bdm)
# resp.bdm <- resp.bdm[Response != 0, ]
# rm(resp.bdmf, resp.bdms)
# 
# # # Prepare taxa/D2 statistics
# # results$deviance <- arrange(results$deviance, -n)
# # taxa <- results$deviance$Taxon
# # d <- results$deviance$D2; names(d) <- taxa
# 
# for (f in 1:length(families)){
#   family <- families[f]
#   taxa <- dt[Family == family, ]
# 
#   J <- taxa$Taxon
# 
#   for (j in 1:length(J)){
#     taxon <- J[j]
#     responses <- resp.bdm[Taxon==taxon, ]
#     k <- unique(responses[["Variable"]])
#     if(any(c("Temp", "Temp2") %in% k)){
#       k <- k[!(k %in% c("Temp", "Temp2"))]
#       k <- c(k, "Temp")
#     }
#     k <- k[k != "Temp2"]
#     par(mar=c(0.2,0.2,0.2,0.2))
#     pdf(paste('significant_responses/', taxon, '.pdf', sep=''), height = 6, width = 8, onefile = TRUE)
#     if (length(k)>0){
#       for (i in 1:length(k)){
#         variable <- k[i]
#         plot.data <- prob.bdm[Taxon==taxon & Variable==variable,]
# 
#         g <- ggplot(data = plot.data, aes(x = Value, y = Pred, color = Obs))
#         g <- g + geom_point(size=2, alpha = 0.25)
#         g <- g + facet_wrap(~ Label, scales = "free_x", labeller=label_parsed, strip.position="bottom")
#         g <- g + theme_bw(base_size=14)
#         g <- g + labs(title = "",
#                       x = "",
#                       y = "Probability of occurrence",
#                       color = "Observation")
#         g <- g + scale_color_manual(name = "Observation", values=c("#FF0000", "#0077FF"), labels=c("Absence", "Presence"))
#         g <- g + guides(color = FALSE)
#         g <- g + scale_y_continuous(limits = c(0,1), breaks=(seq(0,1,0.2)))
#         g <- g + theme(strip.background = element_blank(), strip.placement = "outside",
#                        plot.title = element_blank(),
#                        axis.text.y = element_text(size=18),
#                        axis.title.x = element_blank(),
#                        axis.text.x = element_text(size=18),
#                        strip.text = element_text(size=24))
#         cat("Plotting: ", taxon, "k:", variable, "\n")
#         print(g)
#       }
#     }
#     dev.off()
#   }
# }
# 
# response.bdm <- bind_rows(response.bdmf, response.bdms) # bind together, since none of the families of interest are duplicates
# response.bdm <- unique(response.bdm) # remove duplicate rows
# dt <- left_join(dt, response.bdm, by="Taxon")


