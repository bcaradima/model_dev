### RESULTS: Paper 1 ####
# Variable selection ####
# K <- c("A10m", "A100m", "A1km", "A.EDO", "IAR", "TU.Dm", "TU.Cr", "LUD", "Urban", "UI", "F10m", "F100m", "F1km", "F.EDO","Forest", "FRI", "Temp", "FV", "WV", "BM", "Morph", "WW", "HP", "Temp2")
# 
# predictors <- prepare.inputs(K, sample.bdms, center = TRUE)
# 
# cor.plot(select(predictors, -SiteId, -SampId), numbers=TRUE)
# write.csv(predictors, 'outputs/predictors.csv', row.names=FALSE)
# 
# 
# predictors.trans <- prepare.inputs(K, sample.bdms, center = FALSE)
# 
# predictors.trans <- mutate(predictors.trans, 
#                      A10m=log10(A10m+1),
#                      A100m=log10(A100m+1),
#                      A1km=log10(A1km+1),
#                      A.EDO=log10(A.EDO+1),
#                      IAR=log10(IAR+1),
#                      TU.Dm=log10(TU.Dm+1),
#                      TU.Cr=log10(TU.Cr+1),
#                      LUD=log10(LUD+1),
#                      Urban=log10(Urban+1),
#                      UI=log10(UI+1),
#                      FRI=log10(FRI+1)
#                      )
# 
# predictors.trans <- rename(predictors.trans,
#                            "log10.A10m" = "A10m",
#                            "log10.A100m" = "A100m",
#                            "log10.A1km" = "A1km",
#                            "log10.A.EDO" = "A.EDO",
#                            "log10.IAR" = "IAR",
#                            "log10.TU.Dm" = "TU.Dm",
#                            "log10.TU.Cr" = "TU.Cr",
#                            "log10.LUD" = "LUD",
#                            "log10.Urban" = "Urban",
#                            "log10.UI" = "UI",
#                            "log10.FRI" = "FRI"
#                            )
# 
# write.csv(predictors.trans, 'outputs/predictors_transformed.csv', row.names=FALSE)
# Variable selection ####
vs.bdms <- extract.vsr("variable_selection/bdms_untransformed", sample.bdms)

# Get mean relative deviance for each model over entire community (excluding rare taxa)
vs.bdms.mean <- vs.bdms %>%
  filter(!is.infinite(rdev.test) & n > 56) %>%
  group_by(Model) %>%
  summarise(mrd.train = mean(rdev.train, na.rm=TRUE), mrd.test = mean(rdev.test, na.rm=TRUE)) %>%
  arrange(mrd.test) %>%
  mutate(Parameters = vapply(strsplit(Model, " "), length, integer(1)))

plot.data <- vs.bdms.mean %>%
  arrange(Parameters) %>%
  group_by(Parameters) %>%
  mutate(Label = paste(Parameters, " (", formatC(uniqueN(Model), big.mark=",")," models)", sep="")) %>%
  ungroup() %>%
  mutate(Label = factor(Label, levels = unique(Label)))

g <- ggplot(plot.data, aes(x = mrd.train, y = mrd.test, color = as.factor(Label)))
g <- g + geom_point(alpha=0.35)
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25)
g <- g + theme_bw(base_size = 18)
# g <- g + theme_gray(base_size = 17) # gray improves contrast
g <- g + labs(y = "Mean standardized deviance during prediction",
              x = "Mean standardized deviance during calibration",
              color = "Number of\nparameters")
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + scale_colour_brewer(palette = "Set1")
print(g)

# PREDICTORS ####
K <- c("Temp", "Temp2", "FV", "F100m", "LUD", "IAR")
predictors <- prepare.inputs(K, sample.bdms, center = TRUE)
write.csv(predictors, 'outputs/predictors.csv', row.names = FALSE)

cor(predictors[,K])
write.csv(cor(predictors[,K]), 'outputs/paper 1 extensions/correlations.csv')
### INDIVIDUAL MODELS ####
# Description: this script selects influence factors, fits and validates the individual models, and outputs equivalent joint model workspaces using the deploy.jsdm() function. Most of these operations are purely for paper 1.

# > Select influence factors ####
# Predictors for paper 1 (BDM species)
# K <- c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD")
# predictors <- prepare.inputs(K, sample.bdms, center=TRUE)
# write.csv(predictors, "outputs/predictors_p1.csv", row.names = F)

# P2 variable selection ####
# Predictors for variable selection per dataset (paper 2)
# predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.bdmf, center = TRUE)
# write.csv(predictors.vs, 'outputs/predictors_vs_bdmf.csv', row.names = F)
# 
# predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.bdms, center = TRUE)
# write.csv(predictors.vs, 'outputs/predictors_vs_bdms.csv', row.names = F)
# 
# predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.invf, center = TRUE)
# write.csv(predictors.vs, 'outputs/predictors_vs_invf.csv', row.names = F)
# 
# predictors.vs <- prepare.inputs(c("A10m","A100m","A1km","A.EDO","IAR","TU.Dm","TU.Cr","LUD","Urban","UI","F10m","F100m","F1km", "p_forest", "F.EDO","bFRI","FRI","Temp","FV","WV","BM","Morph","WW","HP","Temp2","Noise"), sample.invf.plat, center = TRUE)
# write.csv(predictors.vs, 'outputs/predictors_vs_invfp.csv', row.names = F)

# > Individual model calibration ####
isdm <- run.isdm(sample.bdms, predictors, trace = FALSE)

# How many models lowered the deviance?
table(sapply(isdm$models, function(m){
  m$deviance < m$null.deviance
}))

# How many taxa used optim?
table(sapply(isdm$models, function(m){
  m$optim
}))

### > k-fold cross-validation ####
cv.iSDM <- cv.isdm(predictors)
# Transformations ####
# Description: this script prepares the results for paper 1. It processes the joint model workspaces and combines the results with those of the individual models before producing various plots.

# JOINT MODELS ####
# deploy.jsdm() ####
# # K <- c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD")
# deploy.jsdm(K, sample.bdms, "bdms", center=T, cv=0)
# deploy.jsdm(K, sample.bdms, "bdms_train1", center=T, cv=1)
# deploy.jsdm(K, sample.bdms, "bdms_train2", center=T, cv=2)
# deploy.jsdm(K, sample.bdms, "bdms_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "FF0", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "FF0_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "FF0_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "FF0_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TF0", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TF0_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TF0_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TF0_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT0", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT0_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT0_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT0_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT1", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT1_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT1_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT1_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT2", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT2_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT2_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT2_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT3", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT3_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT3_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT3_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT4", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT4_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT4_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT4_train3", center=T, cv=3)


fit.FF0 <- extract.jsdm('outputs/paper 1 extensions', 'FF0') # simplest joint model
fit.TF0 <- extract.jsdm('outputs/paper 1 extensions', 'TF0') # correlated beta.comm
fit.TT0 <- extract.jsdm('outputs/paper 1 extensions', 'TT0') # correlated beta.comm, site effects
fit.TT1 <- extract.jsdm('outputs/paper 1 extensions', 'TT1')
fit.TT2 <- extract.jsdm('outputs/paper 1 extensions', 'TT2')
# fit.TT3 <- extract.jsdm('outputs/paper 1 extensions', 'TT3')
# fit.TT4 <- extract.jsdm('outputs/paper 1 extensions', 'TT4')

### > Combine i/jSDM results ####
# Warning: manually check columns being dropped!
# probability <- bind_rows(isdm.prob, jsdm$probability[, colnames(isdm.prob)])

# Deviance (fit) ####
deviance <- bind_rows(isdm$deviance, fit.FF0$deviance, fit.TF0$deviance, fit.TT0$deviance, fit.TT1$deviance, fit.TT2$deviance)
#jsdm.TT3$deviance, jsdm.TT4$deviance

deviance$Model <- ifelse(deviance$Model=="FF0", "UF0", deviance$Model)
deviance$Model <- ifelse(deviance$Model=="FT0", "UT0", deviance$Model) # random effect model
deviance$Model <- ifelse(deviance$Model=="TF0", "CF0", deviance$Model)
deviance$Model <- ifelse(deviance$Model=="TT0", "CT0", deviance$Model)

deviance$Model <- ifelse(deviance$Model=="TT1", "CT1", deviance$Model)
deviance$Model <- ifelse(deviance$Model=="TT2", "CT2", deviance$Model)


# Calculate 
plot.data <- deviance %>%
  mutate(prevalence = n.present / n.samples * 100,
         ModelLabel = factor(Model, levels=c("Null model", "iSDM", "UF0", "UT0", "CF0", "CT0", "CT1", "CT2")),
         TaxonLabel = factor(Taxon, levels=c(names(sort(n.bdms, decreasing=FALSE)))))

null.deviance <- plot.data %>%
  filter(Model=="iSDM") %>%
  mutate(std.deviance = null.deviance/n.samples,
         Model = "Null model",
         ModelLabel = factor(Model, levels=c("Null model", "iSDM", "UF0", "UT0", "CF0", "CT0", "CT1", "CT2"))
         )


# Construct data for connecting iSDM and mSDM points
data1 <- plot.data %>%
  filter(Model %in% c("iSDM", "UF0")) %>%
  select(Taxon, Model, prevalence) %>%
  spread(Model, prevalence) %>%
  rename("Taxon"="Taxon", "x1"="iSDM", "x2"="UF0")

data2 <- plot.data %>%
  filter(Model %in% c("iSDM", "UF0")) %>%
  select(Taxon, Model, std.deviance) %>%
  spread(Model, std.deviance) %>%
  rename("Taxon"="Taxon", "y1"="iSDM", "y2"="UF0")

dev.segments <- left_join(data1, data2, by = "Taxon")

# Construct data for analytical solution of null deviance 
n <- 580
m <- 1:(n-1)

null.dev.line <- tibble(null.deviance = -2*(m/n*log(m/n)+(n-m)/n*log(1-m/n)), n.samples = m)
null.dev.line$prevalence <- (m/n)*100
null.dev.line$Model <- "Null model"

plot.data <- plot.data %>%
  filter(Model %in% c("iSDM", "UF0"))

g1 <- ggplot()
g1 <- g1 + geom_point(data = plot.data, aes(x=prevalence, y=std.deviance, color=Model), size=3, alpha=0.4)
g1 <- g1 + geom_point(data = null.deviance, aes(x=prevalence, y=std.deviance, color=Model), size = 3, alpha=0.4)
g1 <- g1 + geom_line(data = null.dev.line, aes(x=prevalence, y=null.deviance, color=Model), linetype = "dashed", alpha=0.4, show.legend = FALSE)
g1 <- g1 + geom_segment(data = dev.segments, aes(x=x1, y=y1, xend=x2, yend=y2), alpha = 0.3)

g1 <- g1 + scale_colour_manual(values=c("UF0" = "#048504", "iSDM" = "#790FBF", "Null model" = "#000000"))
g1 <- g1 + theme_bw(base_size = 20)
g1 <- g1 + theme(axis.text=element_text(size=14),
                 plot.title = element_blank())
g1 <- g1 + labs(y = "Standardized deviance",
                x = "Prevalence (%)")
g1 <- g1 + guides(colour = guide_legend(override.aes = list(size=6)), shape = FALSE, color = FALSE)

rm(data1, data2, dev.segments, n, m)

# Plot P1 - D2 vs std dev ###
# Point segments
data1 <- spread(select(plot.data, Taxon, Model, D2), Model, D2) # y coordinates
colnames(data1) <- c("Taxon", "y1", "y2")
data2 <- spread(select(plot.data, Taxon, Model, std.deviance), Model, std.deviance) # x coordinates
colnames(data2) <- c("Taxon", "x1", "x2")
dev.segments <- left_join(data1, data2, by = "Taxon")

g2 <- ggplot()
g2 <- g2 + geom_point(data = plot.data, aes(x = std.deviance, y = D2, color = Model, size = n.present), alpha = 0.4)
g2 <- g2 + geom_segment(data = dev.segments, aes(x = x1, y = y1, xend = x2, yend= y2), alpha = 0.2)
g2 <- g2 + scale_colour_manual(values=c("UF0" = "#048504", "iSDM" = "#790FBF"))
g2 <- g2 + theme_bw(base_size = 20)
g2 <- g2 + theme(axis.text=element_text(size=14),
                 plot.title = element_blank())
g2 <- g2 + labs(x = "Standardized deviance",
                y = expression("D"["j"]^2),
                colour = "Model", 
                size = "Number of\noccurrences")
g2 <- g2 + guides(colour = guide_legend(override.aes = list(size=6)))

rm(data1, data2, dev.segments)

plot.data <- deviance %>%
  mutate(prevalence = n.present / n.samples * 100,
         ModelLabel = factor(Model, levels=c("Null model", "iSDM", "UF0", "UT0", "CF0", "CT0", "CT1", "CT2")),
         ModelLabel = factor(Model, levels=rev(levels(ModelLabel))),
         TaxonLabel = factor(Taxon, levels=c(names(sort(n.bdms, decreasing=FALSE))))
         )

g3 <- ggplot(data = plot.data)
g3 <- g3 + geom_boxplot(aes(x=ModelLabel, y=std.deviance))
g3 <- g3 + coord_flip()
g3 <- g3 + theme_bw(base_size = 20)
g3 <- g3 + labs(x="Model",
                y="Standardized deviance")

# Points and LOESS curves
plot.data <- deviance %>%
  mutate(prevalence = n.present / n.samples * 100,
         ModelLabel = factor(Model, levels=c("Null model", "iSDM", "UF0", "UT0", "CF0", "CT0", "CT1", "CT2")),
         TaxonLabel = factor(Taxon, levels=c(names(sort(n.bdms, decreasing=FALSE))))
  )

g4 <- ggplot(data = plot.data)
g4 <- g4 + geom_point(aes(x=prevalence, y=std.deviance, color=ModelLabel), alpha=0.4)
g4 <- g4 + stat_smooth(aes(x=prevalence, y=std.deviance, group=ModelLabel, color=ModelLabel), geom='line', size = 2, alpha=0.50, se=TRUE, method = "loess")
g4 <- g4 + geom_point(data = null.deviance, aes(x=prevalence, y=std.deviance, color=ModelLabel), alpha=0.4)
g4 <- g4 + geom_line(data = null.dev.line, aes(x=prevalence, y=null.deviance), color="black", linetype = "dashed", alpha=0.4, show.legend = FALSE)
g4 <- g4 + labs(y = "Standardized deviance\n(fitted with LOESS curves)",
                x = "Prevalence (%)",
                color = "Model")
g4 <- g4 + theme_bw(base_size = 20)
# g4 <- g4 + scale_color_brewer(palette = "Pastel1")

g4 <- g4 + scale_colour_manual(breaks=levels(plot.data$ModelLabel),
                               values=c("Null model" = "#000000", 
                                        "iSDM" = "#790FBF",
                                        "UF0" = "#948B8B",
                                        "CF0" = "#048504",
                                        "CT0" = "#DB1111",
                                        "CT1" = "#030AE8",
                                        "CT2" = "#A84E05"))
print(g4)

pdf('outputs/paper 1 extensions/P1 - quality of fit [all] extended.pdf', height=9, width=14)
ggarrange(g1,g2,g4,g3, align="h", labels="auto", font.label = list(size = 16))
dev.off()

# How do correlated community parameters improve quality of fit?
plot.data.bar <- plot.data %>%
  filter(Trial %in% c("FF0", "TF0")) %>%
  as.tibble()

levels(plot.data.bar$TrialLabel)[levels(plot.data.bar$TrialLabel)=="FF0"] <- "Uncorrelated beta_comm"
levels(plot.data.bar$TrialLabel)[levels(plot.data.bar$TrialLabel)=="TF0"] <- "Correlated beta_comm"

plot.data.point <- plot.data %>%
  filter(Trial=="Null model")

g <- ggplot()
g <- g + geom_bar(data=plot.data.bar, aes(x=TaxonLabel, y=std.res.dev, fill=TrialLabel), stat="identity", position="identity", alpha=0.6)
g <- g + geom_point(data=plot.data.point, aes(x=TaxonLabel, y=std.res.dev), size=0.75, alpha=0.35)
g <- g + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g <- g + labs(title = "Standardized deviance of taxa in joint model w/ and w/o correlated beta_comm",
              subtitle = "Null model deviance given by points",
              x = "Taxa (increasing prevalence from left to right)",
              y = "Standardized deviance",
              fill = "Model")
g <- g + scale_fill_brewer(palette = "Set1")

pdf('outputs/paper 1 extensions/P1 - Quality of fit of FF0 vs FT0.pdf', width = 15, paper='a4r')
g
dev.off()
# Probabilities ####
# Plot probabilities of occurrence based on model trial and observed presence-absence: not a very effective plot
# isdm$probability$Trial <- "iSDM"
# plot.data <- bind_rows(jsdm.FF0$probability, jsdm.TF0$probability, jsdm.TT0$probability, jsdm.TT1$probability, jsdm.TT2$probability, isdm$probability)
# 
# plot.data$ObsLabel <- factor(plot.data$Obs, labels = c("0"="Absent", "1"="Present"))
# plot.data$TrialLabel <- factor(plot.data$Trial, levels=c("iSDM", "FF0", "TF0", "TT0", "TT1", "TT2"))
# 
# g <- ggplot(plot.data)
# g <- g + geom_boxplot(aes(x=TrialLabel, y=Pred, fill=ObsLabel))
# g <- g + scale_fill_brewer(palette = "Set1")



# Tjur coefficient of discrimination
isdm$tjur$Trial <- "iSDM"
tjur <- bind_rows(isdm$tjur, jsdm.FF0$tjur, jsdm.TF0$tjur, jsdm.TT0$tjur, jsdm.TT1$tjur, jsdm.TT2$tjur, jsdm.TT3$tjur, jsdm.TT4$tjur)

boxplot(data = tjur, D ~ Model, horizontal = TRUE, xlab="Tjur's R-squared")



# Predictive performance ####
cv.iSDM <- cv.isdm(predictors)
cv.FF0 <- cv.jsdm(folder = "FF0")
cv.TF0 <- cv.jsdm(folder = 'TF0')
cv.TT0 <- cv.jsdm(folder = 'TT0')
cv.TT1 <- cv.jsdm(folder = 'TT1')
cv.TT2 <- cv.jsdm(folder = 'TT2')


# Combine all deviance statistics from all models
cv.FF0$deviance$Trial <- "UF0"
cv.iSDM$deviance$Trial <- "iSDM"
cv.iSDM$deviance <- cv.iSDM$deviance[, colnames(cv.FF0$deviance)]

cv.deviance <- bind_rows(cv.iSDM$deviance, cv.FF0$deviance, cv.TF0$deviance, cv.TT0$deviance, cv.TT1$deviance, cv.TT2$deviance)

cv.deviance$Trial <- ifelse(cv.deviance$Trial=="FF0", "UF0", cv.deviance$Trial)
cv.deviance$Trial <- ifelse(cv.deviance$Trial=="FT0", "UT0", cv.deviance$Trial)
cv.deviance$Trial <- ifelse(cv.deviance$Trial=="TF0", "CF0", cv.deviance$Trial)
cv.deviance$Trial <- ifelse(cv.deviance$Trial=="TT0", "CT0", cv.deviance$Trial)

cv.deviance$Trial <- ifelse(cv.deviance$Trial=="TT1", "CT1", cv.deviance$Trial)
cv.deviance$Trial <- ifelse(cv.deviance$Trial=="TT2", "CT2", cv.deviance$Trial)

# P1 CV iSDM vs UF0 ####


plot.data <- cv.deviance %>%
  filter(Trial %in% c("iSDM", "UF0")) %>%
  select(Taxon, Type, Trial, std.deviance) %>%
  group_by(Trial, Type, Taxon) %>%
  summarise(MSD = mean(std.deviance)) %>%
  ungroup() %>%
  spread(Type, MSD)

ymax <- max(plot.data$Training)

plot.data$ymax <- ifelse(plot.data$Testing > ymax, TRUE, FALSE)

plot.data$Shape <- ifelse(plot.data$ymax, 24, 16)
plot.data$Testing <- ifelse(plot.data$ymax, max(plot.data$Training), plot.data$Testing)
plot.data$Shape <- as.factor(plot.data$Shape)
plot.data$n.present <- n.bdms[plot.data$Taxon]

g1 <- ggplot(plot.data)
g1 <- g1 + geom_point(aes(x = Training, y = Testing, colour = Trial, size = n.present, shape = Shape), alpha = 0.3)
g1 <- g1 + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
g1 <- g1 + theme(strip.background=element_rect(fill="black"), 
                 strip.text=element_text(color="white", face="bold"), 
                 axis.text=element_text(size=18))
g1 <- g1 + theme_bw(base_size = 20)
g1 <- g1 + guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))
g1 <- g1 + labs(x = "Mean standardized deviance for calibration",
                y = "Mean standardized deviance\nfor prediction",
                colour = "Model",
                shape = "Deviance",
                size = "Number of\noccurrences")
g1 <- g1 + scale_colour_manual(values=c("UF0" = "#048504", "iSDM" = "#790FBF"))
g1 <- g1 + scale_size_continuous(range = c(2, 8))
g1 <- g1 + scale_shape_discrete(name  = "Deviance",
                                breaks=c("16", "24"),
                                labels=c("In range", "Out of range"))
g1 <- g1 + coord_cartesian(ylim=c(0, ymax))
g1

# P2 CV max posterior ####
plot.data <- cv.deviance %>%
  select(Taxon, Type, Trial, Fold, std.deviance) %>%
  group_by(Trial, Type, Taxon) %>%
  summarise(MSD = mean(std.deviance)) %>%
  ungroup() %>%
  mutate(Type = factor(Type, levels=c("Training", "Testing")),
         Trial = factor(Trial, levels=c("iSDM", "UF0", "CF0", "CT0", "CT1", "CT2")))

g2 <- ggplot(plot.data)
g2 <- g2 + geom_boxplot(aes(x=Trial, y = MSD, fill=Type), show.legend=FALSE)
g2 <- g2 + theme_bw(base_size = 20)
g2 <- g2 + theme(axis.title.x = element_blank())
g2 <- g2 + ylim(c(0,1.5))
g2 <- g2 + labs(y="Mean standardized deviance\n(based on maximum posterior)",
                caption=expression(gamma==0~"and"~"u"==0~"during prediction for CT0, CT1, and CT2"))
g2 <- g2 + scale_fill_brewer(palette = "Set1")
g2 



# P3 CV full posterior ####
cv.FF0.sample <- cv.jsdm.sample(folder = "FF0")
cv.TF0.sample <- cv.jsdm.sample(folder = "TF0")
cv.TT0.sample <- cv.jsdm.sample(folder = "TT0")
cv.TT1.sample <- cv.jsdm.sample(folder = "TT1")
cv.TT2.sample <- cv.jsdm.sample(folder = "TT2")

cv.deviance.sample <- bind_rows(cv.FF0.sample$deviance, cv.TF0.sample$deviance, cv.TT0.sample$deviance, cv.TT1.sample$deviance, cv.TT2.sample$deviance)

cv.deviance.sample$Trial <- ifelse(cv.deviance.sample$Trial=="FF0", "UF0", cv.deviance.sample$Trial)
cv.deviance.sample$Trial <- ifelse(cv.deviance.sample$Trial=="FT0", "UT0", cv.deviance.sample$Trial)
cv.deviance.sample$Trial <- ifelse(cv.deviance.sample$Trial=="TF0", "CF0", cv.deviance.sample$Trial)
cv.deviance.sample$Trial <- ifelse(cv.deviance.sample$Trial=="TT0", "CT0", cv.deviance.sample$Trial)

cv.deviance.sample$Trial <- ifelse(cv.deviance.sample$Trial=="TT1", "CT1", cv.deviance.sample$Trial)
cv.deviance.sample$Trial <- ifelse(cv.deviance.sample$Trial=="TT2", "CT2", cv.deviance.sample$Trial)

cv.deviance.sample$Type[cv.deviance.sample$Type=="testing"] <- "Testing"

# Get the mean standardized deviance
plot.data <- cv.deviance.sample %>%
  select(Taxon, Type, Trial, Fold, std.deviance) %>%
  group_by(Trial, Type, Taxon) %>%
  summarise(MSD = mean(std.deviance)) %>%
  ungroup() %>%
  mutate(Type = factor(Type, levels=c("Training", "Testing")),
         Trial = factor(Trial, levels=c("UF0", "CF0", "CT0", "CT1", "CT2")))


g3 <- ggplot(plot.data)
g3 <- g3 + geom_boxplot(aes(x=Trial, y = MSD, fill=Type))
g3 <- g3 + scale_fill_brewer(palette = "Set1")
g3 <- g3 + theme_bw(base_size = 20)
g3 <- g3 + theme(axis.title.x = element_blank())
g3 <- g3 + labs(y = "Mean standardized deviance\n(based on sampled posterior)")

g3 <- g3 + ylim(c(0,1.5))
g3

# Grid arranged plots ####
pdf('outputs/paper 1 extensions/P1 - predictive performance.pdf', height=11, width=12)
bottom.row <- plot_grid(g2, g3, align="h", labels=c("b", "c"), label_size = 22)
plot_grid(g1, bottom.row, labels=c("a", ""), ncol=1, label_size = 22)
rm(bottom.row)
dev.off()

# P4 CV posterior richness ####
# richness <- bind_rows(cv.FF0.sample$richness, cv.TF0.sample$richness, cv.TT0.sample$richness, cv.TT1.sample$richness, cv.TT2.sample$richness)
# Calculate the *mean* predicted richness based on k-fold cross-validation
richness <- bind_rows(cv.FF0.sample$richness, cv.TF0.sample$richness, cv.TT0.sample$richness, cv.TT1.sample$richness, cv.TT2.sample$richness) %>%
  mutate(Trial = ifelse(Trial=="FF0", "UF0", Trial),
         # Trial = ifelse(Trial=="FT0", "UT0", Trial),
         Trial = ifelse(Trial=="TF0", "CF0", Trial),
         Trial = ifelse(Trial=="TT0", "CT0", Trial),
         Trial = ifelse(Trial=="TT1", "CT1", Trial),
         Trial = ifelse(Trial=="TT2", "CT2", Trial),
         Type = factor(Type, levels=c("Training", "Testing")),
         Trial = factor(Trial, levels=c("UF0", "CF0", "CT0", "CT1", "CT2"))) # "UT0" omitted

# Observed richness must be manually incorporated for the following plot
richness.obs <- gather(sample.bdms, Taxon, Obs, -SiteId, -SampId)

taxonomy.bdms <- lapply(unique(richness.obs$Taxon), function(j){
  taxon.rank(j, print = FALSE)
})
taxonomy.bdms <- bind_rows(taxonomy.bdms)

species <- taxonomy.bdms$Taxon[taxonomy.bdms$Rank=="Species"]
families <- taxonomy.bdms$Taxon[taxonomy.bdms$Rank=="Family" | taxonomy.bdms$Rank=="Genus" | taxonomy.bdms$Rank=="Species"]

# Calculate observed species richness (training)
richness.species <- richness.obs %>%
  as.tibble() %>%
  filter(Taxon %in% species) %>%
  group_by(SampId) %>%
  summarise(richness.obs = sum(Obs, na.rm = TRUE)) %>%
  mutate(Rank = "Species")

# Calculate observed family richness (training)
richness.family <- richness.obs %>%
  as.tibble() %>%
  filter(Obs == 1 & Taxon %in% families) %>%
  left_join(taxonomy.bdms, by="Taxon") %>%
  group_by(SampId) %>%
  summarise(richness.obs = uniqueN(Family)) %>%
  mutate(Rank = "Family")

richness.obs <- bind_rows(richness.family, richness.species)
rm(richness.species, richness.family)

# Calculate and plot richness of all hierarchical multi-species models
plot.data <- richness %>%
  filter(Type=="Testing") %>%
  ungroup() %>%
  left_join(richness.obs, by=c("SampId", "Rank"))

mean.richness <- richness %>%
  filter(Type=="Testing") %>%
  group_by(Type, Trial, Rank, SampId) %>%
  summarise(mean.richness = mean(total.pred.richness)) %>%
  ungroup() %>%
  left_join(richness.obs, by=c("SampId", "Rank"))

g2 <- ggplot()
# Plot richness over all posterior subsamples
# g2 <- g2 + geom_hex(data = plot.data, aes(x = richness.obs, y = total.pred.richness))
g2 <- g2 + geom_point(data = plot.data, aes(x = richness.obs, y = total.pred.richness), size = 0.25, color = "black", alpha = 0.1)
g2 <- g2 + geom_point(data = mean.richness, aes(x = richness.obs, y = mean.richness), size = 0.75, alpha = 0.4, color = "red")
g2 <- g2 + geom_abline(intercept = 0, slope = 1, color="black", size=1.25)
g2 <- g2 + facet_grid(Rank ~ Trial)
g2 <- g2 + theme_bw(base_size = 20)
g2 <- g2 + labs(x = "Observed sample richness",
                y = "Predicted posterior sample richness")
g2 <- g2 + coord_cartesian(ylim=c(0,75))

# pdf('outputs/paper 1 extensions/P1 - richness.pdf', height=10.5, width=13)
# print(g2)
# dev.off()

# Calculate CV richness for 
# # From results_p1.R script:
# K <- c("Temp", "Temp2", "FV", "F100m", "LUD", "IAR")
# predictors <- prepare.inputs(K, sample.bdms, center = TRUE)
# Assumes training and testing data are already in workspace (e.g., train1/test1)
# Assumes observed richness is already in workspace (richness.obs)
cv.rm <- function(predictors){
  output <- tibble()
  f <- as.formula(paste("richness.obs~", paste(K, collapse="+")))
  for (fold in 1:3){
    cat("fold:", fold, "\n")
    train.data <- as.tibble(get(paste("train",fold,sep="")))
    test.data <- as.tibble(get(paste("test",fold,sep="")))
    
    # Observation data only for training/testing of family/species
    train.richness.obs.family <- train.data %>%
      select(SampId) %>%
      left_join(filter(richness.obs, Rank=="Family"), by="SampId")
    
    train.richness.obs.species <- train.data %>%
      select(SampId) %>%
      left_join(filter(richness.obs, Rank=="Species"), by="SampId")
    
    test.richness.obs.family <- test.data %>%
      select(SampId) %>%
      left_join(filter(richness.obs, Rank=="Family"), by="SampId")
    
    test.richness.obs.species <- test.data %>%
      select(SampId) %>%
      left_join(filter(richness.obs, Rank=="Species"), by="SampId")
    
    # Complete observations and inputs for GLM calibration/prediction
    train.family.model.data <-  train.richness.obs.family %>%
      left_join(predictors, by="SampId")
    
    train.species.model.data <-  train.richness.obs.species %>%
      left_join(predictors, by="SampId")
    
    test.family.model.data <- test.richness.obs.family %>%
      left_join(predictors, by="SampId")
    
    test.species.model.data <-  test.richness.obs.species %>%
      left_join(predictors, by="SampId")
    
    # Calibrate the GLMs using the training data
    rmf <- glm.nb(f, data = train.family.model.data, link=log)
    rms <- glm.nb(f, data = train.species.model.data, link=log)
    
    test.pred.family <- tibble(SampId = test.richness.obs.family$SampId,
                               richness.obs = test.richness.obs.family$richness.obs,
                               richness.pred = predict(rmf, newdata = test.family.model.data, type= "response"),
                               Rank = "Family",
                               Type = "Testing",
                               Fold = fold,
                               Trial = "Richness model")
    
    test.pred.species <- tibble(SampId = test.richness.obs.species$SampId,
                                richness.obs = test.richness.obs.species$richness.obs,
                                richness.pred = predict(rms, newdata = test.species.model.data, type= "response"),
                                Rank = "Species",
                                Type = "Testing",
                                Fold = fold,
                                Trial = "Richness model")
    
    train.pred.family <- tibble(SampId = train.family.model.data$SampId, 
                                richness.obs = train.family.model.data$richness.obs,
                                richness.pred = rmf$fitted.values,
                                Rank = "Family",
                                Type = "Training",
                                Fold = fold,
                                Trial = "Richness model")
    
    train.pred.species <- tibble(SampId = train.species.model.data$SampId, 
                                 richness.obs = train.species.model.data$richness.obs,
                                 richness.pred = rms$fitted.values,
                                 Rank = "Species",
                                 Type = "Training",
                                 Fold = fold,
                                 Trial = "Richness model")
    output <- bind_rows(output, train.pred.family, train.pred.species, test.pred.family, test.pred.species)
    
  }
  return(output)
}

data <- spread(richness.obs, Rank, richness.obs)
data <- left_join(data, predictors, by="SampId")
rmf <- glm.nb(Family ~ Temp + Temp2 + FV + F100m + LUD + IAR, data = data, link=log)
rms <- glm.nb(Species ~ Temp + Temp2 + FV + F100m + LUD + IAR, data = data, link=log)

# Calculate predicted richness of richness model
richness.rm <- cv.rm(predictors)

richness.rm <- richness.rm %>%
  group_by(SampId, Rank, Type) %>%
  summarise(richness.pred = mean(richness.pred)) %>%
  ungroup() %>%
  mutate(Trial = "Richness model") %>%
  left_join(richness.obs, by=c("SampId", "Rank"))

# Calculate predicted richness of multi-species model
cv.FF0 <- cv.jsdm(folder = "FF0")

richness.family.FF0 <- cv.FF0$probability %>%
  left_join(taxonomy.bdms, by = "Taxon") %>%
  # Filter the obs/taxonomy and join the joint model probabilities 
  filter(!is.na(Obs), Rank == 'Family' | Rank == 'Genus' | Rank == 'Species') %>%
  group_by(SampId, Taxon, Type) %>%
  summarise(richness.pred = mean(Pred)) %>%
  ungroup() %>%
  left_join(taxonomy.bdms, by = "Taxon") %>%
  # Calculate the predicted richness by family at each sample
  group_by(SampId, Family, Type) %>%
  summarise(richness.pred = 1-prod(1-richness.pred)) %>%
  group_by(SampId, Type) %>%
  summarise(richness.pred = sum(richness.pred)) %>% 
  # Remove grouping, select columns, and add additional group information
  ungroup() %>%
  select(SampId, richness.pred, Type) %>%
  mutate(Rank = "Family", Trial = "UF0")

richness.species.FF0 <- cv.FF0$probability %>%
  left_join(taxonomy.bdms, by = "Taxon") %>%
  filter(!is.na(Obs), Rank == "Species") %>%
  group_by(SampId, Taxon, Type) %>%
  summarise(richness.pred = mean(Pred)) %>%
  group_by(SampId, Type) %>%
  summarise(richness.pred = sum(richness.pred)) %>%
  ungroup() %>%
  select(SampId, richness.pred, Type) %>%
  mutate(Rank = "Species", Trial = "UF0")

richness.FF0 <- bind_rows(richness.family.FF0, richness.species.FF0) %>% left_join(richness.obs, by = c("SampId", "Rank"))
rm(richness.family.FF0, richness.species.FF0)

table(colnames(richness.FF0) %in% colnames(richness.rm))

richness.FF0 <- richness.FF0[, colnames(richness.rm)]

plot.data <- bind_rows(richness.FF0, richness.rm)
plot.data <- filter(plot.data, Type == "Testing")

g1 <- ggplot(plot.data)
g1 <- g1 + geom_point(aes(x = richness.obs, y = richness.pred, color = Trial), alpha = 0.5)
# g1 <- g1 + facet_grid(Type ~ Rank)
g1 <- g1 + facet_grid(. ~ Rank)
g1 <- g1 + geom_abline(intercept = 0, slope = 1, color="black", size=1.25)
g1 <- g1 + theme_bw(base_size = 20)
g1 <- g1 + labs(y = "Predicted sample richness",
                x = "Observed sample richness",
                color = "Model")
g1 <- g1 + guides(colour = guide_legend(override.aes = list(size=6)))
g1 <- g1 + scale_color_brewer(palette = "Set1")
g1



# Grid arranged plots ####
pdf('outputs/paper 1 extensions/P1 - richness.pdf', height=10.5, width=13)
plot_grid(g1, g2, labels=c("a", "b"), ncol=1, label_size = 20)
dev.off()


tiff("outputs/paper 1 extensions/P1 - richness.tiff", height = 10.5, width = 13, units = 'in', compression = "lzw", res = 400)
plot_grid(g1, g2, labels=c("a", "b"), ncol=1, label_size = 20)
dev.off()

# Example taxa ####
# Plot map.jsdm.pred.taxon() examples ####
# Grid arrange prob.taxon and map.taxon plots
g1 <- plot.prob.taxon(fit.FF0, "Gammaridae", legend=FALSE)
g2 <- map.jsdm.pred.taxon(fit.FF0, "Gammaridae", legend=FALSE)
g3 <- map.jsdm.pred.taxon(fit.FF0, "Nemoura_minima", legend=FALSE)
g4 <- map.jsdm.pred.taxon(fit.FF0, "Protonemura_lateralis", legend=FALSE)

pdf('outputs/paper 1 extensions/P1 - maps example taxa.pdf', width=13, height=9.5)
plot_grid(g1,g2,g3,g4, labels=c("a", "b", "c", "d"), align="hv")
dev.off()

# plot.comm() [all taxa] ####
plot.comm(jsdm.FF0, 'P1 - plot_comm [all taxa]')

# plot.comm() [examples] ####
# Modified plot.comm() code for specific taxa.
response.bdms <- extract.resp(fit.FF0)
beta.samples <- extract.beta(fit.FF0)
taxa <- c("Gammaridae", "Nemoura_minima", "Protonemura_lateralis")

pdf('outputs/paper 1 extensions/P1 - plot_comm [example taxa].pdf', onefile = TRUE)
par(cex=1.25)
for (k in 1:length(inf.fact)){
  variable <- inf.fact[k]
  responses <- response.bdms[response.bdms$Taxon %in% taxa, ]
  
  response <- responses[[variable]]
  names(response) <- responses$Taxon
  
  response[response==0] <- "grey55"
  response[response==1] <- "blue"
  response[response==-1] <- "red"
  cat("Plotting: ",variable,"\n")
  # Test temperature as starting variable
  samples <- beta.samples[Variable == variable,]
  
  # Find the maximum density among the posterior taxon-specific posterior parameters
  samples.sd <- samples %>%
    group_by(Taxon) %>%
    summarise(SD = sd(Value)) %>%
    arrange(SD) # arrange from min to max
  
  sd <- samples.sd[1,]$Taxon
  ymax.sample <- samples[Taxon == sd, ]
  ymax <- density(ymax.sample$Value)
  ymax <- max(ymax$y)
  
  # Plot the community parameter distribution
  mu <- jsdm$mu.beta.comm.maxpost[variable]
  sd <- jsdm$sigma.beta.comm.maxpost[variable]
  x <- seq(mu-4*sd,mu+4*sd,length.out=201)
  beta.comm.density <- dnorm(x, mu, sd)
  beta.taxa.maxpost.density <- density(jsdm$beta.taxa.maxpost[variable, ])
  
  # Match expressions to influence factors
  labels <- c("A10m" = expression(paste(beta["A10m"], " (1/%)")),
              "Temp" = expression(paste(beta["Temp"], " (1/", degree, "C)")),
              "Temp2" = expression(paste(beta["Temp"^2], " (1/", degree, "C"^2,")")),
              "FV" = expression(paste(beta["FV"], " (s/m)")),
              "bFRI" = expression(paste(beta["bFRI"], " (1/%)")),
              "FRI" = expression(paste(beta["FRI"], " (1/%)")),
              "F10m" = expression(paste(beta["F10m"], " (1/%)")),
              "F100m" = expression(paste(beta["F100m"], " (1/%)")),
              "IAR" = expression(paste(beta["IAR"], " 1/(w"["c"]%*%"f"["c"],")")),
              "Urban" = expression(paste(beta["Urban"], " (1/%)")),
              "LUD" = expression(paste(beta["LUD"], " (km"^2,"/CE)"))
  )
  
  plot(numeric(0), numeric(0),
       xlim = c(min(x), max(x)),
       ylim=c(0, ymax-(ymax*0.25)),
       xlab = labels[variable],
       ylab = "Density")
  abline(v=0)
  # Plot the taxon-specific parameters
  
  l <- c("G", "N", "P")
  
  for (j in 1:length(taxa)){
    taxon <- taxa[j]
    sample <- samples[Taxon==taxon,]
    sample <- sample$Value
    
    lines(density(sample), type="l", col=alpha(response[taxon], 1), lwd=1.25)
    text(x=mean(sample), y=max(density(sample)$y)+0.1*max(density(sample)$y), labels = l[j])
  }
  # Plot community parameter distribution
  lines(x, beta.comm.density, type="l", col = "grey50", lwd = 5)
  # Plot maximum posterior values over all taxa
  lines(beta.taxa.maxpost.density, type="l", col="black", lwd=2, lty='longdash')
  
}
dev.off()

# Plot P1 - parameter dotplot [4 pages] ####
# Warning: really bad code ahead...
# Extract the entire posterior beta sample
jsdm.beta.samples <- extract.beta(jsdm.FF0)

jsdm.beta.samples.quantiles <- jsdm.beta.samples %>%
  group_by(Variable, Taxon) %>%
  summarise(quantile10=quantile(Value, 0.10), quantile90=quantile(Value, 0.90)) %>%
  setDT()
colnames(jsdm.beta.samples.quantiles)<- c("Variable", "Taxon", "quantile10", "quantile90")

isdm.beta <- isdm$parameters
if ("Model" %in% colnames(isdm.beta)){
  isdm.beta$Model <- NULL
}
colnames(isdm.beta) <- c("Taxon", "Variable", "iSDM.parameter")


# Number of pages per taxa (ntpp)
pages <- 4
ntpp <- length(n.bdms)/pages
tpp <- split(n.bdms, ceiling(seq_along(n.bdms)/ntpp)) # bin the number of taxa per page
pdf(paste("outputs/paper 1 extensions/parameter_dotplot_notnorm.pdf",sep=""), paper = "special", height = 10, width = 9, onefile = TRUE)
for (page in 1:pages){
  
  plot.data <- jsdm.beta.samples.quantiles[jsdm.beta.samples.quantiles$Taxon %in% names(tpp[page][[1]]), ]
  
  # Prepare the jSDM parameter values
  beta.max <- select(jsdm.FF0$parameters, Taxon, Variable, Parameter)
  colnames(beta.max) <- c("Taxon", "Variable", "jSDM.parameter")
  beta.max$Variable <- as.character(beta.max$Variable)
  
  # Join the quantiles, jSDM, and iSDM parameters together
  plot.data <- left_join(plot.data, beta.max, by = c("Taxon", "Variable"))
  plot.data <- left_join(plot.data, isdm.beta, by = c("Taxon", "Variable"))
  
  # Model (variable) and Parameter (value) gather the two columns, jSDM.parameter and iSDM.parameter
  # into a variable/value pair
  plot.data <- gather(plot.data, Model, Parameter, jSDM.parameter:iSDM.parameter)
  
  # Process iSDM outliers
  # Get min and max posterior parameter values within the 10th-90th quantiles
  lim.min <- aggregate(quantile10 ~ Variable, plot.data, min)
  vmin <- lim.min$quantile10; names(vmin) <- lim.min$Variable
  
  lim.max <- aggregate(quantile90 ~ Variable, plot.data, max)
  vmax <- lim.max$quantile90; names(vmax) <- lim.max$Variable
  
  # Set colors for each parameter
  # g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
  plot.data$Col <- "#048504"
  plot.data$Col[plot.data$Model == 'iSDM.parameter'] <- "#790FBF"
  plot.data$Col <- ifelse(plot.data$Parameter < vmin[plot.data$Variable], "#6B6B6B", plot.data$Col)
  plot.data$Parameter <- ifelse(plot.data$Parameter < vmin[plot.data$Variable], vmin[plot.data$Variable], plot.data$Parameter)
  
  plot.data$Col <- ifelse(plot.data$Parameter > vmax[plot.data$Variable], "#6B6B6B", plot.data$Col)
  plot.data$Parameter <- ifelse(plot.data$Parameter > vmax[plot.data$Variable], vmin[plot.data$Variable], plot.data$Parameter)
  rm(lim.min, lim.max, vmin, vmax)
  
  # Order labels by occurrence frequency
  plot.data$Labels <- factor(paste(plot.data$Taxon, ' - ', n.bdms[plot.data$Taxon]), levels = paste(names(sort(n.bdms)), ' - ', sort(n.bdms)))
  
  # Order variable facets and pass expressions for units in facet labels
  plot.data$Variable <- factor(plot.data$Variable, levels = c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD"))
  levels(plot.data$Variable) <- labeller(c("Temp", "Temp2", "FV", "F100m", "IAR", "Urban", "LUD"))
  
  # Build the plot
  g <- ggplot(plot.data)
  g <- g + geom_hline(yintercept=0, alpha=0.4)
  g <- g + geom_linerange(aes(x = Labels, ymin = quantile10, ymax = quantile90), color = 'black', alpha = 0.6)
  g <- g + geom_point(aes(x = Labels, y = Parameter, colour = Col), stroke=0, size = 2.5, alpha = 0.5)
  g <- g + facet_grid(. ~ Variable, scales = "free", labeller=label_parsed)
  g <- g + coord_flip()
  g <- g + theme_bw()
  g <- g + labs(title = "Taxon-specific parameter estimates by model",
                subtitle = "MLEs (iSDM) and maximum posterior (UF0) with 10th-90th percentile interval",
                x = "Taxon and prevalance",
                y = "Value")
  g <- g + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  g <- g + scale_colour_identity()
  print(g)
}
dev.off()

# Plot P1 - parameter SD ####
# Extract the entire posterior beta sample
jsdm.beta.samples <- extract.beta(jsdm.FF0)

# Fast aggregation of 10M row dataset
plot.data <- jsdm.beta.samples %>%
  group_by(Taxon, Variable) %>%
  summarise(SD = sd(Value), Mean = mean(Value)) %>%
  mutate(n = n.bdms[Taxon], Label = factor(Variable, levels=c("Temp", "Temp2", "FV", "F10m",  "IAR",   "Urban", "LUD")), rSD = SD/abs(Mean))

levels(plot.data$Label) <- labeller(levels(plot.data$Label))
setDT(plot.data)
g <- ggplot(data=plot.data, aes(x = n, y = SD, size = n))
g <- g + geom_point(alpha = 0.5)
g <- g + facet_grid(Label ~ ., scales = "free", labeller=label_parsed)
g <- g + theme_bw(base_size = 14)
g <- g + theme(strip.background=element_rect(fill="grey"),strip.text=element_text(color="black", face="bold"),
               plot.title = element_text(hjust = 0.5, size = 12))
g <- g + labs(title = expression(paste("Standard deviation of posterior taxon-specific parameter distributions ", beta["jk"]^"taxa", " in UF0")),
              x = "Occurrence frequency",
              y = expression(paste("Standard deviation (", sigma[beta["jk"]^"taxa"],")")),
              size = "Occurrence frequency")
g <- g + scale_y_continuous(limits=c(0,NA))

pdf('outputs/paper 1 extensions/P1 - parameter SD.pdf', height = 12.5, width = 9)
print(g)
dev.off()


plot.data$Variable <- factor(plot.data$Variable, levels=inf.fact)

g <- ggplot(data=plot.data)
g <- g + geom_boxplot(aes(x=Variable, y=rSD, fill=Variable))
g <- g + coord_cartesian(ylim=c(0,10))
g <- g + theme_bw()
g <- g + theme(plot.title = element_text(hjust = 0.5, size = 12), 
               axis.text = element_text(size = 12),
               axis.text.x = element_text(angle=45, hjust=1, vjust=0.5))
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + guides(fill=FALSE)
g <- g + labs(title = expression(paste("Relative uncertainty of posterior ", beta["jk"]^"taxa", " in UF0")),
              x = "Explanatory variable",
              y = expression(paste("Relative standard deviation (", sigma[beta["jk"]^"taxa"]," / |",mu[beta["jk"]^"taxa"],"|)")))

pdf('outputs/paper 1 extensions/P1 - parameter uncertainty.pdf')
print(g)
dev.off()

# Map inputs ####
x <- prepare.inputs(K, sample.bdms, center = FALSE)
x <- gather(x, Variable, Value, -SiteId, -SampId)
x <- filter(x, Variable != "Temp2")
x$Trial <- 'BDM species'
map.inputs(x, 'P1 - maps influence factors')

# Plot P1 - map gamma ####
dt <- tibble(SiteId = names(jsdm.TT0$gamma.maxpost), gamma.maxpost = jsdm.TT0$gamma.maxpost)
dt <- left_join(dt, inputs$xy, by="SiteId")
dt$color <- ifelse(dt$gamma.maxpost > 0, "Positive", "Negative")

g <- ggplot()
g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)

g <- g + geom_point(data = dt, aes(X, Y, size = abs(gamma.maxpost), color = color), alpha = 0.35)
g <- g + scale_size_continuous(breaks=seq(0, 3.5, 0.5), limits=c(0,3.5), range = c(1, 8))

g <- g + scale_color_manual(values= c("Positive" = "blue", "Negative" = "red"))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title = "Site effect (CT0)",
              size = expression(paste(gamma["i"], " (absolute value)")),
              color = "")
g <- g + theme_void(base_size = 15)
g <- g + theme(plot.title = element_text(hjust = 0.5),
               panel.grid.major = element_line(colour="transparent"))
print(g)

pdf('P1 - map site effect.pdf', paper = 'special', width = 10.5, onefile = TRUE)
print(g)
dev.off()

x <- prepare.inputs(K, sample.bdms, center = FALSE)
plot.data <- left_join(x, dt, by="SiteId")
pairs.panels(select(plot.data, -SiteId, -SampId), density = TRUE, scale=FALSE, hist.col="grey", cex.cor=1.5, cex.labels=1.5)

# Plot P1 - map richness ####
# Gather observations into narrow dataset
obs <- data.table(gather(sample.bdms, Taxon, Obs, -SiteId, -SampId))

rank <- lapply(unique(obs$Taxon), function(j){
  taxon.rank(j) # ID the taxonomic resolution (i.e., "rank") of each taxon
})
rank <- rbindlist(rank)

# Obtain all observations along with full taxonomy and taxonomic level of identification
obs <- merge(obs, rank, by = "Taxon", all.x = TRUE)

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

g <- ggplot()
g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)

g <- g + geom_point(data = obs.site.richness, aes(X, Y, size = ObsFamilies), alpha = 0.35)

g <- g + scale_radius(limits = c(0,40), breaks = seq(10, 40, 10))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(size = "Observed family\nrichness")
g <- g + theme_void(base_size = 15)
g <- g + theme(plot.title = element_blank(),
               panel.grid.major = element_line(colour="transparent"))

pdf('P1 - map family richness.pdf', paper = 'special', width = 11, onefile = TRUE)
print(g)
dev.off()

g <- ggplot()
g <- g + geom_sf(data = inputs$ch, fill=NA, color="black")
g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)

g <- g + geom_point(data = obs.site.richness, aes(X, Y, size = ObsSpecies), alpha = 0.35)
g <- g + labs(size = "Observed species\nrichness")

g <- g + scale_radius(limits = c(0,30), breaks =seq(5, 30, 5))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + theme_void(base_size = 15)
g <- g + theme(plot.title = element_blank(),
               panel.grid.major = element_line(colour="transparent"))

pdf('P1 - map species richness.pdf', paper = 'special', width = 11, onefile = TRUE)
print(g)
dev.off()

# Get the mean standardized deviance
plot.data <- cv.deviance %>%
  select(Taxon, Type, Trial, std.deviance) %>%
  group_by(Trial, Type, Taxon) %>%
  summarise(MSD = mean(std.deviance)) %>%
  ungroup() %>%
  # spread(Type, MSD) %>%
  mutate(n = n.bdms[Taxon])

# plot.data[is.infinite(plot.data)] <- NA
g <- ggplot(plot.data)
# g <- g + geom_point(aes(x=Training, y=Testing, size=n, color=Trial))
g <- g + geom_boxplot(aes(x=Trial, y = MSD, fill=Type))
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + labs(title="predictive performance under k-fold cross-validation",
              y="Mean standardized deviance")
g <- g + ylim(c(0,2.5))
g

plot.data <- cv.deviance %>%
  select(Taxon, Type, Trial, std.deviance) %>%
  group_by(Trial, Type, Taxon) %>%
  summarise(MSD = mean(std.deviance)) %>%
  ungroup() %>%
  spread(Type, MSD) %>%
  filter(Trial %in% c("FF0", "iSDM"))

ymax <- max(plot.data$Training)

# cv.plot$Infinite <- is.infinite(cv.plot$Testing)
plot.data$ymax <- ifelse(plot.data$Testing > ymax, TRUE, FALSE)

plot.data$Shape <- ifelse(plot.data$ymax, 24, 16)
plot.data$Testing <- ifelse(plot.data$ymax, max(plot.data$Training), plot.data$Testing)
plot.data$Shape <- as.factor(plot.data$Shape)
plot.data$n.present <- n.bdms[plot.data$Taxon]

g <- ggplot(plot.data, aes(x = Training, y = Testing, colour = Trial, size = n.present, shape = Shape))
g <- g + geom_point(alpha = 0.3)
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
g <- g + theme(strip.background=element_rect(fill="black"), strip.text=element_text(color="white", face="bold"), axis.text=element_text(size=18))
g <- g + theme_bw(base_size = 18)
g <- g + guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))
g <- g + labs(x = "Mean standardized deviance for calibration",
              y = "Mean standardized deviance for prediction",
              colour = "Model",
              shape = "Deviance",
              size = "Prevalence")
g <- g + scale_colour_manual(values=c("FF0" = "#048504", "iSDM" = "#790FBF"))
g <- g + scale_size_continuous(range = c(2, 8))
g <- g + scale_shape_discrete(name  = "Deviance",
                              breaks=c("16", "24"),
                              labels=c("In range", "Out of range"))
g <- g + coord_cartesian(ylim=c(0, ymax))
print(g)

library(ggalt)
plot.data$TaxonLabel <- factor(plot.data$Taxon, levels=c(names(sort(n.bdms, decreasing=FALSE))))
plot.data$TrialLabel <- as.factor(plot.data$Trial)
levels(plot.data$TrialLabel)[levels(plot.data$TrialLabel)=="FF0"] <- "Uncorrelated"
levels(plot.data$TrialLabel)[levels(plot.data$TrialLabel)=="TF0"] <- "Correlated"

g <- ggplot(plot.data)
g <- g + geom_dumbbell(aes(x=Training, xend=Testing, y=TaxonLabel, group=TrialLabel, color=TrialLabel), size=1.25)
g <- g + labs(title = "Predictive performance w/ and w/o correlated community parameters",
              x = "Mean standardized deviance",
              y = "Taxa (ordered by occurrence frequency)")
g <- g + theme(axis.text.x = element_blank())
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + scale_color_brewer(palette = "Set1")
g <- g + coord_flip()
g


plot.data <- cv.jsdm.deviance %>%
  select(Taxon, Type, Trial, std.deviance) %>%
  group_by(Trial, Type, Taxon) %>%
  summarise(MSD = mean(std.deviance)) %>%
  ungroup()

g <- ggplot(plot.data)
g <- g + geom_boxplot(aes(x = Type, y = MSD, fill=Trial))
g <- g + scale_color_brewer(palette = "Set1")
g <- g + labs(y = "Mean standardized deviance",
              fill = "Model type")
g

# Map x_lat ####
# Problem: for multiple latent variables, multiple maps are generated but they size their points differently. Solution would combine geom_sf with facet_grid, which requires editing the input object
plot.data <- as.tibble(melt(fit.TT2$x.lat.maxpost))
colnames(plot.data) <- c("SampId", "x.lat", "Value")
plot.data$SiteId <- fit.TT2$sites[plot.data$SampId]
plot.data <- left_join(plot.data, inputs$xy, by="SiteId")

plot.data$abs.value <- abs(plot.data$Value)
plot.data$Color <- ifelse(plot.data$Value > 0, 1, 0)
plot.data$Color <- as.factor(plot.data$Color)

plot.data1 <- filter(plot.data, x.lat==1)
g1 <- ggplot()
g1 <- g1 + geom_sf(data = inputs$ch, fill=NA, size=1.25, color="black", show.legend = FALSE)
g1 <- g1 + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
g1 <- g1 + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)

g1 <- g1 + geom_point(data = plot.data1, aes(x=X, y=Y, size = abs.value, color = Color), alpha=0.35)
g1 <- g1 + labs(title = expression(paste(xi[paste("i","t"["i"],"l")],", l = 1")),
                color = "",
                size = "Absolute\nvalue")
g1 <- g1 + theme_void(base_size = 15)
g1 <- g1 + theme(panel.grid.major = element_line(colour="transparent"))
g1 <- g1 + guides(colour = guide_legend(override.aes = list(size=6)))
g1 <- g1 + scale_color_manual(values=c("0" = "#FF0000", "1" = "#0077FF"), labels = c("Negative", "Positive"))
g1 <- g1 + scale_radius(range = c(2, 8))

# plot(c(0,0), main = expression(paste(xi[paste("i","t"["i"],"l")],", l=2")))
plot.data2 <- filter(plot.data, x.lat==2)
g2 <- ggplot()
g2 <- g2 + geom_sf(data = inputs$ch, fill=NA, size=1.25, color="black", show.legend = FALSE)
g2 <- g2 + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
g2 <- g2 + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)

g2 <- g2 + geom_point(data = plot.data2, aes(x=X, y=Y, size = abs.value, color = Color), alpha=0.35)
g2 <- g2 + labs(title = expression(paste(xi[paste("i","t"["i"],"l")],", l = 2")),
                color = "",
                size = "Absolute\nvalue")
g2 <- g2 + theme_void(base_size = 15)
g2 <- g2 + theme(panel.grid.major = element_line(colour="transparent"))
g2 <- g2 + guides(colour = guide_legend(override.aes = list(size=6)))
g2 <- g2 + scale_color_manual(values=c("0" = "#FF0000", "1" = "#0077FF"), labels = c("Negative", "Positive"))
g2 <- g2 + scale_radius(range = c(2, 8))

pdf('outputs/paper 1 extensions/P1 - maps latent variable.pdf', width = 10, paper='a4r', onefile = TRUE)
print(g1)
print(g2)
dev.off()

# Plot beta_lat ####
# Plot the latent variable and its significance per taxon
beta.lat <- jsdm.TT1$TT1$beta.lat
beta.lat <- as.tibble(beta.lat)
row.names(beta.lat) <- 1:nrow(beta.lat)
colnames(beta.lat) <- colnames(jsdm.TT1$TT1$occur.taxa)

# Identify the taxon with the lowest SD for beta_lat
beta.lat.sd <- apply(beta.lat, 2, sd)
beta.lat.sd <- sort(beta.lat.sd, decreasing = FALSE)
beta.lat.sd.taxon <- names(which.min(beta.lat.sd))

beta.lat <- melt(beta.lat)
colnames(beta.lat) <- c("Taxon", "Value")
beta.lat <- as.tibble(beta.lat)

beta.lat.ymax <- filter(beta.lat, Taxon==beta.lat.sd.taxon)
ymax <- max(density(beta.lat.ymax$Value)$y)

pdf('outputs/paper 1 extensions/P1 - latent variables.pdf', onefile = TRUE)
plot(numeric(0), numeric(0),
     xlim = c(min(beta.lat$Value), max(beta.lat$Value)),
     ylim=c(0, ymax),
     xlab = "Marginal posterior beta_lat",
     ylab = "Density")
abline(v=0)
rm(beta.lat.sd, beta.lat.sd.taxon, beta.lat.ymax, ymax)
taxa <- occur.freq(jsdm.TT1$TT1$occur.taxa)

for (j in 1:length(taxa)){
  taxon <- names(taxa)[j]
  sample <- filter(beta.lat, Taxon==taxon)
  
  # Fill matrix: beta.taxa.response
  significance <- quantile(sample$Value, probs = c(0.05, 0.95))
  
  # If 5th quantile greater than 0, set positive
  # If 95th quantile less than 0, set negative
  if (significance[1] > 0){ # significant positive
    sig <- "blue"
  }
  if (significance[2] < 0){ # significant negative
    sig <- "red"
  }
  # If posterior is !(positive) AND !(negative), set grey
  if (!(significance[1] > 0) & !(significance[2] < 0)){
    sig <- "grey55"
  }
  
  lines(density(sample$Value), type="l", col=alpha(sig, 0.30), lwd=1)
  # cat("Taxon: ", taxon, " | significance: ", sig, "\n")
}
dev.off()

# Map number of significant interactions
extract.lat.resp <- function(jsdm){
  n <- occur.freq(jsdm$occur.taxa)
  sites <- jsdm$env.cond$SiteId
  
  x.lat <- jsdm$x.lat
  beta.lat <- jsdm$beta.lat
  
  significance.beta.lat <- apply(beta.lat, 2, function(j){
    quantiles <- quantile(j, probs = c(0.05, 0.95))
    # If 5th quantile greater than 0, set positive
    # If 95th quantile less than 0, set negative
    if (quantiles[1] > 0){ # significant positive
      significance <- 1
    }
    if (quantiles[2] < 0){ # significant negative
      significance <- -1
    }
    # If posterior is !(positive) AND !(negative), set grey
    if (!(quantiles[1] > 0) & !(quantiles[2] < 0)){
      significance <- 0
    }
    return(significance)
  })
  names(significance.beta.lat) <- names(n)
  
  significance.x.lat <- apply(x.lat, 2, function(i){
    quantiles <- quantile(i, probs = c(0.05, 0.95))
    # If 5th quantile greater than 0, set positive
    # If 95th quantile less than 0, set negative
    if (quantiles[1] > 0){ # significant positive
      significance <- 1
    }
    if (quantiles[2] < 0){ # significant negative
      significance <- -1
    }
    # If posterior is !(positive) AND !(negative), set grey
    if (!(quantiles[1] > 0) & !(quantiles[2] < 0)){
      significance <- 0
    }
    return(significance)
  })
  names(significance.x.lat) <- names(jsdm$x.lat.maxpost)
  
  output <- list("significance.x.lat" = significance.x.lat,
                 "significance.beta.lat" = significance.beta.lat)
  return(output)
}

jsdm.TT1$significance <- extract.lat.resp(jsdm.TT1)

# plot.data <- tibble(SiteId = names(jsdm.TT1$significance$significance.x.lat),
#                         x.lat = jsdm.TT1$TT1$x.lat.maxpost,
#                         sig.x.lat = jsdm.TT1$significance$significance.x.lat,
#                         sig.beta.lat.positive = length(jsdm.TT1$significance$significance.beta.lat[jsdm.TT1$significance$significance.beta.lat==1]),
#                         sig.beta.lat.negative = length(jsdm.TT1$significance$significance.beta.lat[jsdm.TT1$significance$significance.beta.lat==-1]),
#                         sig.beta.lat.neutral = length(jsdm.TT1$significance$significance.beta.lat[jsdm.TT1$significance$significance.beta.lat==0])
#                         )

# Correlation matrix ####
library(corrplot)
# Attempt 1: 
# Compute the linear predictor for the latent variable and try to plot it
n <- occur.freq(jsdm.TT1$occur.taxa)
beta.lat.maxpost <- jsdm.TT1$beta.lat.maxpost[names(n)]

lv <- sapply(beta.lat.maxpost, function(j){
  j * jsdm.TT1$x.lat.maxpost
})

lv.corr <- lv[, names(n[n>10])]
lv.corr <- cor(lv.corr)

colnames(lv.corr) <- 1:ncol(lv.corr)
rownames(lv.corr) <- 1:nrow(lv.corr)
corrplot(lv.corr, method="circle", type="upper", order="hclust", tl.pos="n")

# Attempt 2: rereading Warton (2015) shows they plot correlations among the factor loadings (lambda_j), i.e., the beta^lat
# Just plot the correlations among the marginal posterior distributions of the factor loadings!
beta.lat <- jsdm.TT1$TT1$beta.lat
beta.lat.maxpost <- jsdm.TT1$TT1$beta.lat.maxpost

# Identify taxa with significant factor loadings and n > 10
beta.lat.sig <- names(jsdm.TT1$significance$significance.beta.lat[jsdm.TT1$significance$significance.beta.lat != 0])
n <- occur.freq(jsdm.TT1$TT1$occur.taxa)
taxa <- intersect(names(n[n>10]), beta.lat.sig)

beta.lat.taxa <- beta.lat[, taxa]

# Keep taxa only with significant factor loadings AND more than 10 presence observations
lv.corr <- cor(beta.lat.taxa)

colnames(lv.corr) <- 1:ncol(lv.corr)
rownames(lv.corr) <- 1:nrow(lv.corr)

corrplot(lv.corr, method="square", diag = FALSE, type="lower", tl.pos="n")

# # take the max posterior
# beta.lat.taxa <- beta.lat.maxpost[taxa]
# lv.matrix <- matrix(beta.lat.taxa, nrow=length(taxa), ncol=length(taxa), byrow = TRUE)
# lv.corr.maxpost <- cor(lv.matrix)
# corrplot(lv.corr.maxpost)

# Parameter uncertainty ####
beta.samples.FF0 <- extract.beta(jsdm.FF0)
beta.samples.TF0 <- extract.beta(jsdm.TF0)

plot.data <- bind_rows(beta.samples.FF0, beta.samples.TF0) %>%
  group_by(Taxon, Variable, Trial) %>%
  summarise(SD = sd(Value), Mean = mean(Value)) %>%
  mutate(n = n.bdms[Taxon], 
         VariableLabel = factor(Variable, levels=c("Temp", "Temp2", "FV", "F100m",  "IAR", "LUD")),
         relativeSD = SD/abs(Mean)) %>%
  
  
  
  
  # levels(plot.data$VariableLabel) <- labeller(levels(plot.data$VariableLabel))
g <- ggplot(data=plot.data)
g <- g + geom_boxplot(aes(x=Variable, y=relativeSD, fill=Variable))
g <- g + coord_cartesian(ylim=c(0,10))
g <- g + facet_grid(.~Trial)
g <- g + theme_bw()
# g <- g + theme(plot.title = element_text(hjust = 0.5, size = 12), 
#                axis.text = element_text(size = 12),
#                axis.text.x = element_text(angle=45, hjust=1, vjust=0.5))
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + guides(fill=FALSE)
g <- g + labs(title = expression(paste("Relative uncertainty of posterior taxon-specific parameter distributions ", beta["jk"]^"taxa")),
              x = "Explanatory variable",
              y = expression(paste("Relative standard deviation (", sigma[beta["jk"]^"taxa"]," / |",mu[beta["jk"]^"taxa"],"|)")))

# Linear predictor ####
z.FF0 <- linear.predictor(jsdm.FF0)
z.FF0$Trial <- "UF0"

z.TF0 <- linear.predictor(jsdm.TF0)
z.TF0$Trial <- "CF0"

z.TT0 <- linear.predictor(jsdm.TT0)
z.TT0$Trial <- "CT0"

z.TT1 <- linear.predictor(jsdm.TT1)
z.TT1$Trial <- "CT1"

z.TT2 <- linear.predictor(jsdm.TT2)
z.TT2$Trial <- "CT2"


# z.TT3 <- linear.predictor(jsdm.TT3)
# z.TT3$Trial <- "TT3"
# 
# z.TT4 <- linear.predictor(jsdm.TT4)
# z.TT4$Trial <- "TT4"

z.jsdm <- bind_rows(z.FF0, z.TF0, z.TT0, z.TT1, z.TT2)
# rm(z.FF0, z.TF0, z.TT0, z.TT1, z.TT2)

# Get the max and min z-values by variable and by trial
plot.data <- z.jsdm %>%
  group_by(Taxon, Variable, Trial) %>%
  summarise(z.min = min(z, na.rm=T), z.max = max(z, na.rm=T)) %>%
  ungroup() %>%
  mutate(z.range = z.max-z.min,
         # Variable = factor(Variable, levels=c("Temp", "FV", "F100m", "IAR", "LUD", "Site effect", "TT1", "TT2")),
         Trial = factor(Trial, levels=c("UF0", "CF0", "CT0", "CT1", "CT2"))
         )

z.jsdm %>% filter(Trial=="CT0") %>% count(Variable)
test <- filter(z.jsdm, Variable=="Site effect", Trial=="CT0")
test <- filter(plot.data, Variable=="Site effect", Trial=="CT0")

g <- ggplot(data = plot.data) 
g <- g + geom_boxplot(aes(y=z.range, fill=Variable))
g <- g + facet_grid(. ~ Trial)
g <- g + coord_cartesian(ylim=c(0,16)) 
g <- g + theme_bw(base_size = 24)
g <- g + labs(title="Range of linear predictor for all taxa",
              y = expression(paste("z"["range"])),
              fill = "Variable")
g <- g + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank())
g <- g + scale_fill_brewer(palette = "Set1")
g

pdf('outputs/paper 1 extensions/P1 - slopes.pdf', width = 15)
g
dev.off()

# > i/jSDM cross-validation ####

# RESULTS ####
# Plot P1 - link function ####
pdf('P1 - link function.pdf',height=6, width=8)
par(cex=1.5, mar=c(5,4,2,2)+0.1)
plot(x=seq(from=-4, to=4, by=0.25), y=link.function(seq(from=-4, to=4, by=0.25)), type="l", yaxs="i",
     ylab="P(Y=1)",
     xlab="z")
dev.off()

par(cex=3)
plot(x=seq(from=-4, to=4, by=0.25), y=link.function(seq(from=-4, to=4, by=0.25)), type="l", lwd=2, yaxs="i",
     ylab="P(Y=1)",
     xlab="z")
axis(side = 1, lwd = 2)
axis(side = 2, lwd = 2)
box(lwd=2)

# Plot P1 - BDM occurrence frequency ####
pdf('P1 - BDM prevalence.pdf', height = 10, width = 13)
par(cex=2.5)
plot(x = 1:length(n.bdms), y = n.bdms, axes=F, 
     xlab="Taxa (ordered by decreasing prevalence)", 
     ylab="Prevalence (number of presence data points)")
axis(1, c(0,50,100,150,200,245))
axis(2, c(0,100,200,300,400,500,580))
dev.off()





# Plot P1 - cross-validation ####
# pre-process plotting data for ggplot
cv.plot <- deviance.cv %>%
  select(Taxon, Type, Model, std.deviance) %>%
  group_by(Taxon, Type, Model) %>%
  summarise(MRD = mean(std.deviance)) %>%
  spread(Type, MRD) %>%
  ungroup() %>% # remove grouping information
  mutate(n = n.bdms[Taxon]) # add occurrence frequency again

# ymax <- max(cv.plot$Testing[!is.infinite(cv.plot$Testing)])
ymax <- max(cv.plot$Training)

# cv.plot$Infinite <- is.infinite(cv.plot$Testing)
cv.plot$ymax <- ifelse(cv.plot$Testing > ymax, TRUE, FALSE)

cv.plot$Shape <- ifelse(cv.plot$ymax, 24, 16)
cv.plot$Testing <- ifelse(cv.plot$ymax, max(cv.plot$Training), cv.plot$Testing)
# # Set infinite deviances values to the plot edge
# cv.plot$Testing <- ifelse(is.infinite(cv.plot$Testing), max(cv.plot$Testing[!is.infinite(cv.plot$Testing)]), cv.plot$Testing)
# 
# # Format shape based on is.infinite deviance
# cv.plot$Shape <- ifelse(is.infinite(cv.plot$Testing), 24, 16)
# 
# cv.plot$n <- ifelse(cv.plot$Infinite, NA, cv.plot$n)
cv.plot$Shape <- as.factor(cv.plot$Shape)
cv.plot <- filter(cv.plot, !(Taxon %in% c("Silo", "Stactobia")))

# Plot data
g <- ggplot(cv.plot, aes(x = Training, y = Testing, colour = Model, size = n, shape = Shape))
g <- g + geom_point(alpha = 0.3)
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
g <- g + theme(strip.background=element_rect(fill="black"), strip.text=element_text(color="white", face="bold"), axis.text=element_text(size=18))
g <- g + theme_bw(base_size = 18)
g <- g + guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))
g <- g + labs(x = "Mean standardized deviance for calibration",
              y = "Mean standardized deviance for prediction",
              colour = "Model",
              shape = "Deviance",
              size = "Prevalence")
g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
g <- g + scale_size_continuous(range = c(2, 8))
g <- g + scale_shape_discrete(name  = "Deviance",
                              breaks=c("16", "24"),
                              labels=c("In range", "Out of range"))
g <- g + coord_cartesian(ylim=c(0, ymax))
print(g)

g <- ggplot(cv.plot)
g <- g + geom_point(aes(x=n, y=Testing, color=Model, shape=Shape, size = n))
g <- g + theme(strip.background=element_rect(fill="black"), strip.text=element_text(color="white", face="bold"), axis.text=element_text(size=18))
g <- g + theme_bw(base_size = 18)
g <- g + guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))
g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
g <- g + scale_size_continuous(range = c(2, 8))
g <- g + scale_shape_discrete(name  = "Deviance",
                              breaks=c("16", "24"),
                              labels=c("In range", "Out of range"))
g <- g + coord_cartesian(ylim=c(0, ymax))
print(g)

# pdf('P1 - cross-validation.pdf', height=8, width=10.5)
# print(g)
# dev.off()

# Plot P1 - slopes jSDM ####
slopes <- linear.predictor(jsdm.p1)

# Match expressions to influence factors
slopes.plot <- slopes
v <- jsdm.p1$bdms$inf.fact[jsdm.p1$bdms$inf.fact != "Temp2"]
slopes.plot$Label <- factor(slopes.plot$Variable, levels = v)
levels(slopes.plot$Label) <- labeller(levels(slopes.plot$Label))

# Create data to plot influence factors
x <- prepare.inputs(v, sample.bdms, center = FALSE)
x <- gather(x, Variable, Value, -SiteId, -SampId)
x <- filter(x, Variable != "Temp2")
x$Label <- factor(x$Variable, levels = v)
levels(x$Label) <- labeller(levels(x$Label))

g <- ggplot(data=slopes.plot)
g <- g + geom_line(data = slopes.plot, aes(x = x, y = z, group = Taxon, alpha = alpha))
g <- g + geom_rug(data = x, aes(x = Value, y = -15), sides="b")
g <- g + facet_wrap(~ Label, scales = "free_x", labeller=label_parsed, strip.position = "bottom")
g <- g + theme_bw(base_size = 18)
g <- g + theme(axis.title.y = element_text(size = 22), strip.text.x = element_text(size = 17))
g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
               strip.background = element_blank(), strip.placement = "outside")
g <- g + labs(x = "", y = expression(z[it[i]*jk]))
g <- g + guides(alpha = FALSE)
print(g)

pdf('P1 - slopes jSDM.pdf', width = 13, height = 9, onefile=TRUE)
print(g)
dev.off()

# Plot P1 - plot.comm() ####
plot.comm(jsdm.p1, 'P1 - plot_comm [all taxa]')

# Plot P1 - plot.comm() [example taxa] ####
# Modified plot.comm() code for specific taxa.
pdf('P1 - plot_comm [example taxa].pdf', onefile = TRUE)
par(cex=1.25)
for (k in 1:length(inf.fact)){
  variable <- inf.fact[k]
  responses <- response.bdms[response.bdms$Taxon %in% t, ]
  
  response <- responses[[variable]]
  names(response) <- responses$Taxon
  
  response[response==0] <- "grey55"
  response[response==1] <- "blue"
  response[response==-1] <- "red"
  cat("Plotting: ",variable,"\n")
  # Test temperature as starting variable
  samples <- taxon.samples[Variable == variable,]
  
  # Find the maximum density among the posterior taxon-specific posterior parameters
  samples.sd <- samples %>%
    group_by(Taxon) %>%
    summarise(SD = sd(Value)) %>%
    arrange(SD) # arrange from min to max
  
  sd <- samples.sd[1,]$Taxon
  ymax.sample <- samples[Taxon == sd, ]
  ymax <- density(ymax.sample$Value)
  ymax <- max(ymax$y)
  
  # Plot the community parameter distribution
  mu <- bdms.stan$mu.beta.comm.maxpost[variable]
  sd <- bdms.stan$sigma.beta.comm.maxpost[variable]
  x <- seq(mu-4*sd,mu+4*sd,length.out=201)
  beta.comm.density <- dnorm(x, mu, sd)
  beta.taxa.maxpost.density <- density(bdms.stan$beta.taxa.maxpost[variable, ])
  
  # Match expressions to influence factors
  labels <- c("Temp" = expression(paste(beta["Temp"], " (1/", degree, "C)")),
              "Temp2" = expression(paste(beta["Temp"^2], " (1/", degree, "C"^2,")")),
              "FV" = expression(paste(beta["FV"], " (s/m)")),
              "F10m" = expression(paste(beta["F10m"], " (1/%)")),
              "IAR" = expression(paste(beta["IAR"], " 1/(spray treatments * fraction cropland)")),
              "Urban" = expression(paste(beta["Urban"], " (1/%)")),
              "LUD" = expression(paste(beta["LUD"], " (km"^2,"/CE)"))
  )
  
  plot(numeric(0), numeric(0),
       xlim = c(min(x), max(x)),
       ylim=c(0, ymax-(ymax*0.25)),
       xlab = labels[variable],
       ylab = "Density")
  abline(v=0)
  # Plot the taxon-specific parameters
  t <- c("Gammaridae", "Nemoura_minima", "Protonemura_lateralis")
  l <- c("G", "N", "P")
  
  for (j in 1:length(t)){
    taxon <- t[j]
    sample <- samples[Taxon==taxon,]
    sample <- sample$Value
    
    lines(density(sample), type="l", col=alpha(response[taxon], 1), lwd=1.25)
    text(x=mean(sample), y=max(density(sample)$y)+0.1*max(density(sample)$y), labels = l[j])
  }
  # Plot community parameter distribution
  lines(x, beta.comm.density, type="l", col = "grey50", lwd = 5)
  # Plot maximum posterior values over all taxa
  lines(beta.taxa.maxpost.density, type="l", col="black", lwd=2, lty='longdash')
  
}
dev.off()


# Plot P1 - significant responses ####
bdms.stan <- select.jsdm(jsdm.p1) # retrieve Stan.fit object
inf.fact <- bdms.stan$inf.fact
beta.samples <- bdms.stan[["beta_taxa"]]
dimnames(beta.samples) <- list(1:dim(beta.samples)[1], inf.fact, colnames(bdms.stan$occur.taxa)) # name the dimensions
taxon.samples <- melt(beta.samples)
colnames(taxon.samples) <- c("Sample", "Variable", "Taxon", "Value")
taxon.samples <- setDT(taxon.samples)
taxon.samples$Variable <- as.character(taxon.samples$Variable)
taxon.samples$Taxon <- as.character(taxon.samples$Taxon)

response.bdms <- extract.resp(jsdm.p1)
# Based on quality-of-fit and predictive performance, find taxa of interest that:
temp <- cv.plot %>%
  filter(Model=="jSDM") %>%
  select(Taxon, Training, Testing)

# Combine goodness-of-fit and predictive performance
taxa.int <- dev %>%
  filter(Model=='jSDM', D2 > 0.25 & D2 < 0.60, std.deviance > 0.15 & std.deviance < 0.90) %>%
  select(Taxon, D2) %>%
  left_join(temp, by = "Taxon") %>%
  mutate(n=n.bdms[Taxon], TrTe = Training/Testing)
rm(temp)

# Select taxa for heatmap
d.select <- spread(select(dev, Taxon, Model, D2), Model, D2)
colnames(d.select) <- c("Taxon", "iSDM.D2", "jSDM.D2")

std.dev.select <- spread(select(dev, Taxon, Model, std.deviance), Model, std.deviance)
colnames(std.dev.select) <- c("Taxon", "iSDM.dev", "jSDM.dev")

dt <- left_join(d.select, std.dev.select, by = "Taxon")
setDT(dt)
rm(d.select, std.dev.select)

# # joint model will always have higher residual deviance; what % higher?
# dt$rdev.diff <- (dt$iSDM.dev/dt$jSDM.dev)*100
# 
# # joint model will always have lower D2; what % 
# dt$d.diff <- (dt$iSDM.D2/dt$jSDM.D2)*100

# selection criteria;
# good fit in joint model,
# iSDM deviance >95% of jSDM deviance,
# iSDM D2 not more than 10% of jSDM D2
# nrow(filter(dt, jSDM.D2 > 0.30 & rdev.diff > 0.95 & d.diff < 110))

# Format density plot output for heatmap

response.bdms$n <- n.bdms[response.bdms$Taxon]
colnames(response.bdms) <- c(jsdm.p1$bdms$inf.fact, "Taxon", "n")
response.bdms <- arrange(response.bdms, desc(n))


# Color coding for all taxa
# Subset the taxa with good model performance
btr.int <- filter(response.bdms, Taxon %in% taxa.int$Taxon)
btr.int <- btr.int[, c("Taxon", "n", inf.fact)]
btr.plot <- btr.int[btr.int$n > 56, c(inf.fact)]
btr.plot <- apply(btr.plot, 2, as.numeric)
row.names(btr.plot) <- paste(sub("_", " ", btr.int$Taxon[btr.int$n > 56]), " - ", n.bdms[btr.int$Taxon[btr.int$n > 56]])
btr.plot <- btr.plot[, c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD")]

# Select the number of clusters
d <- dist(btr.plot, method="euclidean")
pfit <- hclust(d, method="ward.D")
plot(pfit, labels=row.names(btr.plot))
rect.hclust(pfit, k=4)

library(pheatmap)
# Save manually
pheatmap(btr.plot, col=c("tomato3", "snow3", "royalblue4"), cellwidth=30, cellheight=11, cluster_rows=T, cutree_rows=5, cluster_cols=F, clustering_distance_rows = "euclidean", legend=F, show_rownames=TRUE, fontsize=13, fontsize_col=15, labels_col = labeller.names(c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD")))


# Plot P1 - significant responses [all taxa] ####
btr.all <- response.bdms
btr.all <- as.data.table(btr.all)
btr.all <- arrange(btr.all, desc(n))
btr.all$TaxonLabel <- factor(paste(btr.all$Taxon, " - ", btr.all$n), levels=paste(names(sort(n.bdms)), " - ", sort(n.bdms)))

# plot.data$Labels <- factor(paste(plot.data$Taxon, ' - ', n.bdms[plot.data$Taxon]), levels = paste(names(sort(n.bdms)), ' - ', sort(n.bdms)))

plot.data <- gather(btr.all, Variable, Value, -Taxon, -TaxonLabel, -n)
plot.data$VariableLabel <- factor(plot.data$Variable, levels=c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD"))

g <- ggplot(plot.data, aes(x = VariableLabel, y = TaxonLabel))
g <- g + geom_tile(aes(fill = as.factor(Value)), colour = "white")
g <- g + scale_fill_manual(values=c("tomato3", "snow3", "royalblue4"))
g <- g + theme_minimal(base_size = 15)
g <- g + labs(fill = "Response",
              y = "Taxon and prevalence",
              x = "Explanatory Variable")
pdf('jSDM heatmap ALL.pdf', paper='special', height=56, width=8.5)
print(g)
dev.off()



# Plot P1 - map.ijSDM() ####
map.isdm(isdm.p1, 'P1 - maps iSDM')
map.jsdm(jsdm.p1, 'P1 - maps jSDM')

# Plot P1 - plot.prob() [labelled] ####
plot.prob(isdm.p1, 'P1 - prob vs inputs iSDM')
plot.prob(jsdm.p1, 'P1 - prob vs inputs jSDM')

# Plot map.jsdm.pred.taxon() examples ####
# Grid arrange prob.taxon and map.taxon plots
g1 <- plot.prob.taxon(jsdm.p1, "Gammaridae", legend=FALSE)
g2 <- map.jsdm.taxon(jsdm.p1, "Gammaridae", legend=FALSE)
g3 <- map.jsdm.taxon(jsdm.p1, "Nemoura_minima", legend=FALSE)
g4 <- map.jsdm.taxon(jsdm.p1, "Protonemura_lateralis", legend=FALSE)

plot_grid(g1,g2,g3,g4, labels=c("(a)", "(b)", "(c)", "(d)"), align="hv")

# Grid arrange prob.taxon and map.taxon plots
g1 <- plot.prob.taxon(jsdm.p1, "Gammaridae", legend=FALSE)
g2 <- map.jsdm.pred.taxon(jsdm.p1, "Gammaridae", legend=FALSE)
g3 <- map.jsdm.pred.taxon(jsdm.p1, "Nemoura_minima", legend=FALSE)
g4 <- map.jsdm.pred.taxon(jsdm.p1, "Protonemura_lateralis", legend=FALSE)

pdf('P1 - example taxa.pdf', width=13, height=9.5)
plot_grid(g1,g2,g3,g4, labels=c("(a)", "(b)", "(c)", "(d)"), align="hv")
dev.off()

# pdf('P1 - map jSDM [Gammaridae].pdf', width=13, height=9.5)
g <- map.jsdm.pred.taxon(jsdm.p1, "Gammaridae", legend=TRUE)
# SVG produced to transfer legend to grid-arrange plot above.
ggsave(file="P1 - map jSDM [Gammaridae].svg", plot=g, width=13, height=9.50)

# dev.off()
# pdf('P1 - example taxa [gammaridae map].pdf', width=13, height=9.5)
# g2 <- fit.jsdm.pred.taxon(jsdm.p1, "Gammaridae", legend=TRUE)
# g2
# dev.off()


# library(gridExtra)
# library(grid)
# library(ggplot2)
# library(lattice)
# grid.arrange(g1,g2,g3,g4, nrow=2)


# library(gtable)
# g1 <- ggplotGrob(g1)
# g2 <- ggplotGrob(g2)
# g3 <- ggplotGrob(g3)
# g4 <- ggplotGrob(g4)
# g <- rbind(g1, g2, g3, g4)
# g$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
# grid.newpage()
# grid.draw(g)

# Plot P1 - maps influence factors ####
x <- prepare.inputs(K, sample.bdms, center = FALSE)
x <- gather(x, Variable, Value, -SiteId, -SampId)
x <- filter(x, Variable != "Temp2")
x$Trial <- 'BDM species'
map.inputs(x, 'P1 - maps influence factors')



# Plot P1 - parameter SD ####
# Extract the entire posterior beta sample
jsdm.beta.samples <- extract.beta(jsdm.p1)

# Fast aggregation of 10M row dataset
plot.data <- jsdm.beta.samples %>%
  group_by(Taxon, Variable) %>%
  summarise(SD = sd(Value), Mean = mean(Value)) %>%
  mutate(n = n.bdms[Taxon], Label = factor(Variable, levels=c("Temp", "Temp2", "FV", "F10m",  "IAR",   "Urban", "LUD")), rSD = SD/abs(Mean))

levels(plot.data$Label) <- labeller(levels(plot.data$Label))
setDT(plot.data)
# g <- ggplot(data=plot.data, aes(x = n, y = SD, size = n))
# g <- g + geom_point(alpha = 0.5)
# g <- g + facet_grid(Label ~ ., scales = "free", labeller=label_parsed)
# g <- g + theme_bw(base_size = 14)
# g <- g + theme(strip.background=element_rect(fill="grey"),strip.text=element_text(color="black", face="bold"),
#                plot.title = element_text(hjust = 0.5, size = 12))
# g <- g + labs(title = expression(paste("Standard deviation of posterior taxon-specific parameter distributions ", beta["jk"]^"taxa")),
#               x = "Prevalence",
#               y = expression(paste("Standard deviation (", sigma[beta["jk"]^"taxa"],")")),
#               size = "Prevalence")
# g <- g + scale_y_continuous(limits=c(0,NA))
# pdf('P1 - parameter SD.pdf', height = 12.5)
# print(g)
# dev.off()

# # Plot P1 - parameter uncertainty ####
# g <- ggplot(data=plot.data, aes(x = n, y = rSD, size = n))
# g <- g + geom_point(alpha = 0.5)
# g <- g + facet_wrap(~ Variable, labeller=label_parsed)
# g <- g + theme_bw(base_size = 13)
# g <- g + theme(strip.background=element_rect(fill="black"), strip.text=element_text(color="white", face="bold"), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g <- g + labs(title = expression(paste(sigma,"/",mu)),
#               x = expression(paste("Mean (", mu[beta["jk"]^"taxa"],")")),
#               y = expression(paste("Standard deviation/mean (", sigma[beta["jk"]^"taxa"],"/", mu[beta["jk"]^"taxa"],"))")),
#               size = "Total occurrences")
# print(g)

plot.data$Variable <- factor(plot.data$Variable, levels=inf.fact)
# levels(plot.data$Variable) <- labeller(levels(plot.data$Variable))

g <- ggplot(data=plot.data)
g <- g + geom_boxplot(aes(x=Variable, y=rSD, fill=Variable))
g <- g + coord_cartesian(ylim=c(0,10))
g <- g + theme_bw()
g <- g + theme(plot.title = element_text(hjust = 0.5, size = 12), 
               axis.text = element_text(size = 12),
               axis.text.x = element_text(angle=45, hjust=1, vjust=0.5))
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + guides(fill=FALSE)
g <- g + labs(title = expression(paste("Relative uncertainty of posterior taxon-specific parameter distributions ", beta["jk"]^"taxa")),
              x = "Explanatory variable",
              y = expression(paste("Relative standard deviation (", sigma[beta["jk"]^"taxa"]," / |",mu[beta["jk"]^"taxa"],"|)")))
pdf('P1 - parameter uncertainty.pdf')
print(g)
dev.off()

# Plot P1 - predictive uncertainty ####
pred <- cv.jsdm.pred("bdms")
gc()

# what about binding the probabilities manually?
taxon <- "Baetis_alpinus"
f1 <- pred$fold1[Taxon=="taxon",]
f2 <- pred$fold2[Taxon=="taxon",]
f3 <- pred$fold3[Taxon=="taxon",]

dt <- bind_rows(f1,f2,f3)
dt <- left_join(dt, sample.bdms[, c("SampId", taxon)], by="SampId")
colnames(dt) <- c("Taxon", "SampId", "Pred", "Obs")

rm(f1,f2,f3)

# Order samples by decreasing mean predicted probability
test <- dt %>%
  group_by(SampId) %>%
  summarise(mean.Pred = mean(Pred)) %>%
  arrange(-mean.Pred)


plot.data <- na.omit(dt)
plot.data$SampId <- factor(plot.data$SampId, levels=test$SampId)
plot.data$Obs <- as.factor(plot.data$Obs)
rm(test)

g <- ggplot(plot.data)
g <- g + stat_density_ridges(aes(x = Pred, y = SampId, fill = Obs), size = 10, color=NA, rel_min_height=0.01, scale=10, alpha = 0.5)
# g <- g + theme_bw(base_size=17)
g <- g + theme_pubr()
g <- g + theme(#plot.title=element_blank(),
  axis.ticks.x = element_blank(),
  axis.text.x = element_blank(),
  panel.background = element_rect(color = "white"))

g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(title="Predictions from k-fold cross-validation",
              x="Posterior probability",
              y="Sample sites (ordered by mean probability)",
              fill="Observation")
g <- g + scale_fill_manual(values= c("1" = "blue", "0" = "red"), labels=c("Presence", "Absence"))
g <- g + coord_flip()
pdf('Predictive uncertainty.pdf', width = 10, height = 7)
g
dev.off()

# Plot P1 - map predictive uncertainty ####
ch <- fortify(inputs$ch)

dt <- extract.jsdm.pred(jsdm.p1, get.quantiles=TRUE) # Get predicted probabilities with default quantiles c(0.05, 0.95)
dt <- left_join(dt, inputs$xy, by="SiteId")
setDT(dt)

# Format data for ggplot() aesthetics
dt <- na.omit(dt)
dt$Obs <- as.factor(dt$Obs)
dt$Alpha <- ifelse(dt$Quantile==0.05, 0.65, 0.35)
dt$Alpha <- as.factor(dt$Alpha)
dt$Shape <- ifelse(dt$Quantile==0.05, 19, 21)
dt$Stroke <- ifelse(dt$Quantile==0.05, 0, 0.75)

taxa <- occur.freq(jsdm.p1$bdms$occur.taxa)

pdf(paste("P1 - maps jSDM.pdf", sep=''), paper = 'special', width = 10.5, onefile = TRUE)
for (j in 1:length(taxa)){
  taxon <- names(taxa[j])
  
  plot.data <- dt[Taxon==taxon, ]
  
  # Map geometries
  g <- ggplot()
  g <- g + geom_polygon(data = ch, aes(x = long, y = lat, group = group), fill=NA, color="black")
  g <- g + geom_point(data = plot.data, aes(X, Y, size = Pred, alpha = Alpha, color = Obs, stroke = Stroke, shape = Shape))
  
  # Configure themes and labels
  g <- g + labs(title = paste("Probability of occurrence vs observations of", taxon),
                subtitle = paste("jSDM:", paste(jsdm.p1$bdms$inf.fact, collapse = "+", sep = " "), "- page", j),
                x = "",
                y = "",
                size = "Probability of\noccurrence",
                alpha = "Posterior",
                color = "Observation")
  g <- g + theme_minimal(base_size = 15)
  g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  # Configure legends and scales
  g <- g + guides(size = guide_legend(override.aes = list(color="black", stroke=0), order=1),
                  alpha = guide_legend(override.aes = list(size=6, shape=c(19,21), stroke=c(0,0.75), color="black"), order=2),
                  color = guide_legend(override.aes = list(size=6, stroke=0), order=3))
  g <- g + scale_y_continuous(breaks=NULL)
  g <- g + scale_x_continuous(breaks=NULL)
  g <- g + scale_radius(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 6))
  g <- g + scale_color_manual(values=c("0" = "#FF0000", "1" = "#0077FF"), labels=c("Absence", "Presence"))
  g <- g + scale_alpha_manual(values=c("0.65"="0.65", "0.35"="0.35"), labels=c("5th quantile", "95th quantile"))
  g <- g + scale_shape_identity()

  cat("Plotting taxon: ", taxon, "\n")
  print(g)
}
dev.off()

# Table P1 - Count significant responses ####
# Count significant parameters
p1.uncertainty <- jsdm.beta.samples %>%
  group_by(Variable, Taxon) %>%
  do(data.frame(t(quantile(.$Value, probs = c(0.05, 0.95))))) %>%
  rename(quant5 = X5., quant95 = X95.) %>%
  setDT()

p1.uncertainty$response <- ifelse(p1.uncertainty$quant5 > 0, "positive", "")
p1.uncertainty$response <- ifelse(p1.uncertainty$quant95 < 0, "negative", p1.uncertainty$response)
p1.uncertainty$response <- ifelse(!(p1.uncertainty$quant5 > 0) & !(p1.uncertainty$quant95 < 0), "non-significant", p1.uncertainty$response)

p1.uncertainty <- p1.uncertainty %>%
  group_by(Variable, response) %>%
  summarise(n = n()) %>%
  ungroup()

p1.uncertainty <- spread(p1.uncertainty, response, n)
p1.uncertainty[is.na(p1.uncertainty)] <- 0

write.csv(p1.uncertainty, 'jSDM - significant responses.csv', row.names = F)


pdf('outputs/jSDM SD vs mean.pdf', paper='a4r', onefile = TRUE)
print(g)
dev.off()

# Plot P1 - map epsilon ####
epsilon <- data.table(SiteId = names(jsdm.re$eps.maxpost), epsilon.maxpost = jsdm.re$eps.maxpost)
epsilon <- left_join(epsilon, inputs$xy, by="SiteId")
epsilon$color = ifelse(epsilon$epsilon.maxpost > 0, "Positive", "Negative")

ch <- fortify(inputs$ch)
g <- ggplot()
g <- g + geom_polygon(data = ch, aes(x = long, y = lat, group = group), fill=NA, color="black")
g <- g + geom_point(data = epsilon, aes(X, Y, size = abs(epsilon.maxpost), color = color), alpha = 0.35)
g <- g + labs(title = "",
              subtitle = "",
              x = "",
              y = "")
g <- g + scale_y_continuous(breaks=NULL)
g <- g + scale_x_continuous(breaks=NULL)
# g <- g + scale_size_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 7))
g <- g + scale_size_continuous(breaks=seq(0, 3.5, 0.5), limits=c(0,3.5), range = c(1,8))
g <- g + scale_color_manual(values= c("Positive" = "blue", "Negative" = "red"))
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + labs(size = expression(paste(epsilon["i"], " (absolute value)")),
              color = "")
g <- g + theme_minimal(base_size = 15)
print(g)

pdf('P1 - map epsilon.pdf', paper = 'special', width = 10.5, onefile = TRUE)
print(g)
dev.off()

# Plot P1 pairwise scatterplot epsilon ####
x <- prepare.inputs(c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD"), sample.bdms, center=FALSE)
x$Temp2 <- (p1.predictors$Temp-mean(x$Temp))^2
epsilon.plot <- left_join(epsilon, x,by="SiteId")
epsilon.plot <- select(epsilon.plot, Temp, Temp2, FV, F10m, IAR, Urban, LUD, epsilon.maxpost)
colnames(epsilon.plot)[8] <- "RandomEffect"
setDT(epsilon.plot)
pairs.panels(epsilon.plot, density = TRUE, scale=FALSE, hist.col="grey", cex.cor=1.5, cex.labels=1.5)


# Plot P1 - parameter dotplot [4 pages] ####
# Warning: really bad code ahead...
# Extract the entire posterior beta sample
jsdm.beta.samples <- extract.beta(jsdm.p1)

jsdm.beta.samples.quantiles <- jsdm.beta.samples %>%
  group_by(Variable, Taxon) %>%
  summarise(quantile10=quantile(Value, 0.10), quantile90=quantile(Value, 0.90)) %>%
  setDT()
colnames(jsdm.beta.samples.quantiles)<- c("Variable", "Taxon", "quantile10", "quantile90")

isdm.beta <- isdm.p1$parameters
if ("Model" %in% colnames(isdm.beta)){
  isdm.beta$Model <- NULL
}
colnames(isdm.beta) <- c("Taxon", "Variable", "iSDM.parameter")
isdm.beta$Variable <- as.character(isdm.beta$Variable)

# Number of pages per taxa (ntpp)
pages <- 4
ntpp <- length(n.bdms)/pages
tpp <- split(n.bdms, ceiling(seq_along(n.bdms)/ntpp)) # bin the number of taxa per page
pdf(paste("outputs/jsdm_p1/parameter_dotplot_notnorm.pdf",sep=""), paper = "special", height = 10, width = 9, onefile = TRUE)
for (page in 1:pages){
  
  plot.data <- jsdm.beta.samples.quantiles[jsdm.beta.samples.quantiles$Taxon %in% names(tpp[page][[1]]), ]
  
  # Prepare the jSDM parameter values
  beta.max <- select(jsdm.p1$parameters, Taxon, Variable, Parameter)
  colnames(beta.max) <- c("Taxon", "Variable", "jSDM.parameter")
  beta.max$Variable <- as.character(beta.max$Variable)
  
  # Join the quantiles, jSDM, and iSDM parameters together
  plot.data <- left_join(plot.data, beta.max, by = c("Taxon", "Variable"))
  plot.data <- left_join(plot.data, isdm.beta, by = c("Taxon", "Variable"))
  
  # Model (variable) and Parameter (value) gather the two columns, jSDM.parameter and iSDM.parameter
  # into a variable/value pair
  plot.data <- gather(plot.data, Model, Parameter, jSDM.parameter:iSDM.parameter)
  
  # Process iSDM outliers
  # Get min and max posterior parameter values within the 10th-90th quantiles
  lim.min <- aggregate(quantile10 ~ Variable, plot.data, min)
  vmin <- lim.min$quantile10; names(vmin) <- lim.min$Variable
  
  lim.max <- aggregate(quantile90 ~ Variable, plot.data, max)
  vmax <- lim.max$quantile90; names(vmax) <- lim.max$Variable
  
  # Set colors for each parameter
  # g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
  plot.data$Col <- "#048504"
  plot.data$Col[plot.data$Model == 'iSDM.parameter'] <- "#790FBF"
  plot.data$Col <- ifelse(plot.data$Parameter < vmin[plot.data$Variable], "#6B6B6B", plot.data$Col)
  plot.data$Parameter <- ifelse(plot.data$Parameter < vmin[plot.data$Variable], vmin[plot.data$Variable], plot.data$Parameter)
  
  plot.data$Col <- ifelse(plot.data$Parameter > vmax[plot.data$Variable], "#6B6B6B", plot.data$Col)
  plot.data$Parameter <- ifelse(plot.data$Parameter > vmax[plot.data$Variable], vmin[plot.data$Variable], plot.data$Parameter)
  rm(lim.min, lim.max, vmin, vmax)
  
  # Order labels by occurrence frequency
  plot.data$Labels <- factor(paste(plot.data$Taxon, ' - ', n.bdms[plot.data$Taxon]), levels = paste(names(sort(n.bdms)), ' - ', sort(n.bdms)))
  
  # Order variable facets and pass expressions for units in facet labels
  plot.data$Variable <- factor(plot.data$Variable, levels = c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD"))
  levels(plot.data$Variable) <- labeller(c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD"))

  # Build the plot
  g <- ggplot(plot.data)
  g <- g + geom_hline(yintercept=0, alpha=0.4)
  g <- g + geom_linerange(aes(x = Labels, ymin = quantile10, ymax = quantile90), color = 'black', alpha = 0.6)
  g <- g + geom_point(aes(x = Labels, y = Parameter, colour = Col), stroke=0, size = 2.5, alpha = 0.5)
  g <- g + facet_grid(. ~ Variable, scales = "free", labeller=label_parsed)
  g <- g + coord_flip()
  g <- g + theme_bw()
  g <- g + labs(title = "Taxon-specific parameter estimates by model",
                subtitle = "MLEs (iSDM) and maximum posterior (jSDM) with 10th-90th percentile interval",
                x = "Taxon and prevalance",
                y = "Value")
  g <- g + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  g <- g + scale_colour_identity()
  print(g)
}
dev.off() 

# > Plots: k-fold CV PROB ###
# Plot overall probability ###
glm.prob$Model <- "iSDM"
glm.prob <- glm.prob[, colnames(j$probability)]

prob <- rbind(j$probability, glm.prob)
prob <- na.omit(prob)
prob$Obs <- ifelse(prob$Obs == 1, "Present", "Absent")
prob$Obs <- as.factor(prob$Obs)
prob$Type <- factor(prob$Type, levels = c("Training", "Testing"))

g <- ggplot(prob, aes(x=Obs, y = Pred))
g <- g + geom_boxplot(aes(fill=Obs))
# g <- g + stat_boxplot(geom = "errorbar", width = 0.25)
g <- g + facet_grid(. ~ Model + Fold + Type)
g <- g + labs(x="")
g <- g + theme_bw()
g <- g + theme(strip.background=element_rect(fill="black"))
g <- g + theme(strip.text=element_text(color="white", face="bold"))
g <- g + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
g <- g + theme(plot.title = element_text(hjust = 0.5))
g <- g + labs(title = "Training/testing probabilities for jSDM vs iSDM (50-50)",
              x = "Observed",
              y = "Predicted Probability")
g <- g + scale_fill_brewer(palette = "Set1")
print(g)
ggsave('../CV - D1 vs GLM probability.pdf', device = 'pdf', width = 11, height = 8.5, units = "in")

# Plot probability ~ species ###
page <- 1
pdf(paste('../CV - D1 vs GLM probability by species', ".pdf", sep=''), paper = 'special', width = 14, height = 7, onefile = TRUE)
prob$Type <- factor(prob$Type, levels = c("Training", "Testing"))

# Get unique taxa that occur in each fold
spf <- aggregate(Species ~ Fold, prob[prob$Type == 'Training', ], unique)
# Get taxa that intersect the folds
# taxa <- intersect(spf$Species[[1]], intersect(spf$Species[[2]], spf$Species[[3]]))
taxa <- sapply(spf$Species, function(x){unique(x)})
taxa <- unique(unlist(taxa))
taxa <- n.bdms[taxa]
taxa <- sort(taxa, decreasing = T)
taxa <- names(taxa)
rm(spf)
for (i in 1:length(taxa)){
  tname <- taxa[i]
  pdata <- prob[prob$Species == tname, ]
  pdata$Trial <- factor(pdata$Trial, levels = c("iSDM", "jSDM"))
  pdata$Fold <- factor(pdata$Fold)
  g <- ggplot(pdata, aes(x=Obs, y=Pred, fill=Fold))
  g <- g + geom_boxplot()
  g <- g + facet_grid(. ~ Trial + Type, scales = "free")
  
  g <- g + theme_bw(base_size = 14)
  g <- g + ggtitle(paste("Page ", page, " - probability for training/testing of iSDM vs jSDM for ", paste(tname), ' (n = ', n.bdms[tname], ')', sep = ''))
  g <- g + scale_fill_brewer(palette = "Set1")
  
  print(g)
  page <- page + 1
  cat("Plotting taxa: ", tname, "\n")
}
dev.off()


# Plot probability vs influence factor ####
# For each taxon, plot the probabilities against each influence factor
# First dataset needed: probabilities

# Second dataset needed: predictors (raw data)
pdf(paste('outputs/P1 - prob iSDM.pdf', sep=''), paper = 'special', height = 9, width = 12, onefile = TRUE)

d <- sapply(isdm.bdms$models, function(m){
  D2(m)
})

for (j in 1:length(names(isdm.bdms$taxa))){
  taxon <- names(isdm.bdms$taxa)[j]
  d.squared <- round(d[j], 2)
  
  taxon.prob <- filter(isdm.bdms$, Species == taxon)
  influence.factors <- gather(x, Variable, Value, -SiteId, -SampId)
  
  # Create plot data
  plot.data <- left_join(taxon.prob, influence.factors, by = c("SiteId", "SampId"))
  plot.data <- filter(plot.data, Variable != "Temp2")
  plot.data$Obs <- as.factor(plot.data$Obs)
  plot.data$Label <- factor(plot.data$Variable, levels = c("Temp", "FV", "F10m", "IAR", "Urban", "LUD"))
  levels(plot.data$Label) <- c(expression(paste("Temp (", degree, "C)")),
                               expression(paste("FV (m/s)")),
                               expression(paste("F10m (%)")),
                               expression(paste("IAR (spray treatments * fraction cropland)")),
                               expression(paste("Urban (%)")),
                               expression(paste("LUD (CE/km"^2,")")))
  
  g <- ggplot(data = plot.data, aes(x = Value, y = Pred, color = Obs))
  g <- g + geom_point(alpha = 0.25)
  g <- g + theme_bw(base_size=15)
  g <- g + facet_wrap(~ Label, scales = "free_x", labeller=label_parsed, strip.position="bottom")
  g <- g + labs(title = paste(pg, " - ", taxon),
                x = "Influence factor",
                y = "Probability",
                color = "Observation")
  g <- g + theme(strip.background = element_blank(), strip.placement = "outside")
  g <- g + scale_color_manual(name = "Observation", values=c("#FF0000", "#0077FF"), labels=c("Absence", "Presence"))
  g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
  print(g)
  cat("Plotting taxa: ", taxon, "\n")
}
dev.off()

# Joint model

jsdm$deviance <- arrange(jsdm$deviance, desc(n))
d <- jsdm$deviance$D2; names(d) <- jsdm$deviance$Taxon

x <- prepare.inputs(c("Temp", "FV", "F10m", "IAR", "Urban", "LUD"), sample.bdms, center=FALSE)
x <- gather(x, Variable, Value, -SiteId, -SampId)
setDT(x)

plot.data <- left_join(jsdm$probability, x, by=c("SiteId", "SampId"))
plot.data <- filter(plot.data, Variable != "Temp2")

plot.data$Obs <- as.factor(plot.data$Obs)
plot.data$Label <- factor(plot.data$Variable, levels = c("Temp", "FV", "F10m", "IAR", "Urban", "LUD"))
levels(plot.data$Label) <- c("Temp" = expression(paste("Temp (", degree, "C)")),
                             "FV" = expression(paste("FV (m/s)")),
                             "F10m" = expression(paste("F10m (%)")),
                             "IAR" = expression(paste("IAR (spray treatments * fraction cropland)")),
                             "Urban" = expression(paste("Urban (%)")),
                             "LUD" = expression(paste("LUD (CE/km"^2,")")))
setDT(plot.data)

pdf(paste('jSDM prob vs X.pdf', sep=''), paper = 'special', height = 9, width = 12, onefile = TRUE)
for (j in 1:length(names(d))){
  taxon <- names(d)[j]
  d.squared <- round(d[j], 2)
  
  taxon.prob <- plot.data[Taxon==taxon,]
  
  g <- ggplot(data = taxon.prob, aes(x = Value, y = Pred, color = Obs))
  g <- g + geom_point(alpha = 0.25)
  g <- g + theme_bw(base_size=15)
  g <- g + facet_wrap(~ Label, scales = "free_x", labeller=label_parsed, strip.position="bottom")
  g <- g + labs(title = "", # paste(j, " - ", taxon)
                x = "",
                y = "Probability of occurrence")
  g <- g + theme(axis.title.y = element_text(size = 20), strip.text.x = element_text(size = 12), strip.background = element_blank(), strip.placement = "outside")
  # g <- g + scale_colour_manual(values=c("1" = "blue", "0" = "red"))
  g <- g + scale_color_manual(name = "Observation", values=c("#FF0000", "#0077FF"), labels=c("Absence", "Presence"))
  g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
  print(g)
  cat("Plotting taxa: ", taxon, "\n")
}
dev.off()

# Additional analyses (depreciated) ####
# ### Quantile regression ###
# # Regression of quantiles of D2 against occurrence frequency (D2 ~ n)
# iqr <- rq(D2 ~ n, data = deviance[Model=='iSDM',], tau = c(1:9/10))
# jqr <- rq(D2 ~ n, data = deviance[Model=='jSDM',], tau = c(1:9/10))
# 
# ### quantile regression with ggplot2
# g <- ggplot(data = deviance, aes(x = n, y = D2))
# g <- g + geom_point(alpha = 3/10, aes(x = n, y = D2, size = n))
# g <- g + geom_quantile(quantiles = c(1:9/10))
# g <- g + facet_grid(~ Model)
# # g <- g + scale_y_continuous(limits = c(0, 1))
# g <- g + theme_bw(base_size = 18)
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g <- g + labs(title = expression(paste("Calibration D"^2, " vs occurrence frequency by model")),
#               subtitle = expression(paste("Quantile regression of D"^2, " ~ occurrence frequency (10% quantile intervals)")),
#               x = "Total occurrences",
#               y = expression(paste("D"^2)))
# g <- g + labs(size = "Total occurrences")
# print(g)
# 
# # Residual deviance: compare models ###
# deviance.plot <- select(d.squared.all, Model, Species, res.dev)
# deviance.plot <- spread(deviance.plot, Model, res.dev)
# deviance.plot$n <- n[deviance.plot$Species]

# g <- ggplot(data = deviance.plot, aes(x = jSDM, y = iSDM, size = n))
# g <- g + geom_point(alpha = 0.3)
# g <- g + labs(title = "iSDM vs jSDM deviance (proposed model, entire BDM)")
# g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
# print(g)

# Plot D2 by model ###
# d.squared.all.plot <- spread(d.squared.all, Model, D2)
# 
# # d.squared.all.plot$Shape <- ifelse(d.squared.all.plot$iSDM < 0 | d.squared.all.plot$jSDM < 0, 17, 16)
# # d.squared.all.plot$iSDM[d.squared.all.plot$iSDM < 0] <- 0
# # d.squared.all.plot$jSDM[d.squared.all.plot$jSDM < 0] <- 0
# 
# g <- ggplot(d.squared.all.plot, aes(x = iSDM, y = jSDM))
# g <- g + geom_point(aes(size = n), alpha = 0.4, color = 'black')
# g <- g + scale_size_continuous(range = c(2, 10))
# g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
# g <- g + theme_bw(base_size = 18)
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# 
# g <- g + labs(title = expression(paste("Calibration D"^2," of joint vs individual model")),
#               subtitle = "Goodness-of-fit for all BDM sites",
#               x = expression(paste("iSDM D"^2)),
#               y = expression(paste("jSDM D"^2)))
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g <- g + guides(colour = guide_legend(override.aes = list(size=6)), shape = FALSE)
# g <- g + labs(size = "Total occurrences")
# print(g)

### Plot D2 ratio ###
# fit.plot <- d.squared.all.plot
# fit.plot$ratio <- fit.plot$jSDM/fit.plot$iSDM
# plot(fit.plot$n, fit.plot$ratio)

# g <- ggplot(data = fit.plot, aes(x = n, y = ratio, size = n))
# g <- g + geom_point(alpha = 0.3)
# g <- g + geom_vline(xintercept=50, color = "black", linetype = "longdash")
# g <- g + theme_bw(base_size = 16)
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g <- g + labs(title = expression(paste("D"^2, " ratio for jSDM vs iSDM")),
#               x = "Total occurrences",
#               # y = expression(paste("D"["jSDM"]^2,"/", "D"["iSDM"]^2)))
#               y = expression(paste(frac("D"["jSDM"]^2,"D"["iSDM"]^2))))
# # y = "")
# g <- g + labs(size = "Total occurrences")
# print(g)

### Plot hist(dev ratios) ###
# dt.hist <- d
# dt.hist$ratio <- dt.hist$D2.j/dt.hist$D2.i

# ggplot
# bw <- 0.5
# 
# g <- ggplot(data = dt.hist, aes(x = ratio))
# g <- g + geom_histogram(binwidth = bw, fill = "white", color = "black")
# g <- g + stat_bin(binwidth= bw, geom="text", aes(label=..count..), vjust=-1.5)
# g <- g + geom_vline(xintercept = 1, color = "black", linetype = "longdash")
# g <- g + theme_bw(base_size = 14)
# g <- g + theme(plot.title=element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g <- g + labs(title = expression(paste("Frequency of D"^2," ratios between jSDM and iSDM")),
#               # subtitle = expression(paste("D"^2, " ratio = D"["jSDM"]^2, "/ ", "D"["iSDM"]^2)),
#               x = expression(paste("D"^2," ratio", " (D"["jSDM"]^2, " / ", "D"["iSDM"]^2,")")),
#               y = "Frequency")
# g <- g + scale_x_continuous(breaks=c(-3,-2,-1,0,1,2,3))
# print(g)

# pdf('beta_comm_sensitivetaxa_plot.pdf', onefile = TRUE)
# par(cex=1.25)
# for (k in 1:length(inf.fact)){
#   variable <- inf.fact[k]
#   cat("Plotting: ",variable,"\n")
#   # Test temperature as starting variable
#   samples <- taxon.samples[Variable == variable,]
#   
#   # Find the maximum density among the posterior taxon-specific posterior parameters
#   taxon.specific.sd <- samples %>%
#     group_by(Taxon) %>%
#     summarise(SD = sd(Value)) %>%
#     arrange(SD) # arrange from min to max
#   
#   sd <- taxon.specific.sd[1,]$Taxon
#   ymax.sample <- samples[Taxon == sd, ]
#   ymax <- density(ymax.sample$Value)
#   ymax <- max(ymax$y)
#   
#   # Plot the community parameter distribution
#   mu <- mu.beta.comm.maxpost[variable]
#   sd <- sigma.beta.comm.maxpost[variable]
#   x <- seq(mu-4*sd,mu+4*sd,length.out=201)
#   beta.comm.density <- dnorm(x, mu, sd)
#   
#   # Match expressions to influence factors
#   
#   labels <- c("Temp" = expression(paste(beta["Temp"], " (1/", degree, "C)")),
#               "Temp2" = expression(paste(beta["Temp"^2], " (1/", degree, "C)")),
#               "FV" = expression(paste(beta["FV"], " (s/m)")),
#               "F10m" = expression(paste(beta["F10m"], " (1/%)")),
#               "IAR" = expression(paste(beta["IAR"], " (spray treatments * fraction cropland)")),
#               "Urban" = expression(paste(beta["Urban"], " (1/%)")),
#               "LUD" = expression(paste(beta["LUD"], " (km"^2,"/CE)"))
#   )
#   x.label <- labels[variable]
#   
#   plot(numeric(0), numeric(0),
#        xlim = c(min(x), max(x)),
#        ylim=c(0, ymax),
#        xlab = x.label,
#        ylab = "Density")
#   
#   # Plot the taxon-specific parameters
#   for (j in 1:length(names(n))){
#     taxon <- names(n)[j]
#     sample <- samples[Taxon==taxon,]
#     sample <- sample$Value
#     # Fill matrix: beta.taxa.response
#     jsig <- quantile(sample, probs = c(0.05, 0.95))
#     
#     # If 5th quantile greater than 0, set positive
#     # If 95th quantile less than 0, set negative
#     if (jsig[1] > 0){ # significant positive
#       c <- "blue"
#     }
#     if (jsig[2] < 0){ # significant negative
#       c <- "red"
#     }
#     # If posterior is !(positive) AND !(negative), set grey
#     if (!(jsig[1] > 0) & !(jsig[2] < 0)){
#       c <- "grey55"
#     }
#     
#     # beta.taxa.response[j, k] <- c
#     if (taxon %in% btr.int$Taxon){
#       lines(density(sample), type="l", col=alpha(c, 0.30), lwd=1)
#     }
#   }
#   # Plot community parameter distribution
#   lines(x, beta.comm.density, type="l", col = "grey50", lwd = 5)
#   # Plot maximum posterior values over all taxa
#   lines(density(beta.taxa.maxpost[variable, ]), type="l", col="black", lwd=2, lty='longdash')
# }
# dev.off()



# # Gather quantile values - requires normalized parameters
# beta.taxa.quantiles <- matrix(nrow=length(n.bdms), ncol=length(inf.fact))
# for (k in 1:length(inf.fact)){
#   variable <- inf.fact[k]
#   # Test temperature as starting variable
#   samples <- taxon.samples[Variable == variable,]
#   
#   # Plot the taxon-specific parameters
#   for (j in 1:length(names(n))){
#     taxon <- names(n)[j]
#     sample <- samples[Taxon==taxon,]
#     sample <- sample$Value
#     
#     jsig <- quantile(sample, probs = c(0.05, 0.95))
#     if (jsig[1] > 0){ # If 5th quantile greater than 0, set positive
#       response <- jsig[1]
#     }
#     if (jsig[2] < 0){ # If 95th quantile less than 0, set negative
#       response <- jsig[2]
#     }
#     # If posterior is !(positive) AND !(negative), set grey
#     if (!(jsig[1] > 0) & !(jsig[2] < 0)){
#       response <- 0
#     }
#     beta.taxa.quantiles[j, k] <- response
#   }
# }
# 
# beta.taxa.quantiles <- as.data.table(beta.taxa.quantiles)
# beta.taxa.quantiles$Taxon <- names(n)
# 
# btq.int <- filter(beta.taxa.quantiles, Taxon %in% beta.int$Taxon)
# colnames(btq.int) <- c(inf.fact, "Taxon")
# btq.plot <- btq.int[, c("Taxon", inf.fact)]
# 
# row.names(btq.plot) <- paste(sub("_", " ", btq.plot$Taxon), " - ", n[btq.plot$Taxon])
# btq.plot$Taxon <- NULL
# 
# pheatmap(btq.plot, cellwidth=30, cellheight=11, cluster_rows=T, kmeans_k=5, cluster_cols=F, legend=T, show_rownames=TRUE, fontsize=13, fontsize_col=15)


# # Deviance w/o RE ###
# # Combine deviances for jSDM w/o random effect
# bdm.ch.re$deviance$RE <- "RE"
# bdm.ch$deviance$RE <- "noRE"
# 
# jsdm.deviance <- bind_rows(bdm.ch$deviance, bdm.ch.re$deviance)
# 
# test <- jsdm.deviance %>%
#   select(Species, std.deviance, RE) %>%
#   spread(RE, std.deviance) %>%
#   mutate(n = n[test$Species])
# 
# plot(test$noRE, test$RE)
# g <- ggplot(test, aes(x = noRE, y = RE, size = n))
# g <- g + geom_point(alpha = 0.3)
# g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
# 
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# # Adjust legends
# g <- g + guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))
# 
# g <- g + theme_bw(base_size = 15)
# g <- g + labs(x = "Without random effect",
#               y = "With random effect",
#               title = "Standardized deviance of joint model taxa")
# print(g)

# bdm.glm.prob <- ExtractFV(bdm.glm.ch)
# # Second dataset needed: predictors (raw data)
# pg <- 1
# 
# pdf(paste('outputs/iSDM_prob.pdf', sep=''), paper = 'special', height = 9, width = 12, onefile = TRUE)
# 
# d <- sapply(bdm.glm.ch$models, function(m){
#   (m$null.deviance-m$deviance)/m$null.deviance
# })
# 
# for (t in 1:length(names(bdm.glm.ch$taxa))){
#   taxon <- names(bdm.glm.ch$taxa)[t]
#   d.squared <- round(d[t], 2)
#   
#   taxon.prob <- bdm.glm.prob[bdm.glm.prob$Species == taxon, ]
#   influence.factors <- gather(x, Variable, Value, -SiteId, -SampId)
#   plot.data <- left_join(taxon.prob, influence.factors, by = c("SiteId", "SampId"))
#   plot.data <- filter(plot.data, Variable != "Temp2")
#   plot.data$Obs <- as.factor(plot.data$Obs)
#   plot.data$Variable <- factor(plot.data$Variable, levels = c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD"))
#   levels(plot.data$Variable) <- c(expression(paste("Temp (", degree, "C)")),
#                                   expression(paste("Temp"^2," (", degree, "C)")),
#                                   expression(paste("FV (m/s)")),
#                                   expression(paste("F10m (%)")),
#                                   expression(paste("IAR")),
#                                   expression(paste("Urban (%)")),
#                                   expression(paste("LUD (CE/km"^2,")")))
#   
#   g <- ggplot(data = plot.data, aes(x = Value, y = Pred, color = Obs))
#   g <- g + geom_point(alpha = 0.25)
#   g <- g + facet_wrap(~ Variable, scales = "free", labeller=label_parsed)
#   
#   g <- g + labs(title = paste("Predicted probabilities vs influence factors of",paste(taxon)),
#                 subtitle = paste("iSDM: y ~", paste(bdm.glm.ch$inf.fact, collapse = " + ", sep = " "), "- page", pg, "\n D2 = ", d.squared),
#                 x = "",
#                 y = "Probability",
#                 color = "Observation")
#   g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   g <- g + scale_colour_manual(values=c("1" = "blue", "0" = "red"))
#   g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
#   g <- g + theme_bw(base_size = 15)
#   cat("Plotting taxa: ", taxon, "\n")
#   pg <- pg + 1
#   print(g)
# }
# dev.off()

# Plots: basic deviance ###
# pdf('dev vs n.pdf')
# par(mfrow=c(1,2))
# # individual model
# plot(bdm.glm.deviance$n, bdm.glm.deviance$null.dev, col="black", cex=1.5,
#      main="Individual models",
#      ylab="Null (black) and proposed (red) model deviance",
#      xlab="Occurrence frequency")
# points(bdm.glm.deviance$n, bdm.glm.deviance$res.dev, col="red", cex=1.5)
# 
# # joint model
# plot(bdm.ch$deviance$n, bdm.ch$deviance$null.dev, col="black", cex=1.5,
#      main="Joint model",
#      ylab="Null deviance (black), proposed model deviance (red)",
#      xlab="Occurrence frequency")
# points(bdm.ch$deviance$n, bdm.ch$deviance$res.dev, col="red", cex=1.5)
# dev.off()

# pdf('outputs/deviance_plots.pdf', paper='special', width = 14, height = 9, onefile = TRUE)
# par(mfrow=c(2,3), oma=c(0,0,0,0),  mgp=c(1.5, 0.5, 0))
# # Calibrated iSDM: null and residual deviance
# plot(bdm.glm.deviance$n, bdm.glm.deviance$null.dev, col="black", cex=1.5,
#      main="Calibration iSDM: null dev vs n",
#      ylab="Null deviance (black), proposed model deviance (red)",
#      xlab="Occurrence frequency")
# points(bdm.glm.deviance$n, bdm.glm.deviance$res.dev, col="red", cex=1.5)
# 
# # Calibrated iSDM: null - residual deviance
# plot(bdm.glm.deviance$n, bdm.glm.deviance$null.dev-bdm.glm.deviance$res.dev, col="black", cex=1.5,
#      main="Calibration iSDM: null-res vs n",
#      ylab="Null minus proposed model deviance",
#      xlab="Occurrence frequency")
# 
# # Calibrated iSDM: D2
# plot(bdm.glm.deviance$n, bdm.glm.deviance$D2, col="black", cex=1.5, ylim=c(-0.2,1),
#      main="Calibration iSDM: D2 vs n",
#      ylab="D2",
#      xlab="Occurrence frequency")
# 
# # Prediction iSDM: null and residual deviance
# data <- filter(cv.glm.deviance, Type=='Testing')
# plot(data$n, data$null.dev, col="black", cex=1.5,
#      main="Prediction iSDM: null dev vs n",
#      ylab="Null deviance (black), proposed model deviance (red)",
#      xlab="Occurrence frequency")
# points(data$n, data$residual.deviance, col="red", cex=1.5)
# 
# # Prediction iSDM: null - residual deviance
# plot(data$n, data$null.dev-data$residual.deviance, col="black", cex=1.5,
#      main="Prediction iSDM: null-res dev vs n",
#      ylab="Null deviance - proposed model deviance",
#      xlab="Occurrence frequency")
# 
# # Prediction iSDM: D2
# plot(data$n, data$D2, col="black", cex=1.5, ylim=c(-10,1),
#      main="Prediction iSDM: D2 vs n",
#      ylab="D2",
#      xlab="Occurrence frequency")
# 
# 
# 
# # Calibration jSDM: null deviance
# data <- bdm.ch$deviance
# plot(data$n, data$null.dev, col="black", cex=1.5,
#      main="Calibration jSDM: null dev vs n",
#      ylab="Null deviance (black), proposed model deviance (red)",
#      xlab="Occurrence frequency")
# points(data$n, data$res.dev, col="red", cex=1.5)
# 
# # Calibration jSDM: null - residual deviance
# plot(data$n, data$null.dev-data$res.dev, col="black", cex=1.5, ylim=c(0,500),
#      main="Calibration jSDM: null-res dev vs n",
#      ylab="Null deviance - proposed model deviance",
#      xlab="Occurrence frequency")
# 
# # Calibration jSDM: D2
# plot(data$n, data$D2, col="black", cex=1.5, ylim=c(-0.2,1),
#      main="Calibration jSDM: D2 vs n",
#      ylab="D2",
#      xlab="Occurrence frequency")
# 
# # Prediction jSDM: null deviance
# data <- filter(cv.jsdm$deviance, Type == 'Testing')
# plot(data$n, data$null.dev, col="black", cex=1.5,
#      main="Prediction jSDM: null dev vs n",
#      ylab="Null deviance (black), proposed model deviance (red)",
#      xlab="Occurrence frequency")
# points(data$n, data$residual.dev, col="red", cex=1.5)
# 
# # Prediction jSDM: null - residual deviance
# plot(data$n, data$null.dev-data$residual.deviance, col="black", cex=1.5,
#      main="Prediction jSDM: null-res vs occurrence frequency",
#      ylab="Null deviance minus proposed model deviance",
#      xlab="Occurrence frequency")
# 
# # Prediction jSDM: D2
# plot(data$n, data$D2, col="black", cex=1.5, ylim=c(-10,1),
#      main="Prediction jSDM: D2 vs occurrence frequency",
#      ylab="D2",
#      xlab="Occurrence frequency")
# dev.off()

# # Combine parameters #
# c <- intersect(colnames(bdm.glm.ch$parameters), colnames(bdm.ch$parameters))
# 
# beta.all <- bind_rows(bdm.glm.beta[, c], bdm.ch$parameters[, c])
# 
# # Plot parameters: omit outliers #
# beta.all.plot <- beta.all
# 
# temp <- aggregate(Parameter ~ Variable, beta.all.plot[beta.all.plot$Model == 'jSDM'], IQR)
# viqr <- temp$Parameter
# names(viqr) <- temp$Variable
# 
# temp <- aggregate(Parameter ~ Variable, beta.all.plot[beta.all.plot$Model == 'jSDM'], function(x){
#   quantile(x, 0.25)
# })
# q1 <- temp$Parameter
# names(q1) <- temp$Variable
# temp <- aggregate(Parameter ~ Variable, beta.all.plot[beta.all.plot$Model == 'jSDM'], function(x){
#   quantile(x, 0.75)
# })
# q3 <- temp$Parameter
# names(q3) <- temp$Variable
# rm(temp)
# 
# # Set parameters values greater or less than 2.5*sd to NA (i.e., outliers)
# beta.all.plot$Parameter[beta.all.plot$Parameter > q3[beta.all.plot$Variable] + (viqr[beta.all.plot$Variable]*5)] <- NA # maximum
# beta.all.plot$Parameter[beta.all.plot$Parameter < q1[beta.all.plot$Variable] - (viqr[beta.all.plot$Variable]*5)] <- NA # minimum
# 
# g <- ggplot(beta.all.plot, aes(x = Parameter, fill = Model))
# g <- g + geom_density(alpha = 0.3, na.rm = TRUE)
# g <- g + facet_wrap(~ Variable, scales = "free")
# g <- g + theme_bw(base_size = 14)
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g <- g + labs(title = "Distributions of parameters for joint and individual SDM by influence factor",
#               subtitle = "Maximum posterior and maximum likelihood estimates for BDM dataset",
#               x = expression(paste(beta, " parameter value")),
#               y = "Density")
# g <- g + scale_fill_manual(values= c("jSDM" = "red", "iSDM" = "blue"))
# g <- g + labs(fill = "Model")
# print(g)
# 
# # Plot parameters: keep outliers #
# colnames(bdm.glm.beta)[1] <- "Species"
# bdm.glm.beta$Model <- "iSDM"
# beta.all <- rbind(bdm.glm.beta, bdm.ch$parameters) # MLE and MP estimates
# beta.plot <- beta.all
# 
# # pdf('outputs/beta_plot.pdf', paper='special',width=16, height=12, onefile = TRUE)
# # par(mfrow = c(2, 4), mar=c(4,3,5,0.5), oma = c(2,2,2,2), cex.main = 1.25, cex.lab = 1.25, cex.sub = 1.25)
# pdf('outputs/ijSDM beta parameters.pdf', paper='special', width = 9, height = 11, onefile = TRUE)
# par(mfrow = c(4,2),  mar = c(3, 1.5, 2.5, 1.5), oma = c(2, 3, 2.5, 2))
# letters <- c("a", "b", "c", "d", "e", "f", "g")
# for (k in 1:uniqueN(beta.plot$Variable)){
#   # subset data for analysis
#   variable <- unique(beta.plot$Variable)[k]
#   letter <- letters[k]
#   # parameters by variable
#   plot.data <- filter(beta.plot, Variable == variable)
#   
#   # parameters by variable and model
#   isdm.p <- filter(plot.data, Model == "iSDM")
#   isdm.p <- isdm.p$Parameter
#   jsdm.p <- filter(plot.data, Model == "jSDM")
#   jsdm.p <- jsdm.p$Parameter
#   
#   # calculate densities with all data
#   
#   jsdm.density <- density(jsdm.p)
#   isdm.density <- density(isdm.p)
#   isdm.density <- density(isdm.p, bw = isdm.density$bw*9)
# 
#   # get plot limits
#   xmin <- 1.25*min(jsdm.p)
#   xmax <- 1.25*max(jsdm.p)
#   
#   ymin <- 0
#   ymax.isdm <- max(isdm.density$y)
#   ymax.jsdm <- max(jsdm.density$y)
#   ymax <- ifelse(ymax.isdm > ymax.jsdm, ymax.isdm, ymax.jsdm)
#   
#   # setup empty plot based on calculated XY limits
#   plot.title <- c(expression(paste("Temp (", degree, "C)")),
#                   expression(paste("Temp"^2," (", degree, "C"^2,")")),
#                   expression(paste("FV (m/s)")),
#                   expression(paste("F10m (%)")),
#                   expression(paste("IAR")),
#                   expression(paste("Urban (%)")),
#                   expression(paste("LUD (CE/km"^2,")")))
#   
#   # expression(paste("Standard deviation of ", beta["jk"]^"taxa"))
#   plot(numeric(0), numeric(0), xlim=c(xmin,xmax), ylim=c(ymin,ymax),
#        main = plot.title[k],
#        # main = expression(paste(beta[bquote(.(variable))]^"taxa")),
#        ylab = "",
#        xlab = "", 
#        cex.main=1.5, cex.sub=1.5)
#   mtext(letter, side = 3, line = -3, adj = 0.1, cex = 1.5)
#   # Plot densities
#   lines(jsdm.density, type="l", col="red", lwd = 3) # jsdm; e.g. transparency: alpha("red", 0.3)
#   lines(isdm.density, type="l", col="blue", lwd = 3) # isdm
#   # Plot rug lines
#   rug(jsdm.p, ticksize = 0.03, side = 3, lwd = 0.5, col = "red")
#   rug(isdm.p, ticksize = 0.03, side = 1, lwd = 0.5, col = "blue")
# }
# # Labels
# mtext(text="Maximum posterior and maximum likelihood parameter distributions for joint vs individual model",side=3,line=1,outer=TRUE)
# mtext(text="Parameter value",side=1,line=1,outer=TRUE)
# mtext(text="Density",side=2,line=1.5,outer=TRUE)
# par(fig = c(0, 1, 0, 1), mar = c(2.5, 0, 0, 2), new = TRUE)
# plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
# legend("bottomright",legend = c("iSDM", "jSDM"),
#        col=c("blue", "red"), lwd=3, cex=1.5, horiz = FALSE)
# dev.off()
# 
# plot.data <- beta.plot[beta.plot$Variable == "Temp",]
# g <- ggplot(data = plot.data, aes(x = Parameter, fill = Model))
# g <- g + geom_density(alpha = 0.3, position="identity")
# g <- g + coord_cartesian(xlim = c(xmin, xmax), ylim = c(ymin, ymax))
# g <- g + theme_bw(base_size = 14)
# g <- g + scale_fill_manual(values= c("jSDM" = "red", "iSDM" = "blue"))
# g <- g + labs(fill = "Model")
# print(g)
# plots[[k]] <- g

# # Combine CV deviances ###
# cv <- bind_rows(cv.jsdm$deviance, cv.glm.deviance)
# 
# # Get the mean deviance over all folds and arrange the data for plotting
# test <- cv %>%
#   select(Species, Type, Fold, Model, D2) %>%
#   filter(Type == 'Testing') %>%
#   group_by(Species, Model) %>%
#   summarise(mean.D2 = mean(D2)) %>%
#   spread(Model, mean.D2)
# 
# test <- cv %>%
#   select(Species, Type, Fold, Model, residual.deviance) %>%
#   filter(Type == 'Testing') %>%
#   group_by(Species, Model) %>%
#   summarise(mean.dev = mean(residual.deviance)) %>%
#   spread(Model, mean.dev)
# 
# test$n <- n[test$Species]
# 
# test$Shape <- ifelse(is.infinite(test$iSDM), 24, 16)
# # cv.plot$n <- ifelse(is.infinite(cv.plot$Testing), NA, cv.plot$n)
# test$Shape <- as.factor(test$Shape)
# test$iSDM <- ifelse(is.infinite(test$iSDM), max(test$iSDM[!is.infinite(test$iSDM)]), test$iSDM)
# 
# g <- ggplot(data = test, aes(x = jSDM, y = iSDM, size = n, shape = Shape))
# g <- g + geom_point(alpha = 0.3)
# g <- g + labs(title = "iSDM vs jSDM mean predictive deviance",
#               x = "Joint model",
#               y = "Individual model",
#               size = "Occurrence frequency",
#               shape = "Deviance")
# g <- g + theme_bw(base_size = 13)
# g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
# g <- g + scale_y_continuous(limits=c(0,NA))
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# # Adjust legends
# g <- g + guides(colour = guide_legend(override.aes = list(size=6)), shape = guide_legend(override.aes = list(size=6)))
# g <- g + labs(colour = "Model",shape = "Deviance", size = "Total occurrences")
# g <- g + scale_shape_discrete(name  = "Deviance",
#                               breaks=c("16", "24"),
#                               labels=c("Finite", "Infinite"))
# print(g)
# 
# # Plot test/train dev ratios ###
# dev.ratio.plot <- cv.deviance.mean.all
# dev.ratio.plot$Ratio <- dev.ratio.plot$Testing/dev.ratio.plot$Training
# # dev.ratio.plot$Ratio <- log10(dev.ratio.plot$Testing/dev.ratio.plot$Training)
# 
# # compute lower and upper whiskers
# ylim1 = boxplot.stats(dev.ratio.plot$Ratio)$stats[c(1, 5)]
# 
# g <- ggplot(dev.ratio.plot, aes(Model, Ratio))
# g <- g + geom_boxplot(alpha = 0.3, width = 0.6)
# g <- g + geom_hline(yintercept = 1, color = "black", linetype = "longdash")
# g <- g + theme_bw(base_size = 18)
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# g <- g + labs(title = "Prediction and calibration deviance by model",
#               x = "Model",
#               y = expression(paste("mean relative deviance (prediction / calibration)")))
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g.final <- g + coord_flip(ylim = ylim1*1.05)
# # g.final <- g.final + coord_cartesian(ylim = ylim1*1.05) # scale y limits based on ylim1
# print(g.final)

# Plot overall deviance ###
# model.dev.plot <- model.dev[!is.infinite(model.dev$D2) & model.dev$D2 > 0,] # remove Inf and D2 < 0.
# model.dev.plot <- model.dev[!is.infinite(model.dev$D2), ]
# dplot <- cv.dev
# dplot$Type <- factor(dplot$Type, levels = c("Training", "Testing"))
# 
# g <- ggplot(dplot, aes(x = Type, y = D2))
# g <- g + geom_boxplot(aes(fill=Type))
# g <- g + facet_grid(. ~ Trial + Fold + Type)
# g <- g + labs(x="")
# g <- g + theme_bw()
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# g <- g + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# g <- g + theme(plot.title = element_text(hjust = 0.5))
# g <- g + labs(title = "Training/testing relative deviance for iSDM vs jSDM in overall community (50-50)
#                 \n relative deviance = (resid.dev^2) / (number of training sites)",
#               x = "",
#               y = "Relative deviance")
# g <- g + scale_fill_brewer(palette = "Set1")
# print(g)

# g <- ggplot(d, aes(x = Type, y = D2))
# g <- g + geom_point(aes(fill=Trial))
# g <- g + facet_grid(. ~ Fold + Type)
# g <- g + labs(x="")
# g <- g + theme_bw()
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# g <- g + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# g <- g + theme(plot.title = element_text(hjust = 0.5))
# g <- g + scale_fill_brewer(palette = "Set1")
# print(g)



### How much better is jSDM at prediction? ###
# test <- cv.deviance.mean.all[, c("Species", "Model", "Testing")]
# test.plot <- spread(test, Model, Testing)
# test.plot$n <- n[test.plot$Species]
# 
# g <- ggplot(data = test.plot, aes(x = jSDM, y = iSDM, size = n))
# g <- g + geom_point(alpha = 0.3)
# g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
# print(g)
# 
# jadv <- spread(cv.deviance.mean.all[, c("Species", "Model", "Testing")], Model, Testing)
# jadv$n <- n[jadv$Species]
# jadv$ratio <- 1-(jadv$jSDM/jadv$iSDM)
# 
# hist(jadv$ratio, col = "grey",
#      main = "Predictive deviance ratio for jSDM vs iSDM",
#      xlab = "Predictive deviance ratio (1 - jSDM / iSDM)"
# )
# 
# jadv.plot <- jadv[jadv$ratio > -1,] # drop (two) species >100% worse in jSDM predictive performance
# g <- ggplot(data = jadv.plot, aes(x = n, y = ratio, size = n))
# g <- g + geom_point(alpha = 0.3)
# g <- g + theme_bw(base_size = 16)
# 
# g <- g + theme(strip.background=element_rect(fill="black"))
# g <- g + theme(strip.text=element_text(color="white", face="bold"))
# g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
# g <- g + labs(title = "Predictive deviance ratio for jSDM vs iSDM",
#               x = "Total occurrences",
#               y = "Predictive deviance ratio (1 - jSDM / iSDM)")
# g <- g + labs(size = "Total occurrences")
# print(g)
# 
# test <- jadv.plot
# hist(test$n[test$ratio>0], col = "grey", labels=TRUE,
#      main = "Occurrence frequency of taxa with predictive deviance ratio > 0",
#      xlab = "Occurrence frequency")