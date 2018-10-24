# ANALYSES
# 
# - run three JSDMs:
# (1) joint model WITHOUT correlated community parameters (in progress)
# (2) joint model with correlated community parameters, WITHOUT site effects or latent variable (in progress)
# (3) joint model with correlated community parameters, WITH site effects and latent variable (in progress)
# 
# - compare the quality of fit and predictive performance of (1) and (2)
# - compare how the community parameters change in (1) and (2)
# - compare how the taxon-specific parameters change between (1) and (2)
# 
# - map the latent variable (x_lat) and the marginal posterior factor loading for each taxon (beta_lat) (DONE)
# - calculate linear predictor of the latent variable
# - plot maximum posterior correlations of community parameters (DONE)

# Paper 1 extended results
# Based on results of untransformed variables, prepare variables for model run
# deploy.jsdm() ####
K <- c("Temp", "Temp2", "FV", "F100m", "LUD", "IAR")
predictors <- prepare.inputs(K, sample.bdms, center = TRUE)
write.csv(predictors, 'outputs/predictors.csv', row.names = FALSE)

deploy.jsdm(K, sample.bdms, "FF0", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TF0", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT0", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT1", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT2", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT3", center=T, cv=0)
deploy.jsdm(K, sample.bdms, "TT4", center=T, cv=0)

deploy.jsdm(K, sample.bdms, "FF0_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "FF0_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "FF0_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TF0_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TF0_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TF0_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT0_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT0_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT0_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT1_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT1_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT1_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT2_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT2_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT2_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT3_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT3_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT3_train3", center=T, cv=3)

deploy.jsdm(K, sample.bdms, "TT4_train1", center=T, cv=1)
deploy.jsdm(K, sample.bdms, "TT4_train2", center=T, cv=2)
deploy.jsdm(K, sample.bdms, "TT4_train3", center=T, cv=3)

# Extract JSDMs ####
jsdm.FF0 <- extract.jsdm('outputs/paper 1 extensions', 'FF0') # simplest joint model
jsdm.TF0 <- extract.jsdm('outputs/paper 1 extensions', 'TF0') # correlated beta.comm
jsdm.TT0 <- extract.jsdm('outputs/paper 1 extensions', 'TT0') # correlated beta.comm, site effects
jsdm.TT1 <- extract.jsdm('outputs/paper 1 extensions', 'TT1')
jsdm.TT2 <- extract.jsdm('outputs/paper 1 extensions', 'TT2')

# CAUTION: chain merging is not yet implemented for more than two latent variables
jsdm.TT3 <- extract.jsdm('outputs/paper 1 extensions', 'TT3')
jsdm.TT4 <- extract.jsdm('outputs/paper 1 extensions', 'TT4')
# jsdm.p1 <- extract.jsdm.extensions('outputs/paper 1 results', 'jsdm_bdms')

# Quality of fit ####
deviance.jsdm <- bind_rows(jsdm.FF0$deviance, jsdm.TF0$deviance, jsdm.TT0$deviance, jsdm.TT1$deviance, jsdm.TT2$deviance)
deviance.isdm <- isdm$deviance
deviance.isdm$Trial <- "iSDM"
deviance.isdm <- select(deviance.isdm, Taxon, Model, null.dev, res.dev, std.res.dev, D2, n.samples, n, Trial)

deviance <- bind_rows(deviance.jsdm, deviance.isdm)

plot.data <- deviance
plot.data$Prevalence <- plot.data$n / plot.data$n.samples
plot.data$TrialLabel <- factor(plot.data$Trial, levels=c("Null model", "iSDM", "FF0", "TF0", "TT0", "TT1", "TT2"))
plot.data$TaxonLabel <- factor(plot.data$Taxon, levels=c(names(sort(n.bdms, decreasing=FALSE))))

null.data <- plot.data %>%
  as.tibble() %>%
  filter(Model=="iSDM") %>%
  mutate(std.res.dev = null.dev/n.samples,
         Trial = "Null model",
         TrialLabel = factor(Trial, levels=c("Null model", "iSDM", "FF0", "TF0", "TT0", "TT1", "TT2")))

plot.data <- bind_rows(plot.data, null.data)
rm(null.data)

g1 <- ggplot(data = plot.data)
g1 <- g1 + geom_point(aes(x=TaxonLabel, y=std.res.dev, color=TrialLabel, size=Prevalence), alpha=0.35)
g1 <- g1 + stat_smooth(aes(x=TaxonLabel, y=std.res.dev, group=TrialLabel, color=TrialLabel), geom='line', size = 2, alpha=0.50, se=FALSE, method = "loess")
# g1 <- g1 + stat_smooth(data = null.data, aes(x=TaxonLabel, y=std.res.dev, group=TrialLabel), color="black", linetype = "dashed", geom='line', size = 2, alpha=0.50, se=FALSE, method = "loess")
g1 <- g1 + labs(title = "Quality of fit for entire community by model",
             y = "Standardized deviance",
             x = "Taxa (increasing prevalence from left to right)",
             color = "Model",
             subtitle = "Curves fitted with LOESS method")
g1 <- g1 + theme(axis.text.x = element_blank(),
               axis.ticks.x = element_blank())
g1 <- g1 + scale_color_brewer(palette = "Set1")

g2 <- ggplot(data = plot.data)
g2 <- g2 + geom_boxplot(aes(x=TrialLabel, y=std.res.dev, fill=TrialLabel))
g2 <- g2 + scale_fill_brewer(palette = "Set1")

pdf('outputs/paper 1 extensions/P1 - quality of fit.pdf', width = 20, paper='a4r')
ggarrange(g1, g2, ncol=2, common.legend = FALSE, legend="right", align="hv")
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
# jsdm.FF0$tjur <- extract.tjur(jsdm.FF0$probability)
# jsdm.TF0$tjur <- extract.tjur(jsdm.TF0$probability)
# jsdm.TT0$tjur <- extract.tjur(jsdm.TT0$probability)
# jsdm.TT1$tjur <- extract.tjur(jsdm.TT1$probability)
# jsdm.TT2$tjur <- extract.tjur(jsdm.TT2$probability)
# 
# jsdm.FF0$tjur$Model <- "FF0"
# jsdm.TF0$tjur$Model <- "TF0"
# jsdm.TT0$tjur$Model <- "TT0"
# jsdm.TT1$tjur$Model <- "TT1"
# jsdm.TT2$tjur$Model <- "TT2"

# isdm$tjur <- extract.tjur(isdm$probability)
# isdm$tjur$Model <- "iSDM"
tjur <- bind_rows(isdm$tjur, jsdm.FF0$tjur, jsdm.TF0$tjur, jsdm.TT0$tjur, jsdm.TT1$tjur, jsdm.TT2$tjur)

boxplot(data = tjur, D ~ Model, horizontal = TRUE, xlab="Tjur's R-squared")



# Predictive performance ####
cv.FF0 <- cv.jsdm('outputs/paper 1 extensions', 'FF0')
cv.TF0 <- cv.jsdm('outputs/paper 1 extensions', 'TF0')
cv.TT0 <- cv.jsdm('outputs/paper 1 extensions', 'TT0')
cv.TT1 <- cv.jsdm('outputs/paper 1 extensions', 'TT1')
cv.TT2 <- cv.jsdm('outputs/paper 1 extensions', 'TT2')
cv.TT3 <- cv.jsdm('outputs/paper 1 extensions', 'TT3')
cv.TT4 <- cv.jsdm('outputs/paper 1 extensions', 'TT4')

# Combine all deviance statistics from all models
isdm.cv$deviance$Trial <- "iSDM"
cv.deviance <- bind_rows(cv.FF0$deviance, cv.TF0$deviance, cv.TT0$deviance, cv.TT1$deviance, cv.TT2$deviance, cv.TT3$deviance, cv.TT4$deviance, isdm.cv$deviance[, colnames(cv.FF0$deviance)])

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
map.latent.variable <- function(jsdm){
  if (jsdm$extensions$n.latent==0){
    stop("Latent variables not found in joint model.")
    
  } 
  else if (jsdm$extensions$n.latent==1){
    dt <- tibble(SiteId = names(jsdm$x.lat.maxpost), LatentVariable = jsdm$x.lat.maxpost)
    dt <- left_join(dt, inputs$xy, by="SiteId")
    
    plot.data <- dt
    plot.data$LatentVariableAbs <- abs(plot.data$LatentVariable)
    plot.data$Color <- ifelse(plot.data$LatentVariable > 0, 1, 0)
    plot.data$Color <- as.factor(plot.data$Color)
    
    g <- ggplot()
    g <- g + geom_sf(data = inputs$ch, fill=NA, size=1.25, color="black", show.legend = FALSE)
    g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
    g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)
    
    g <- g + geom_point(data = plot.data, aes(x=X, y=Y, size = LatentVariableAbs, color = Color), alpha=0.35)
    g <- g + labs(title = paste("Latent variable distribution (",unique(jsdm$tjur$Model),")", sep=""),
                  color = "Direction",
                  size = "Magnitude")
    g <- g + theme_void(base_size = 15)
    g <- g + theme(panel.grid.major = element_line(colour="transparent"))
    g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
    g <- g + scale_color_manual(values=c("0" = "#FF0000", "1" = "#0077FF"), labels = c("Negative", "Positive"))
    g <- g + scale_radius(range = c(2, 8))
    
  } 
  else if (jsdm$extensions$n.latent > 1){
    g <- lapply(1:jsdm$extensions$n.latent, function(lv){
      # cat(lv,"\n")
      dt <- tibble(SiteId = rownames(jsdm$x.lat.maxpost), LatentVariable = jsdm$x.lat.maxpost[, lv])
      dt <- left_join(dt, inputs$xy, by="SiteId")

      plot.data <- dt
      plot.data$LatentVariableAbs <- abs(plot.data$LatentVariable)

      g <- ggplot()
      g <- g + geom_sf(data = inputs$ch, fill=NA, size=1.25, color="black", show.legend = FALSE)
      g <- g + geom_sf(data = inputs$rivers.major, fill=NA, color="lightblue", show.legend = FALSE)
      g <- g + geom_sf(data = inputs$lakes.major, fill="lightblue", color="lightblue", show.legend = FALSE)

      g <- g + geom_point(data = plot.data, aes(x=X, y=Y, size = LatentVariableAbs), alpha=0.35)
      g <- g + labs(title = paste("Latent variable ",lv," (",unique(jsdm$tjur$Model),")", sep=""),
                    color = "Direction",
                    size = "Magnitude")
      g <- g + theme_void(base_size = 15)
      g <- g + theme(panel.grid.major = element_line(colour="transparent"))
      g <- g + scale_radius(range = c(2, 8))
      g
    })
  }
  # Return either a ggplot object or a list of ggplot objects
  return(g)
}
g <- map.latent.variable(jsdm.TT2)
ggarrange(g[[1]], g[[2]], common.legend = TRUE)

pdf('outputs/paper 1 extensions/P1 - maps latent variable.pdf', width = 10, paper='a4r')
g
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
# (orangutan attempts to use a telescope)
n <- occur.freq(jsdm.TT1$TT1$occur.taxa)
beta.lat.maxpost <- jsdm.TT1$TT1$beta.lat.maxpost[names(n)]

lv <- sapply(beta.lat.maxpost, function(j){
j * x.lat.maxpost
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
z.FF0$Trial <- "FF0"

z.TF0 <- linear.predictor(jsdm.TF0)
z.TF0$Trial <- "TF0"

z.TT0 <- linear.predictor(jsdm.TT0)
z.TT0$Trial <- "TT0"

z.TT1 <- linear.predictor(jsdm.TT1)
z.TT1$Trial <- "TT1"

z.TT2 <- linear.predictor(jsdm.TT2)
z.TT2$Trial <- "TT2"

z.TT2.mc <- linear.predictor(jsdm.TT2.mc)
z.TT2.mc$Trial <- "TT2 MC"

# z.TT3 <- linear.predictor(jsdm.TT3)
# z.TT3$Trial <- "TT3"
# 
# z.TT4 <- linear.predictor(jsdm.TT4)
# z.TT4$Trial <- "TT4"

z.jsdm <- bind_rows(z.FF0, z.TF0, z.TT0, z.TT1, z.TT2, z.TT2.mc)
rm(z.FF0, z.TF0, z.TT0, z.TT1)

# Get the max and min z-values by variable and by trial
plot.data <- z.jsdm %>%
  group_by(Taxon, Variable, Trial) %>%
  summarise(z.min = min(z, na.rm=T), z.max = max(z, na.rm=T)) %>%
  mutate(z.range = z.max-z.min) %>%
  ungroup()

plot.data$Variable <- factor(plot.data$Variable, levels=c("Temp", "FV", "F100m", "IAR", "Urban", "LUD", "Site effect", "TT1", "TT2", "TT3", "TT4"))

g <- ggplot(data = plot.data) 
g <- g + geom_boxplot(aes(y=z.range, fill=Variable))
g <- g + facet_grid(. ~ Trial)
g <- g + coord_cartesian(ylim=c(0,16)) 
g <- g + theme_bw(base_size = 24)
g <- g + labs(title="Distribution of range of linear predictor for each taxon by JSDM trial",
              y = expression(paste("z"["range,kj"])),
              fill = "Variable")
g <- g + theme(axis.title.x=element_blank(),
               axis.text.x=element_blank(),
               axis.ticks.x=element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.minor = element_blank())
g <- g + scale_fill_brewer(palette = "Set1")

pdf('outputs/paper 1 extensions/P1 - slopes sample-based LV.pdf', width = 15)
g
dev.off()

# Test if latent variable (max. post) is correlated with the influence factors
# plot.data <- p[, c("SiteId", inf.fact[inf.fact != "Temp2"])]
# plot.data <- plot.data[plot.data$SiteId %in% names(jsdm.TT1$x.lat.maxpost),]
plot.data <- predictors

plot.data %>% group_by(SiteId) %>% summarise(n=n()) %>% View()
plot.data <- unique(plot.data)

plot.data$x.lat.maxpost <- jsdm.TT1$x.lat.maxpost[plot.data$SiteId]

plot.data <- as.matrix(select(plot.data, -SiteId, -SampId))

pairs.panels(plot.data, main="Environmental variables (mean-centered) and latent variable")


# Correlated beta_comm ####
K <- intersect(jsdm.FF0$FF0$inf.fact, jsdm.TF0$TF0$inf.fact)

beta.comm <- lapply(K, function(k){
  mu <- jsdm.FF0$FF0$mu.beta.comm.maxpost[k]
  sigma <- jsdm.FF0$FF0$sigma.beta.comm.maxpost[k]
  beta.comm.FF0 <- rnorm(1000, mu, sigma)
  
  data1 <- tibble(Trial = "FF0", Variable = k, Value = beta.comm.FF0)

  mu <- jsdm.TF0$TF0$mu.beta.comm.maxpost[k]
  sigma <- jsdm.TF0$TF0$sigma.beta.comm.maxpost[k]
  beta.comm.TF0 <- rnorm(1000, mu, sigma)
  data2 <- tibble(Trial = "TF0", Variable = k, Value = beta.comm.TF0)
  
  data <- bind_rows(data1, data2)
})

beta.comm <- rbindlist(beta.comm)

# See Github thread below for justification of this code:
# https://github.com/tidyverse/tidyr/issues/426
beta.comm <- beta.comm %>%
  group_by(Trial, Variable) %>%  # group by everything other than the value column. 
  mutate(row_id=1:n()) %>% 
  ungroup() %>%  # build group index
  spread(key=Variable, value=Value) %>%    # spread
  select(-row_id)  # drop the index

colors <- c("FF0"="red", "TF0"="blue")

beta.comm <- select(beta.comm, Trial, Temp, Temp2, FV, F100m, IAR, LUD)
# colnames(beta.comm)[3] <- "Temp^2"

pairs(beta.comm[, K], pch=19,cex=0.5, col = colors[beta.comm$Trial], lower.panel = NULL)
legend("bottomleft", fill = c("red", "blue"), col = c("red", "blue"), legend = c("FF0", "TF0"))

panel.ellipse <- function(x, y, ...){
  ellipse()
  # abline(v=0,h=0)
}
# panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y))
#   txt <- format(c(r, 0.123456789), digits = digits)[1]
#   txt <- paste0(prefix, txt)
#   if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor * r)
# }

# library(MASS)
# mu.TF0 <- c(jsdm.TF0$TF0$mu.beta.comm.maxpost["Temp"], jsdm.TF0$TF0$mu.beta.comm.maxpost["FV"])
# sigma.TF0 <- matrix(c(jsdm.TF0$TF0$sigma.beta.comm.maxpost["Temp"],
#                       jsdm.TF0$TF0$sigma.beta.comm.maxpost["FV"],
#                       jsdm.TF0$TF0$sigma.beta.comm.maxpost["FV"],
#                       jsdm.TF0$TF0$sigma.beta.comm.maxpost["Temp"]), 2)
# # sigma.TF0 <- cov(jsdm.TF0$TF0$sigma.beta.comm.maxpost["Temp"], jsdm.TF0$TF0$sigma.beta.comm.maxpost["FV"])
# bivn <- mvrnorm(1000, mu = mu.TF0, Sigma = sigma.TF0)
# 
# bivn <- mvrnorm(1000, mu = c(0, 0),
#                 Sigma = matrix(c(1, .5, .5, 1), 2))
# bivn2 <- mvrnorm(1000, mu = c(-2, 2),
#                  Sigma = matrix(c(1.5, 1.5, 1.5, 1.5), 2))
# 
# # mu.TF0 <- jsdm.TF0$TF0$mu.beta.comm.maxpost
# # sigma.TF0 <- jsdm.TF0$TF0$sigma.beta.comm.maxpost




# x <- jsdm.FF0$FF0$env.cond[, inf.fact]
# y <- as.matrix(jsdm.FF0$FF0$occur.taxa)
# 
# grfo.act <- inf.fact
# mu.beta.grfo.maxpost <- mu.FF0
# sigma.beta.grfo.maxpost <- sigma.FF0
# beta.taxa.maxpost <- jsdm.FF0$FF0$beta.taxa.maxpost
# 
# data <- cbind(1:ncol(x),rep(ncol(x)),rep(length(grfo.act),ncol(x)),mu.beta.grfo.maxpost,sigma.beta.grfo.maxpost)
# for ( i in 1:length(grfo.act) ) data <- cbind(data,R.beta.grfo.maxpost[i,,])
# data <- cbind(data/1e6,beta.taxa.maxpost)
# 
# panel.ellipse <- function(x,y,...)
# {
#   i      <- floor(1e6*x[1]+0.1)
#   j      <- floor(1e6*y[1]+0.1)
#   n.pred <- floor(1e6*x[2]+0.1)
#   n.grfo <- 1
#   n      <- length(x)
#   points(x[(3+2*n.grfo+n.grfo*n.pred+1):n],y[(3+2*n.grfo+n.grfo*n.pred+1):n],...)
#   abline(v=0,h=0)
#   for ( k in 1:n.grfo )
#   {
#     e <- ellipse(x      = x[3+2*n.grfo+(k-1)*n.pred+j]*1e6,
#                  scale  = c(x[3+n.grfo+k],y[3+n.grfo+k])*1e6,
#                  centre = c(x[3+k],y[3+k])*1e6,
#                  level  = 0.9)
#     lines(x=e[,1],y=e[,2],col=c("blue","red","green","orange","violet","black")[k])
#   }
# }
# 
# pairs(t(data),
#       main="FF0 vs TF0 community parameters",
#       pch=19,cex=0.3,
#       lower.panel=panel.ellipse,upper.panel=panel.ellipse)