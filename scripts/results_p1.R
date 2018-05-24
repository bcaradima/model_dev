### RESULTS: Paper 1 ####
# Description: this script prepares the results for paper 1. It processes the joint model workspaces and combines the results with those of the individual models before producing various plots.

# Calibration results ####
# Results: paper 1 BDM species
jsdm.re <- extract.jsdm('outputs/jsdm_p1', 'RE', epsilon=TRUE)
jsdm <- extract.jsdm('outputs/jsdm_p1', 'noRE', epsilon=FALSE)

### > Combine i/jSDM results ####
# Warning: manually check columns being dropped!
probability <- bind_rows(isdm.prob, jsdm$probability[, colnames(isdm.prob), with=FALSE])
deviance <- bind_rows(isdm.dev, select(jsdm$deviance, -Trial))

# > i/jSDM cross-validation ####
jsdm.cv <- cv.jsdm('noRE')
# Combine CV deviance ####
# Warning: manually check columns being dropped
deviance.cv <- bind_rows(isdm.cv.dev, select(jsdm.cv$deviance, -Trial))

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
pdf('P1 - BDM occurrence frequency.pdf', height = 13, width = 15)
par(cex=2.5)
plot(x = 1:length(n.bdms), y = n.bdms, axes=F, xlab="Taxa (sorted from highest to lowest occurrence frequency)", ylab="Occurrence frequency")
axis(1, c(0,50,100,150,200,247))
axis(2, c(0,100,200,300,400,500,581))
dev.off()




# Plot P1 - deviance vs occurrence ####
deviance$rel.freq <- deviance$n/deviance$n.samples

test1 <- spread(select(deviance, Taxon, Model, std.res.dev), Model, std.res.dev) # y coordinates
colnames(test1) <- c("Taxon", "iSDM.dev", "jSDM.dev")

test2 <- spread(select(deviance, Taxon, Model, rel.freq), Model, rel.freq) # x coordinates
colnames(test2) <- c("Taxon", "iSDM.rf", "jSDM.rf")

deviance.segments.fig2 <- left_join(test1, test2, by = "Taxon")

null.dev <- deviance[Model == "iSDM", ]
null.dev$Model <- "Null model"

n <- 581
m <- 1:(n-1)
par(cex=2.25)
plot(m/n*100,-2*(m/n*log(m/n)+(n-m)/n*log(1-m/n)), 
     type="l", lwd=2,
     xlab="Occurrence frequency (%)",
     ylab="Null model deviance")

null.dev.line <- data.table(null.dev = -2*(m/n*log(m/n)+(n-m)/n*log(1-m/n)), n.samples = m)
null.dev.line$rel.freq <- m/n
null.dev.line$Model <- "Null model"

pdf('P1 - deviance vs occurrence.pdf', paper='special', width = 9, height = 6)
g <- ggplot(data=deviance)
g <- g + geom_point(data = null.dev, aes(x=rel.freq, y=null.dev/n.samples, color=Model), size = 3, alpha=0.4)
g <- g + geom_line(data = null.dev.line, aes(x=rel.freq, y=null.dev, color=Model), alpha=0.4, show.legend = FALSE)
g <- g + geom_point(aes(x=rel.freq, y=std.res.dev, color = Model), size=3, alpha=0.4)
g <- g + geom_segment(data = deviance.segments.fig2, aes(x = iSDM.rf, y = iSDM.dev, xend = jSDM.rf, yend= jSDM.dev), alpha = 0.3)

g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF", "Null model" = "#000000"))
g <- g + theme_bw(base_size = 18)
g <- g + theme(axis.text=element_text(size=18))
g <- g + theme(plot.title = element_blank())
g <- g + labs(title = "Model deviance vs relative occurrence frequency",
              y = "Standardized deviance",
              x = "Relative occurrence frequency")
g <- g + guides(colour = guide_legend(override.aes = list(size=6)), shape = FALSE, color = FALSE)
print(g)
dev.off()
rm(test1, test2, deviance.segments.fig2, null.dev, null.dev.line, n, m)

# Plot P1 - D2 vs std dev ####
# Point segments
test1 <- spread(select(deviance, Taxon, Model, D2), Model, D2) # y coordinates
colnames(test1) <- c("Taxon", "iSDM.D2", "jSDM.D2")
test2 <- spread(select(deviance, Taxon, Model, std.res.dev), Model, std.res.dev) # x coordinates
colnames(test2) <- c("Taxon", "iSDM.resdev", "jSDM.resdev")
deviance.segments <- left_join(test1, test2, by = "Taxon")

pdf('P1 - D2 vs std dev.pdf', paper='special', width = 13, height = 9)
g <- ggplot(data = deviance)
g <- g + geom_segment(data = deviance.segments, aes(x = iSDM.resdev, y = iSDM.D2, xend = jSDM.resdev, yend= jSDM.D2), alpha = 0.2)
g <- g + geom_point(data = deviance, aes(x = std.res.dev, y = D2, color = Model, size = n), alpha = 0.4)
g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
g <- g + theme_bw(base_size = 18)
g <- g + theme(axis.text=element_text(size=18))
g <- g + theme(plot.title = element_blank())

g <- g + labs(title = "",
              x = expression(paste("Standardized d"["j"]^"proposed")),
              y = expression("D"["j"]^2),
              colour = "Model", 
              size = "Occurrence\nfrequency")
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
print(g)
dev.off()
rm(test1, test2, deviance.segments)

# Plot P1 - cross-validation ####
# Use dplyr verbs to shape plotting data for ggplot
cv.plot <- deviance.cv %>%
  select(Taxon, Type, Model, std.deviance) %>%
  group_by(Taxon, Type, Model) %>%
  summarise(MRD = mean(std.deviance)) %>%
  spread(Type, MRD) %>%
  ungroup() %>% # remove grouping information
  mutate(n = n.bdms[Taxon]) # add occurrence frequency again


# Format infinite predictive deviance taxa for plotting
cv.plot$Shape <- ifelse(is.infinite(cv.plot$Testing), 24, 16)
# cv.plot$n <- ifelse(is.infinite(cv.plot$Testing), NA, cv.plot$n)
cv.plot$Shape <- as.factor(cv.plot$Shape)
cv.plot$Testing <- ifelse(is.infinite(cv.plot$Testing), max(cv.plot$Testing[!is.infinite(cv.plot$Testing)]), cv.plot$Testing)

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
              size = "Occurrence\nfrequency")
g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
g <- g + scale_size_continuous(range = c(2, 8))
g <- g + scale_shape_discrete(name  = "Deviance",
                              breaks=c("16", "24"),
                              labels=c("In range", "Out of range"))

pdf('P1 - cross-validation.pdf', height=8, width=10.5)
print(g)
dev.off()

# Plot P1 - slopes jSDM ####
slopes <- linear.predictor(jsdm)

# Match expressions to influence factors
slopes.plot <- slopes
v <- jsdm$noRE$inf.fact[jsdm$noRE$inf.fact != "Temp2"]
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
g <- g + labs(x = "", y = expression(paste(Delta, "z")))
g <- g + guides(alpha = FALSE)

pdf('P1 - slopes jSDM.pdf', width = 13, height = 9, onefile=TRUE)
print(g)
dev.off()

# Plot P1 - parameters jSDM [all taxa] ####
plot.comm(jsdm, 'P1 - parameters jSDM [all taxa]')

# Plot P1 - significant responses ####
bdms.stan <- select.jsdm(jsdm) # retrieve Stan.fit object
beta.samples <- bdms.stan[["beta_taxa"]]
dimnames(beta.samples) <- list(1:dim(beta.samples)[1], jsdm$noRE$inf.fact, colnames(bdms.stan$occur.taxa)) # name the dimensions
taxon.samples <- melt(beta.samples)
colnames(taxon.samples) <- c("Sample", "Variable", "Taxon", "Value")
taxon.samples <- setDT(taxon.samples)
taxon.samples$Variable <- as.character(taxon.samples$Variable)
taxon.samples$Taxon <- as.character(taxon.samples$Taxon)

beta.taxa.response <- matrix(nrow=length(n.bdms), ncol=length(inf.fact))

for (k in 1:length(inf.fact)){
  variable <- inf.fact[k]
  samples <- taxon.samples[Variable == variable,]
  
  for (j in 1:length(names(n.bdms))){
    taxon <- names(n.bdms)[j]
    sample <- samples[Taxon==taxon,]
    sample <- sample$Value
    
    # Fill matrix: beta.taxa.response
    jsig <- quantile(sample, probs = c(0.05, 0.95))
    
    # If 5th quantile greater than 0, set positive
    # If 95th quantile less than 0, set negative
    if (jsig[1] > 0){ # significant positive
      response <- "blue"
    }
    if (jsig[2] < 0){ # significant negative
      response <- "red"
    }
    # If posterior is !(positive) AND !(negative), set grey
    if (!(jsig[1] > 0) & !(jsig[2] < 0)){
      response <- "grey55"
    }
    beta.taxa.response[j, k] <- response
  }
}

# Based on quality-of-fit and predictive performance, find taxa of interest that:
test <- cv.plot %>%
  filter(Model=="jSDM") %>%
  select(Taxon, Training, Testing)

# Combine goodness-of-fit and predictive performance
taxa.int <- deviance %>%
  filter(Model=='jSDM', D2 > 0.25 & D2 < 0.60, std.res.dev > 0.15 & std.res.dev < 0.90) %>%
  select(Taxon, D2) %>%
  left_join(test, by = "Taxon") %>%
  mutate(n=n.bdms[Taxon], TrTe = Training/Testing)
rm(test)

# Select taxa for heatmap
d.select <- spread(select(deviance, Taxon, Model, D2), Model, D2)
colnames(d.select) <- c("Taxon", "iSDM.D2", "jSDM.D2")

std.dev.select <- spread(select(deviance, Taxon, Model, std.res.dev), Model, std.res.dev)
colnames(std.dev.select) <- c("Taxon", "iSDM.dev", "jSDM.dev")

dt <- left_join(d.select, std.dev.select, by = "Taxon")
setDT(dt)
rm(d.select, std.dev.select)

# joint model will always have higher residual deviance; what % higher?
dt$rdev.diff <- (dt$iSDM.dev/dt$jSDM.dev)*100

# joint model will always have lower D2; what % 
dt$d.diff <- (dt$iSDM.D2/dt$jSDM.D2)*100

# selection criteria;
# good fit in joint model,
# iSDM deviance >95% of jSDM deviance,
# iSDM D2 not more than 10% of jSDM D2
nrow(filter(dt, jSDM.D2 > 0.30 & rdev.diff > 0.95 & d.diff < 110))

# Format density plot output for heatmap
beta.taxa.response <- as.data.table(beta.taxa.response)
beta.taxa.response$Taxon <- names(n.bdms)
beta.taxa.response$n <- n.bdms[beta.taxa.response$Taxon]
colnames(beta.taxa.response) <- c(inf.fact, "Taxon", "n")
beta.taxa.response <- arrange(beta.taxa.response, desc(n))

setDT(beta.taxa.response)
beta.taxa.response[beta.taxa.response=="red"] <- -1
beta.taxa.response[beta.taxa.response=="blue"] <- 1
beta.taxa.response[beta.taxa.response=="grey55"] <- 0

# Color coding for all taxa
# Subset the taxa with good model performance
btr.int <- filter(beta.taxa.response, Taxon %in% taxa.int$Taxon)
btr.int <- btr.int[, c("Taxon", "n", inf.fact)]
btr.plot <- btr.int[btr.int$n > 56, c(inf.fact)]
btr.plot <- apply(btr.plot, 2, as.numeric)
row.names(btr.plot) <- paste(sub("_", " ", btr.int$Taxon[btr.int$n > 56]), " - ", n.bdms[btr.int$Taxon[btr.int$n > 56]])

# Select the number of clusters
d <- dist(btr.plot, method="euclidean")
pfit <- hclust(d, method="ward.D")
plot(pfit, labels=row.names(btr.plot))
rect.hclust(pfit, k=4)

library(pheatmap)
# Save manually
pheatmap(btr.plot, col=c("tomato3", "snow3", "royalblue4"), cellwidth=30, cellheight=11, cluster_rows=T, cutree_rows=5, cluster_cols=F, clustering_distance_rows = "euclidean", legend=F, show_rownames=TRUE, fontsize=13, fontsize_col=15)


# Plot P1 - significant responses [all taxa] ####
btr.all <- beta.taxa.response
btr.all <- as.data.table(btr.all)
btr.all <- arrange(btr.all, desc(n))
btr.all$Label <- factor(paste(sub("_", " ", btr.all$Taxon), " - ", btr.all$n), levels=paste(sub("_", " ", btr.all$Taxon), " - ", btr.all$n))
plot.data <- gather(btr.all, Variable, Value, -Taxon, -Label, -n)

g <- ggplot(plot.data, aes(x = Variable, y = Label))
g <- g + geom_tile(aes(fill = Value), colour = "white")
g <- g + scale_fill_manual(values=c("tomato3", "snow3", "royalblue4"))
g <- g + theme_minimal(base_size = 15)
pdf('jSDM heatmap ALL.pdf', paper='special', height=56, width=8.5)
print(g)
dev.off()


# Plot P1 - parameters jSDM [example taxa] ####
pdf('jSDM parameter distributions [example taxa].pdf', onefile = TRUE)
par(cex=1.25)
for (k in 1:length(jsdm$noRE$inf.fact)){
  variable <- jsdm$noRE$inf.fact[k]
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
              "Temp2" = expression(paste(beta["Temp"^2], " (1/", degree, "C)")),
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
    
    lines(density(sample), type="l", col=alpha(c, 1), lwd=1.25)
    text(x=mean(sample), y=max(density(sample)$y)+0.1*max(density(sample)$y), labels = l[j])
  }
  # Plot community parameter distribution
  lines(x, beta.comm.density, type="l", col = "grey50", lwd = 5)
  # Plot maximum posterior values over all taxa
  lines(beta.taxa.maxpost.density, type="l", col="black", lwd=2, lty='longdash')
  
}
dev.off()

# Plot P1 - prob jSDM [labelled]
plot.prob(isdm, 'P1 - prob vs inputs iSDM')
plot.prob(jsdm, 'P1 - prob vs inputs jSDM')

# Plot P1 - maps ijSDM ####
map.isdm(isdm, 'P1 - maps iSDM')
map.jsdm(jsdm, 'P1 - maps jSDM')

# Plot P1 - maps influence factors ####
x <- prepare.inputs(K, sample.bdms, center = FALSE)
x <- gather(x, Variable, Value, -SiteId, -SampId)
x <- filter(x, Variable != "Temp2")
x$Trial <- 'BDM species'
map.inputs(x, 'P1 - maps influence factors')



# Plot P1 - parameter SD ####
# Extract the entire posterior beta sample
beta.samples <- jsdm$noRE$beta_taxa
bs <- melt(beta.samples, varnames = c("Sample", "Variable", "Taxa"), value.name = "Value", as.is = TRUE)
colnames(bs) <- c("Sample", "Variable", "Taxon", "Parameter")
bs$Variable <- as.character(bs$Variable)
bs$Taxon <- as.character(bs$Taxon)
bs <- tbl_df(bs)
rm(beta.samples)

# Fast aggregation of 10M row dataset
bssd <- bs %>%
  group_by(Taxon, Variable) %>%
  summarise(SD = sd(Parameter), Mean = mean(Parameter))

bssd$n <- n.bdms[bssd$Taxon]

bssd$Label <- factor(bssd$Variable, levels = inf.fact)
levels(bssd$Label) <- labeller(levels(bssd$Label))

g <- ggplot(data=bssd, aes(x = n, y = SD, size = n))
g <- g + geom_point(alpha = 0.5)
g <- g + facet_grid(Label ~ ., scales = "free", labeller=label_parsed)
g <- g + theme_bw(base_size = 14)
g <- g + theme(strip.background=element_rect(fill="grey"),strip.text=element_text(color="black", face="bold"),
               plot.title = element_text(hjust = 0.5, size = 12))
g <- g + labs(title = expression(paste("Standard deviation of posterior taxon-specific parameter distributions ", beta["jk"]^"taxa")),
              x = "Occurrence frequency",
              y = expression(paste("Standard deviation (", sigma[beta["jk"]^"taxa"],")")),
              size = "Occurrence\nfrequency")
g <- g + scale_y_continuous(limits=c(0,NA))
pdf('P1 - parameter SD.pdf', height = 12.5)
print(g)
dev.off()

# Plot P1 - parameter uncertainty ####
bssd$rSD <- bssd$SD/bssd$Mean
g <- ggplot(data=bssd, aes(x = n, y = rSD, size = n))
g <- g + geom_point(alpha = 0.5)
g <- g + facet_wrap(~ Variable, labeller=label_parsed)
g <- g + theme_bw(base_size = 13)
g <- g + theme(strip.background=element_rect(fill="black"), strip.text=element_text(color="white", face="bold"), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
g <- g + labs(title = expression(paste(sigma,"/",mu)),
              x = expression(paste("Mean (", mu[beta["jk"]^"taxa"],")")),
              y = expression(paste("Standard deviation/mean (", sigma[beta["jk"]^"taxa"],"/", mu[beta["jk"]^"taxa"],"))")),
              size = "Total occurrences")
print(g)


bssd.cv <- bssd %>%
  mutate(rSD = SD/abs(Mean))
bssd.cv$Label <- factor(bssd.cv$Variable, levels = c("Temp", "Temp2", "FV", "F10m", "IAR", "Urban", "LUD"))

g <- ggplot(data=bssd.cv)
g <- g + geom_boxplot(aes(Label, rSD, fill=Label))
g <- g + coord_cartesian(ylim=c(0,10))
g <- g + theme_bw(base_size = 15)
g <- g + theme(strip.background=element_rect(fill="black"), strip.text=element_text(color="white", face="bold"),
               plot.title = element_text(hjust = 0.5, size = 12), axis.text = element_text(size = 12))
g <- g + scale_fill_brewer(palette = "Set1")
g <- g + guides(fill=FALSE)
g <- g + labs(title =  expression(paste("Relative uncertainty in the posterior taxon-specific parameter distributions ", beta["jk"]^"taxa")),
              x = "Influence factor",
              y = expression(paste("Relative standard deviation (", sigma[beta["jk"]^"taxa"]," / |",mu[beta["jk"]^"taxa"],"|)")))
pdf('P1 - parameter uncertainty.pdf')
print(g)
dev.off()

# Table P1 - Count significant responses ####
# Count significant parameters
p1.uncertainty <- bs %>%
  group_by(Variable, Taxon) %>%
  do(data.frame(t(quantile(.$Parameter, probs = c(0.05, 0.95))))) %>%
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

g <- ggplot(inputs$ch, aes(POINT_X, POINT_Y))
g <- g + geom_path(lineend = "round")
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
# Extract the entire posterior beta sample
beta.samples <- jsdm$noRE$beta_taxa
# Subset the parameter sample based on lower and upper percentiles
beta.samples <- apply(beta.samples, c(2,3), function(x){
  lower <- quantile(x, probs = 0.10)
  upper <- quantile(x, probs = 0.90)
  x <- c(lower, upper)
})

isdm.beta <- isdm.bdms$parameters
if ("Model" %in% colnames(isdm.beta)){
  isdm.beta$Model <- NULL
}
colnames(isdm.beta) <- c("Taxon", "Variable", "iSDM.parameter")
isdm.beta$Variable <- as.character(isdm.beta$Variable)

# Number of pages per taxa (ntpp)
pages <- 4
ntpp <- length(n)/pages
tpp <- split(n, ceiling(seq_along(n)/ntpp)) # bin the number of taxa per page
pdf(paste(output_dir, "/", trials[t], "_dotplot_notnorm.pdf",sep=""), paper = "special", height = h, width = w, onefile = TRUE)
for (page in 1:pages){
  # Tidy the jSDM beta quantiles
  # Melt 3D array to 2D array
  bs <- melt(beta.samples, varnames = c("Sample", "Variable", "Taxon"), value.name = "Value", as.is = TRUE)
  bs <- tbl_df(spread(bs, Sample, Value))
  colnames(bs)<- c("Variable", "Taxon", "quantile10", "quantile90")
  
  bs <- bs[bs$Taxon %in% names(tpp[page][[1]]), ]
  
  # Prepare the jSDM parameter values
  beta.max <- data.table(t(beta.taxa.maxpost), stringsAsFactors = F)
  colnames(beta.max) <- inf.fact
  beta.max$Taxon <- colnames(occur.taxa)
  beta.max <- gather(beta.max, Variable, Value, -Taxon)
  colnames(beta.max) <- c("Taxon", "Variable", "jSDM.parameter")
  beta.max$Variable <- as.character(beta.max$Variable)
  
  # Join the quantiles, jSDM, and iSDM parameters together
  bs <- left_join(bs, beta.max, by = c("Taxon", "Variable"))
  bs <- left_join(bs, isdm.beta, by = c("Taxon", "Variable"))
  
  # Model (variable) and Parameter (value) gather the two columns, jSDM.parameter and iSDM.parameter
  # into a variable/value pair
  bs <- gather(bs, Model, Parameter, jSDM.parameter:iSDM.parameter)
  
  # Process iSDM outliers
  # Get min and max posterior parameter values within the 10th-90th quantiles
  lim.min <- aggregate(quantile10 ~ Variable, bs, min)
  vmin <- lim.min$quantile10; names(vmin) <- lim.min$Variable
  
  lim.max <- aggregate(quantile90 ~ Variable, bs, max)
  vmax <- lim.max$quantile90; names(vmax) <- lim.max$Variable
  
  # Set colors for each parameter
  # g <- g + scale_colour_manual(values=c("jSDM" = "#048504", "iSDM" = "#790FBF"))
  bs$Col <- "#048504"
  bs$Col[bs$Model == 'iSDM.parameter'] <- "#790FBF"
  bs$Col <- ifelse(bs$Parameter < vmin[bs$Variable], "#6B6B6B", bs$Col)
  bs$Parameter <- ifelse(bs$Parameter < vmin[bs$Variable], vmin[bs$Variable], bs$Parameter)
  
  bs$Col <- ifelse(bs$Parameter > vmax[bs$Variable], "#6B6B6B", bs$Col)
  bs$Parameter <- ifelse(bs$Parameter > vmax[bs$Variable], vmin[bs$Variable], bs$Parameter)
  rm(lim.min, lim.max, vmin, vmax)
  
  # Order labels by occurrence frequency
  bs$Labels <- factor(paste(bs$Taxon, ' - ', n[bs$Taxon]), levels = paste(names(sort(n)), ' - ', sort(n)))
  
  # Order variable facets and pass expressions for units in facet labels
  bs$Variable <- factor(bs$Variable, levels = c(inf.fact))
  levels(bs$Variable) <- c(expression(paste("Temp (1/", degree, "C)")),
                           expression(paste("Temp"^2," (1/", degree, "C)")),
                           expression(paste("FV (s/m)")),
                           expression(paste("F10m (1/%)")),
                           expression(paste("IAR")), # units don't fit: (spray treatments * fraction cropland)
                           expression(paste("Urban (1/%)")),
                           expression(paste("LUD (CE/km"^2,")")))
  # Build the plot
  g <- ggplot(bs)
  g <- g + geom_hline(yintercept=0, alpha=0.4)
  g <- g + geom_linerange(aes(x = Labels, ymin = quantile10, ymax = quantile90), color = 'black', alpha = 0.6)
  g <- g + geom_point(aes(x = Labels, y = Parameter, colour = Col), stroke=0, size = 2.5, alpha = 0.5)
  g <- g + facet_grid(. ~ Variable, scales = "free", labeller=label_parsed)
  g <- g + coord_flip()
  g <- g + theme_bw()
  g <- g + labs(title = "Taxon-specific parameters in the individual and joint models",
                subtitle = "Maximum likelihood estimates (iSDM) and maximum posterior values with 10th-90th percentile interval (jSDM)",
                x = "Occurrence frequency and taxon",
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
taxa <- n[taxa]
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
  g <- g + ggtitle(paste("Page ", page, " - probability for training/testing of iSDM vs jSDM for ", paste(tname), ' (n = ', n[tname], ')', sep = ''))
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
#   select(Species, std.res.dev, RE) %>%
#   spread(RE, std.res.dev) %>%
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