# r <- 56 # 56, 180, 102 for BDMs/BDMf, Invf, Invf_plateau
# Extract VS results ####
vs.bdms <- extract.vsr("variable_selection/bdms_transformed", sample.bdms)

# Get mean relative deviance for each model over entire community (excluding rare taxa)
vs.bdms.mean <- vs.bdms %>%
  filter(!is.infinite(rdev.test) & n > 56) %>%
  group_by(Model) %>%
  summarise(mrd.train = mean(rdev.train, na.rm=TRUE), mrd.test = mean(rdev.test, na.rm=TRUE)) %>%
  arrange(mrd.test) %>%
  mutate(Parameters = vapply(strsplit(Model, " "), length, integer(1)))

View(vs.bdms.mean)

# write.csv(vs.bdms.mean, 'outputs/variable_selection/bdms_transformed/vs_bdms_mean.csv', row.names = FALSE)

vs.bdms.ut <- extract.vsr("variable_selection/bdms_untransformed", sample.bdms)

# Get mean relative deviance for each model over entire community (excluding rare taxa)
vs.bdms.mean.ut <- vs.bdms.ut %>%
  filter(!is.infinite(rdev.test) & n > 56) %>%
  group_by(Model) %>%
  summarise(mrd.train = mean(rdev.train, na.rm=TRUE), mrd.test = mean(rdev.test, na.rm=TRUE)) %>%
  arrange(mrd.test) %>%
  mutate(Parameters = vapply(strsplit(Model, " "), length, integer(1)))

View(vs.bdms.mean.ut)

# write.csv(vs.bdms.mean.ut, 'outputs/variable_selection/bdms_untransformed/vs_bdms_mean_ut.csv', row.names = FALSE)

test <- filter(vs.bdms.mean.ut, !grepl('SS1|SS2|SS3', Model))

# d.mean <- filter(d.mean, Parameters >= 4 & Parameters <= 6)
# d.mean.invf46 <- d.mean

# write.csv(d.mean, paste('outputs',vs,'meandev_models.csv',sep='/'))

# rare.occurrences <- c(7, 20, 30, 40, 56)
# pdf('variable_frequency_rarity_BDMSp.pdf', paper='special', width = 14, height = 8.5)
# for (r in 1:length(rare.occurrences)){
#   rarity <- rare.occurrences[r]
#   d.mean <- d %>%
#     filter(!is.infinite(rdev.test), n > rarity) %>%
#     group_by(Model) %>%
#     summarise(mrd.train = mean(rdev.train, na.rm=TRUE), mrd.test = mean(rdev.test, na.rm=TRUE)) %>%
#     arrange(mrd.test) %>%
#     mutate(Parameters = vapply(strsplit(Model, " "), length, integer(1)), n_excl = rarity)
#   models <- d.mean[1:20, ]$Model
#   v <- unlist(strsplit(models, " "), use.names = FALSE)
#   barplot(sort(table(v), decreasing = TRUE),
#           main = paste("Variable frequency in top", length(models),"models for BDM species"))
#   # write.csv(d.mean[1:20, ], paste('top_models_BDMSp', '_', rarity, '.csv', sep=''))
# }
# dev.off()

### DEVIANCES ####
plot.data <- vs.bdms.mean %>%
  arrange(Parameters) %>%
  group_by(Parameters) %>%
  mutate(Label = paste(Parameters, " parameters (", formatC(uniqueN(Model), big.mark=",")," models)", sep="")) %>%
  ungroup() %>%
  mutate(Label = factor(Label, levels = unique(Label)))

plot.data.points <- plot.data %>%
  group_by(Parameters) %>%
  summarise(x = mean(mrd.train), y = mean(mrd.test))
label <- levels(plot.data$Label)
names(label) <- 1:uniqueN(plot.data$Label)
plot.data.points$Label <- label[as.character(plot.data.points$Parameters)]
plot.data.points$Label <- factor(plot.data.points$Label, levels = unique(plot.data.points$Label))
rm(label)

# Other plots ####
pdf(paste('outputs',vs,'Model_selection_plots.pdf',sep='/'), width = 14, height = 9, onefile = TRUE)
m <- sapply(output, function(m){nrow(m$models.to.run)})
plot(np, m, type = "b",
     main="Number of models to cross-validate versus model complexity \n (models with collinear variables are not cross-validated)",
     xlab="Number of parameters",
     ylab=paste("Number of models run (total = ", sum(m),")", sep=""))

d.mean$PercentDecrease <- 1-(d.mean$mrd.train/d.mean$mrd.test)

g <- ggplot(d.mean, aes(x = Label, y = PercentDecrease))
g <- g + geom_boxplot()
g <- g + theme_bw(base_size = 17)
g <- g + labs(title = "Predictive performance versus model complexity",
          y = "% change in MSD for prediction over calibration",
          x = "Number of parameters")
g <- g + theme(axis.text.x=element_text(angle=25,hjust=1))
print(g)

g <- ggplot(d.mean)
g <- g + geom_point(aes(x = mrd.train, y = mrd.test), alpha = 0.20)
g <- g + geom_density2d(aes(x = mrd.train, y = mrd.test), color= "red")
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
g <- g + geom_point(data = d.mean.points, aes(x = x, y = y), color = "blue", size = 4)
# g <- g + facet_grid(Parameters ~ .)
g <- g + facet_wrap(~ Label, ncol = 4)
g <- g + theme_bw(base_size = 17)
g <- g + labs(title = "Mean standardized deviance (MSD) of prediction vs calibration per model",
              subtitle = "Grouped by parameters (MSD based on 3-fold CV per taxon)",
              x = "Mean standardized deviance for calibration",
              y = "Mean standardized deviance for prediction")
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
print(g)

plot.data <- filter(plot.data, Parameters < 11)
g <- ggplot(plot.data, aes(x = mrd.train, y = mrd.test, color = as.factor(Label)))
g <- g + geom_point(alpha=0.35)
# g <- g + geom_density2d(aes(x=mrd.train, y=mrd.test, color = as.factor(Label)))
g <- g + geom_abline(intercept = 0, slope = 1, color="black", size=1.25, alpha = 0.4)
g <- g + theme_bw(base_size = 17)
g <- g + labs(title = "Predictive deviance vs calibration deviance",
              y = "Mean standardized deviance during prediction",
              x = "Mean standardized deviance during calibration",
              color = "Number of Parameters")
g <- g + theme(plot.title=element_blank())
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + scale_colour_brewer(palette = "Set1")
print(g)

# How often do individual influence factors appear in the top 20 models?
# Barplot: variable frequency in top models ###
n.models <- 10
d.mean <- arrange(d.mean, mrd.test)
models <- d.mean$Model[1:n.models]
v <- unlist(strsplit(models, " "), use.names = FALSE)
barplot(sort(table(v), decreasing = TRUE),
        main = paste("Variable frequency in top", n.models,"models"))

g <- ggplot(d.mean, aes(x = mrd.train, y = PercentDecrease, color = as.factor(Label)))
g <- g + geom_point(alpha=0.35)
# g <- g + geom_smooth(method="lm")
g <- g + theme_bw(base_size = 17)
g <- g + labs(title = "Percent decrease of predictive deviance vs calibration deviance",
              y = "% change in mean standardized deviance for prediction",
              x = "Mean standardized deviance for calibration",
              color = "Number of Parameters")
g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
g <- g + scale_colour_brewer(palette = "Set1")
print(g)
dev.off()

# Group rows by number of parameters and get bottom n rows by lowest mean relative deviance for predictions
top.models. <- d.mean %>%
  group_by(Parameters) %>%
  arrange(mrd.test)

top.models.variety <- d.mean %>%
  group_by(Parameters) %>%
  top_n(-3, mrd.test) %>%
  arrange(mrd.test)

# What is the increase in deviance from calibration to prediction?
top.models.all$ratio <- 1-(top.models.all$mrd.train/top.models.all$mrd.test)
top.models.excl$ratio <- 1-(top.models.excl$mrd.train/top.models.excl$mrd.test)

# What is the total number of model calibrations over k-folds?
n.models <- sapply(output, simplify = TRUE, function(x){
  # Sum model runs over number of parameters and k-folds (k-fold taxa have models to run)
  n.models <- nrow(x$fold1) + nrow(x$fold2) + nrow(x$fold3)
  return(n.models)
})
sum(n.models)

# Fit top models ####
# Get parameters and deviance of top models
d.mean <- arrange(d.mean, mrd.test)
models <- d.mean$Model[1:10]

write.csv(d.mean[1:100, ], 'top_100_models_invf.csv', row.names=T)

top.models.beta <- data.table()
top.models.deviance <- data.table()
for (m in 1:length(models)){
  model <- strsplit(models[m], " ")[[1]]
  model.inputs <- pre.selection[, c("SiteId", "SampId", model),with=F]

  results <- RunModel(bdm.sample, model.inputs, 100, trace=FALSE)
  results$parameters$Model <- paste(results$inf.fact, collapse=" ")
  results$deviance$Model <- paste(results$inf.fact, collapse=" ")

  top.models.beta <- bind_rows(top.models.beta, results$parameters)
  top.models.deviance <- bind_rows(top.models.deviance, results$deviance)
}

# top models without substrates ####
top.models.SS <- filter(top.models, !grepl("SS", Model))
which(top.models$Model==top.models.SS$Model[1]) # 6 parameter model ranked 62nd in predictive performance

# z-scales for top models ####
# For m in M models, k in K influence factors, for j in J taxa
M <- lapply(models, function(m){
  variables <- strsplit(m, " ")[[1]]
  # Run the model and extract results
  model.inputs <- pre.selection[, c("SiteId", "SampId", variables), with=F]
  model.results <- RunModel(bdm.sample, model.inputs, 100, trace=FALSE)
  model.results$parameters$Model <- paste(model.results$inf.fact, collapse=" ")
  setDT(model.results$parameters)
  
  # Loop through K variables for J taxa, calculate Z
  variables <- variables[!variables %in% "Temp2"]
  
  # Pass untransformed variables
  x <- p[, c("SiteId", "SampId", variables)]
  x <- x[x$SiteId %in% bdm$SiteId, ]
  setDT(x)
  
  K <- lapply(variables, function(k){
    if (k == "Temp"){
      J <- lapply(names(model.results$taxa), function(j){
        beta.jk <- model.results$parameters[Model == m & Variable %in% c("Temp", "Temp2") & Taxon == j,]
        beta.jk <- spread(beta.jk, Variable, Parameter)
        beta.j.temp1 <- beta.jk$Temp
        beta.j.temp2 <- beta.jk$Temp2
        x.k <- select_(x, .dots = k)
        x.k <- x.k[[1]]
        
        z.jk <- beta.j.temp1*(x.k-mean(x.k, na.rm=TRUE)) + beta.j.temp2*(x.k-mean(x.k, na.rm=TRUE))^2
        z.jk <- data.table(z = z.jk, x = x.k, Parameter = beta.j.temp1, Taxon = j, Variable = k, Model = m)
      })
    } else{
      J <- lapply(names(model.results$taxa), function(j){
        beta.jk <- model.results$parameters[Model == m & Variable == k & Taxon == j,]
        beta.jk <- beta.jk$Parameter
        x.k <- select_(x, .dots = k)
        x.k <- x.k[[1]]

        x.k.min <- min(x.k, na.rm=TRUE)
        x.k.mean <- mean(x.k, na.rm=TRUE)
        x.k.max <- max(x.k, na.rm=TRUE)
        # Calculate delta-z at minimum and maximum of influence factor (i.e., two points to plot slope)
        z.jk <- c(beta.jk*(x.k.min - x.k.mean),
                  beta.jk*(x.k.max - x.k.mean))
        
        # z.jk <- beta.jk*(x.k - mean(x.k, na.rm=TRUE))

        z.jk <- data.table(z = z.jk, x = c(x.k.min, x.k.max), Parameter = beta.jk, Taxon = j, Variable = k, Model = m)
      })
    }
    J <- rbindlist(J)
  })
  K <- rbindlist(K)
  
  return(K)
  cat("Model completed: ", m, "\n")
})
slopes <- rbindlist(M)
rm(M)

vimp <- data.table()
for (m in 1:length(models[1:10])){
  model <- models[m]
  variables <- strsplit(model, " ")[[1]]
  
  # Run the model and extract results
  model.inputs <- pre.selection[, c("SiteId", "SampId", variables), with=F]
  model.results <- RunModel(bdm.sample, model.inputs, 100, trace=FALSE)
  model.results$parameters$Model <- paste(model.results$inf.fact, collapse=" ")
  setDT(model.results$parameters)
  model.results$parameters$n <- n[model.results$parameters$Taxon]
  model.results$parameters <- model.results$parameters[n > 56, ]
  
  # Loop through variables
  variables <- variables[!variables %in% c("Temp", "Temp2")]
  for (k in 1:length(variables)){
    variable <- variables[k]
    beta.k <- model.results$parameters[Variable == variable]$Parameter
    
    # 90th quantile of absolute parameter estimate
    beta.k <- quantile(abs(beta.k), probs = 0.9)
    xrange <- range(model.results$mdata[, variable])
    z <- beta.k*(xrange[2]-xrange[1])
    
    mv.result <- data.table(Model = model, Variable = variable, z = z)
    vimp <- bind_rows(vimp, mv.result)
  }
}


g <- ggplot(vimp, aes(x=Model, y=z))
g <- g + geom_col()
g <- g + facet_grid(. ~ Variable, scales="free_y")
# g <- g + theme_bw(base_size = 16)
g <- g + theme(axis.text.x=element_text(angle=90,hjust=1))
g <- g + coord_cartesian(ylim=c(0,6.5))
pdf('outputs/variable_selection_bdms/top10_models_Z.pdf', width = 17, height=9)
print(g)
dev.off()

# Include occurrence frequency
slopes$n=n[slopes$Taxon]

pdf('iSDMs delta-Z [10].pdf', onefile = TRUE, paper = 'special', width = 13.5, height = 9.5)
for (m in 1:length(models)){
  model <- models[m]
  cat("Plot: ", model, "\n")
  variables <- strsplit(model, " ")[[1]]
  
  # Subset Z values based on the model
  slopes.plot <- slopes[Model == model & n > 56, ]
  
  # Set transparency per taxon
  # transparency <- n/max(n)
  # transparency <- rescale(transparency, to = c(0.10, 0.25))
  # taxa <- top.models.beta[Model == model, ]
  # taxa <- unique(taxa$Taxon)
  # slopes.plot$alpha <- transparency[taxa]
  
  # Subset raw input data for plotting
  x.raw <- p[, c("SiteId", "SampId","A10m", "A100m", "A1km", "A.EDO", "IAR", "TU.Dm", "TU.Cr", "LUD", "Urban", "UI", "F10m", "F100m", "F1km", "F.EDO", "bFRI", "FRI", "Temp", "FV", "SS1", "SS2", "SS3", "WV", "BM", "Morph", "WW", "HP")]
  x.raw <- x.raw[x.raw$SiteId %in% bdm$SiteId, ]
  x.raw$Temp2 <- x.raw$Temp^2
  x.raw$FV2 <- x.raw$FV^2
  x.raw$Noise <- runif(nrow(x.raw), 0, 1)
  
  # # Dataset for geom_blank()
  # xlim <- x.raw %>%
  #   select_(.dots=c("SiteId", "SampId", variables)) %>%
  #   filter(SampId %in% bdm$SampId) %>%
  #   gather(Variable, Value, -SiteId, -SampId) %>%
  #   filter(Variable != "Temp2") %>%
  #   group_by(Variable) %>%
  #   summarise(xmin = min(Value, na.rm=TRUE), xmax = max(Value, na.rm=TRUE))
  
  # xlim$Variable_f <- factor(xlim$Variable, levels=c("Temp", "FV", "FRI", "F10m", "F100m", "SS1", "SS2", "SS3", "IAR", "LUD", "Urban", "UI"))
  # levels(xlim$Variable_f) <- c(expression(paste("Temp (", degree, "C)")),
  #                                expression(paste("FV (m/s)")),
  #                                expression(paste("Forest-river intersection (%)")),
  #                                expression(paste("F10m (%)")),
  #                                expression(paste("F100m (%)")),
  #                                # 1. "mobile blocs>250mm +  coarse sediments (25-250mm)", 2. "mosses+hydrophytes+coarse organic matter" as two classes that are favorable, and 3. "sand and silt <2.5 mm + mud <0.1 mm" as probably rather unfavorable classes. 
  #                                expression(paste("Substrate: cobbles+coarse.sed (%)")),
  #                                expression(paste("Substrate: moss/hydrophytes/coarse-OM (%)")),
  #                                expression(paste("Substrate: Sand/silt/mud (%)")),
  #                                expression(paste("IAR (annual spray treatments / % cropland)")),
  #                                expression(paste("LUD (CE/km"^2,")")),
  #                                expression(paste("Urban (%)")),
  #                                expression(paste("Urban index (weighted urban)"))
  # )
  
  # Dataset for geom_rug()
  rug.plot <- x.raw %>%
    select_(.dots=c("SiteId", "SampId", variables)) %>%
    gather(Variable, Value, -SiteId, -SampId) %>%
    filter(Variable %in% variables) %>%
    filter(Variable != "Temp2")
  
  # rug.plot$Variable_f <- factor(rug.plot$Variable, levels=c("Temp", "FV", "FRI", "F10m", "F100m", "SS1", "SS2", "SS3", "IAR", "LUD", "Urban", "UI"))
  # levels(rug.plot$Variable_f) <- c(expression(paste("Temp (", degree, "C)")),
  #                              expression(paste("FV (m/s)")),
  #                              expression(paste("Forest-river intersection (%)")),
  #                              expression(paste("F10m (%)")),
  #                              expression(paste("F100m (%)")),
  #                              # 1. "mobile blocs>250mm +  coarse sediments (25-250mm)", 2. "mosses+hydrophytes+coarse organic matter" as two classes that are favorable, and 3. "sand and silt <2.5 mm + mud <0.1 mm" as probably rather unfavorable classes. 
  #                              expression(paste("Substrate: cobbles+coarse.sed (%)")),
  #                              expression(paste("Substrate: moss/hydrophytes/coarse-OM (%)")),
  #                              expression(paste("Substrate: Sand/silt/mud (%)")),
  #                              expression(paste("IAR (annual spray treatments / % cropland)")),
  #                              expression(paste("LUD (CE/km"^2,")")),
  #                              expression(paste("Urban (%)")),
  #                              expression(paste("Urban index (weighted urban)"))
  # )
  rug.plot.ss <- rug.plot[grepl("SS", rug.plot$Variable),]
  rug.plot <- rug.plot[!grepl("SS", rug.plot$Variable), ]
  rug.plot.ss$y <- -20
  
  g <- ggplot(data = slopes.plot)
  g <- g + geom_line(data = slopes.plot, aes(x=x, y=z, group = Taxon))
  # g <- g + geom_blank(data = xlim, aes(x = xmin))
  # g <- g + geom_blank(data = xlim, aes(x = xmax))
  g <- g + geom_rug(data=rug.plot, aes(x = Value), sides="b")
  g <- g + geom_jitter(data=rug.plot.ss, aes(x = Value, y = y), alpha = 0)
  g <- g + geom_rug(data=rug.plot.ss, aes(x = Value, y = y), position = "jitter")
  g <- g + theme_bw()
  g <- g + facet_wrap(~ Variable, scales = "free", labeller=label_parsed, strip.position = "bottom")
  # g <- g + theme(strip.background = element_blank(), strip.placement = "outside")
  g <- g + guides(alpha = FALSE)

  g <- g + labs(title = paste("Relative effects of variables: ", model),
                x = "",
                y = expression(paste(Delta, "z")))
  
  # pdf('iSDMs delta-Z.pdf', paper = 'special', width = 13.5, height = 9.5)
  g <- g + ylim(-20, 10)
  print(g)
  # dev.off()
}
dev.off()

# P-values ####
# Test significance of parameters for a given model
variables <- strsplit(models[10], " ")[[1]]
model <- RunModel(bdm.sample, pre.selection[, c("SiteId", "SampId", variables), with=F], iterations=100, trace=F)
setDT(model$parameters)
taxa <- model$taxa[model$taxa > 56]
p.values <- data.table()
for (j in 1:length(taxa)){
  taxon <- names(taxa)[j]
  m <- SelectModel(model, taxon)
  
  if(!m$optim){
    pv <- coef(summary(m))[, 4]
    
    beta <- model$parameters[Taxon == taxon, ]
    beta$pv <- pv[beta$Variable]
    
    p.values <- bind_rows(p.values, beta)
  }
}

test <- p.values %>%
  group_by(Variable) %>%
  summarise(sig = length(pv[pv<0.05]), insig = length(pv[pv>0.05])) %>%
  mutate(psig = sig/(sig+insig), ntaxa = sig+insig) %>%
  setDT()

# # Plot deviance distributions ####
# plotDev <- function(dir.file, n.models){
#   title <- paste("Predictive performance of potential models")
#   
#   top.models <- arrange(d.mean, mrd.test)
#   top.models <- top.models[1:n.models, ]
#   top <- top.models$Model
#   
#   # Get deviances for select models and exclude rare taxa
#   dev <- d[Model %in% top & n > 56, ]
#   
#   # Compute the skewness of deviances per model
#   k <- dev %>%
#     group_by(Model) %>%
#     summarise(skewness = skewness(rdev.test)) %>%
#     arrange(skewness)
#   k$scale <- rescale(k$skewness, to = c(0.1,0.3))
#   
#   # Get mean standardized deviance for each taxon
#   mdev <- dev %>%
#     filter(!is.infinite(rdev.test)) %>%
#     group_by(Model, Taxon) %>%
#     summarise(mrd.train = mean(rdev.train, na.rm=TRUE), mrd.test = mean(rdev.test, na.rm=TRUE)) %>%
#     ungroup() %>%
#     arrange(mrd.test)
#   
#   # Set up the plot
#   pdf(paste(dir.file, '.pdf', sep= ''), paper = 'a4r')
#   plot(numeric(0),numeric(0), 
#        main = title,
#        xlab = "Mean relative deviance for prediction per model \n (based on 3-fold cross-validation per taxon)",
#        ylab = "Density",
#        xlim = c(min(mdev$mrd.test, na.rm = TRUE), 2.5),
#        ylim = c(0, max(density(mdev$mrd.test)$y, na.rm = TRUE))
#   )
#   # Plot the mean and median of all selected models
#   abline(v=median(mdev$mrd.test), col=alpha("red", 0.5))
#   abline(v=mean(mdev$mrd.test), col=alpha("black", 0.5))
#   
#   # Loop through deviance distributions
#   for (m in 1:length(top)){
#     model <- top[m]
#     plot.data <- mdev[mdev$Model == model, ]
#     a <- k[k$Model == model, "scale"]
#     lines(density(plot.data$mrd.test), col=alpha("black", a))
#     cat("Progress: ", as.integer((m/length(top))*100), " %","\n")
#   }
#   dev.off()
# }

# plotDev("outputs/variable_selection_bdms/iSDM_mdev_test_100", 100)
# plotDev("outputs/variable_selection_bdms/iSDM_mdev_test_1000", 1000)

# # PROBABILITIES ###
# # For the top models, get the predicted probabilities based on k-fold cross-validation
# output.prob <- data.table()
# folds <- 1:3
# cores <- 4
# # cores <- 4
# # user  system elapsed 
# # 4.15    0.79   61.74
# # cores <- 8
# # user  system elapsed 
# # 5.38    1.54   51.91 
# ptm <- proc.time()
# for (k in folds){
#   fold <- folds[k]
#   # Get the k-fold calibration/prediction data
#   # Only model taxa that occur in calibration data
#   train.data <- as.data.frame(get(paste("train", k, sep = "")))
#   n.train <- apply(train.data[,!(colnames(train.data) %in% c("SiteId", "SampId"))], 2, sum, na.rm = T)
#   n.train <- sort(n, decreasing = TRUE)
#   n.train <- n.train[n.train > 0]
#   train.data <- train.data[, c("SiteId", "SampId", names(n.train))]
#   
#   test.data <- as.data.frame(get(paste("test", k, sep = "")))
#   
#   ### Parallel processing ###
#   # Execute models in parallel by taxa
#   # Taxa are parallelized to CPU cores; models per taxa are executed
#   cl <- makeCluster(cores)
#   registerDoParallel(cl)
#   pred.comm <- foreach(t = 1:length(n.train), .combine='rbind') %dopar% {
#     library(dplyr)
#     library(data.table)
#     # Get taxon
#     taxon <- names(n.train)[t]
#     
#     # Select calibration/prediction observations for taxon
#     train.taxon.obs <- train.data[, c("SiteId", "SampId", taxon)]
#     test.taxon.obs <- test.data[, c("SiteId", "SampId", taxon)]
#     
#     # For each model, calibrate and predict
#     out <- data.table()
#     for (m in 1:nrow(top.mean.models)){
#       model <- top.mean.models$Model[m]
#       variables <- strsplit(model, " ")
#       variables <- variables[[1]]
#       
#       # Select variables as needed
#       inputs <- predictors[, c("SiteId", "SampId", variables)]
#       
#       # Prepare calibration data
#       tt <- left_join(train.taxon.obs, inputs, by = c("SiteId", "SampId"))
#       tt <- na.omit(tt)
#       
#       # Prepare prediction data
#       td <- left_join(test.taxon.obs, inputs, by = c("SiteId", "SampId"))
#       td <- na.omit(td)
#       
#       # Run the model
#       f <- paste(taxon, "~", paste(variables, collapse = "+"))
#       f <- as.formula(f)
#       m <- glm(f, family = binomial(link = "logit"), data = tt, na.action = na.omit)
#       
#       # Relative deviance of prediction
#       prob <- predict.glm(m, newdata = td, type = "response") # get model predictions for test sites
#       y <- td[, taxon]
#       sites <- td[, c("SiteId", "SampId")]
#       mtable <- data.table(SiteId = sites[,c("SiteId")], SampId = sites[,c("SampId")], Taxon = taxon, Pred = prob, Obs = y, Model = model, Fold = fold)
#       out <- rbind(out, mtable)
#     }
#     out
#   }
#   output.prob <- rbind(output.prob, pred.comm)
#   stopCluster(cl)
#   rm(train.data, n.train, test.data, pred.comm)
# }
# proc.time() - ptm
# 
# # Plot: Map uncertainty ###
# # Get the mean predicted probability over the k-folds
# # mean.prob <- aggregate(Pred ~ SiteId + Taxon + Model + Obs, output.prob, mean)
# flow.models <- c("IAR LUD F.EDO Temp FV.old Morph", "IAR LUD F.EDO Temp FV1 Morph", "IAR LUD F.EDO Temp FV2 Morph")
# sd.prob <- aggregate(Pred ~ SiteId + Taxon + Obs, output.prob[output.prob$Model %in% flow.models,], sd)
# xy <- select(inv, SiteId, X, Y)
# xy <- distinct(xy, SiteId, X, Y)
# sd.prob$n <- n[sd.prob$Taxon]
# sd.prob <- arrange(sd.prob, desc(n))
# pdf('outputs/iSDM_SD_flow_maps.pdf', paper = 'special', width = 11, onefile = TRUE)
# for (t in 1:uniqueN(sd.prob$Taxon)){
#   page <- 1
#   taxon <- unique(sd.prob$Taxon)[t]
#   plot.data <- sd.prob[sd.prob$Taxon == taxon, ]
#   plot.data <- left_join(plot.data, xy, by = "SiteId")
#   plot.data$Obs <- as.factor(plot.data$Obs)
#   
#   g <- ggplot(inputs$ch, aes(POINT_X, POINT_Y))
#   g <- g + geom_path(lineend = "round")
#   g <- g + geom_point(data = plot.data, aes(X, Y, color = Obs, size = Pred), alpha = 0.35)
#   # g <- g + scale_size_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 7))
#   g <- g + labs(title = "Std. dev. of predicted probabilities for models with lowest deviance",
#                 subtitle = paste("individual SDM for", taxon,"| Page", page),
#                 x = "",
#                 y = "")
#   g <- g + scale_colour_brewer(palette = "Set1")
#   g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   g <- g + scale_y_continuous(breaks=NULL)
#   g <- g + scale_x_continuous(breaks=NULL)
#   g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
#   g <- g + labs(color = "Observation", size = "Standard deviation of \npredicted probability")
#   g <- g + theme_minimal(base_size = 15)
#   cat("Plotting taxa: ", taxon, "\n")
#   page <- page + 1
#   print(g)
# }
# dev.off()
# 
# # Megaboxplot ###
# mean.prob <- aggregate(Pred ~ SiteId + Taxon + Obs + Model, output.prob, mean)
# pdf('outputs/iSDM_probcomm_excl.pdf', height = 12, width = 20)
# for (m in 1:length(top.mean.models$Model)){
#   # Prepare data
#   model <- top.mean.models$Model[m]
#   plot.data <-  mean.prob[mean.prob$Model == model,]
#   plot.data$n <- n[plot.data$Taxon]
#   plot.data <- plot.data[plot.data$n > 10,]
#   plot.data <- arrange(plot.data, desc(n))
#   plot.data$Taxon <- factor(plot.data$Taxon, levels = names(n))
#   plot.data$Obs <- ifelse(plot.data$Obs == 1, 'Presence', 'Absence')
#   plot.data$Obs <- factor(plot.data$Obs, levels = c('Presence', 'Absence'))
#   # Plot the data
#   g <- ggplot(plot.data, aes(x = Taxon, y = Pred))
#   g <- g + geom_boxplot(aes(fill = Obs), alpha = 1, outlier.alpha = 0.25, outlier.size = 0.5)
#   # g <- g + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
#   g <- g + scale_fill_manual(name = "Observation", values = c("blue", "red"))
#   g <- g + labs(title = "Predicted probabilities versus observations by taxon for top mean iSDMs",
#                 subtitle = paste("Variables: ", model),
#                 y = "Mean probability over k-fold predictions (k = 3)")
#   g <- g + theme_bw()
#   g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#                  axis.text.x = element_text(angle = 45, hjust=1))
#   print(g)
# }
# dev.off()
# 
# # plotSDM for top models ##
# map.data <- mean.prob[mean.prob$Model == top.mean.models[1, "Model"]$Model, ] # Get the top model
# map.data$n <- n[map.data$Taxon]
# map.data <- arrange(map.data, desc(n))
# map.data$Taxon <- factor(map.data$Taxon, levels = names(n))
# map.data$Obs <- as.factor(map.data$Obs)
# # map.data$Obs <- ifelse(map.data$Obs == 1, 'blue', 'red')
# xy <- select(inv, SiteId, SampId, X, Y)
# xy <- distinct(xy, SiteId, SampId, X, Y)
# map.data <- left_join(map.data, xy, by = c("SiteId"))
# m <- unique(map.data$Model)
# page <- 1
# 
# # Map: probability ###
# pdf(paste("outputs/iSDM_topModel_maps", ".pdf", sep=''), paper = 'special', width = 11, onefile = TRUE)
# for (t in 1:length(unique(map.data$Taxon))){
#   taxon <- as.character(unique(map.data$Taxon)[t])
#   plot.data <- map.data[map.data$Taxon == taxon, ]
#   g <- ggplot(inputs$ch, aes(POINT_X, POINT_Y))
#   g <- g + geom_path(lineend = "round")
#   g <- g + geom_point(data = plot.data, aes(X, Y, color = Obs, size = Pred), alpha = 0.35)
#   # g <- g + scale_size_continuous(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 7))
#   g <- g + scale_radius(limits = c(0,1), breaks = seq(0, 1, 0.2), range = c(2, 7))
#   g <- g + labs(title = paste("Mean predictions (3-fold CV) vs observations of",paste(taxon)),
#                 subtitle = paste("individual SDM -",taxon,"~", paste(m, collapse = "+", sep = " "), "- page", page),
#                 x = "",
#                 y = "")
#   g <- g + scale_colour_brewer(palette = "Set1")
#   g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
#   g <- g + scale_y_continuous(breaks=NULL)
#   g <- g + scale_x_continuous(breaks=NULL)
#   g <- g + guides(colour = guide_legend(override.aes = list(size=6)))
#   g <- g + labs(color = "Observation", size = "Predicted probability")
#   # g <- g + labs(names = c("1" = "Present", "0" = "Absent"))
#   g <- g + theme_minimal(base_size = 15)
#   cat("Plotting taxa: ", taxon, "\n")
#   page <- page + 1
#   print(g)
# }
# dev.off()
# 
# 
# ### OVAT ###
# # Loop through d.mean.excl:
# # drop OVAT and do k-fold CV
# # get mean relative deviance for train/test
# top.models.to.test <- arrange(d.mean.excl, mrd.test)
# top.models.to.test <- top.models.to.test[1:5,]
# cores <- 4
# output.ovat <- data.table()
# # Loop through models
# for (m in 1:nrow(top.models.to.test)){
#   model <- top.models.to.test[m, "Model"]
#   variables <- strsplit(model, " ")
#   variables <- variables[[1]]
#   
#   folds <- 1:3
#   # Do k-fold CV on full model
#   for (k in 1:length(folds)){
#     # Get the k-fold calibration/prediction data
#     # Only model taxa that occur in calibration data
#     train.data <- as.data.frame(get(paste("train", k, sep = "")))
#     n.train <- apply(train.data[,!(colnames(train.data) %in% c("SiteId", "SampId"))], 2, sum, na.rm = T)
#     n.train <- sort(n, decreasing = TRUE)
#     n.train <- n.train[n.train > 0]
#     train.data <- train.data[, c("SiteId", names(n.train))]
#     
#     test.data <- as.data.frame(get(paste("test", k, sep = "")))
#     
#     # For each k-fold, fit sub-models to taxa in parallel
#     cl <- makeCluster(cores)
#     registerDoParallel(cl)
#     dev.comm <- foreach(t = 1:length(n.train), .combine='rbind') %dopar% {
#       library(dplyr)
#       # Get taxon
#       taxon <- names(n.train)[t]
#       
#       # Select calibration/prediction observations for taxon
#       train.taxon.obs <- train.data[, c("SiteId", taxon)]
#       test.taxon.obs <- test.data[, c("SiteId", taxon)]
#       # Select variables as needed
#       inputs <- i[, c("SiteId", variables)]
#       
#       # Prepare calibration data
#       tt <- left_join(train.taxon.obs, inputs, by = "SiteId")
#       tt <- na.omit(tt)
#       
#       # Prepare prediction data
#       td <- left_join(test.taxon.obs, inputs, by = "SiteId")
#       td <- na.omit(td)
#       
#       # Run the model
#       f <- paste(taxon, "~", paste(variables, collapse = "+"))
#       f <- as.formula(f)
#       m <- glm(f, family = binomial(link = "logit"), data = tt, na.action = na.omit)
#       
#       # Relative deviance of calibration
#       relative.dev.train <- m$deviance/nrow(tt)
#       
#       # Relative deviance of prediction
#       prob <- predict.glm(m, newdata = td, type = "response") # get model predictions for test sites
#       y <- td[, taxon]
#       dev.resid <- sign(y-0.5)*sqrt(-2*(y*log(prob)+(1-y)*log(1-prob)))
#       r.dev <- sum(dev.resid^2, na.rm = TRUE)
#       relative.dev.test <- r.dev/nrow(td)
#       rm(prob, y, dev.resid, r.dev)
#       
#       # Relative deviance statistics
#       dev.taxon <- c(relative.dev.train, relative.dev.test)
#       
#       # Prepare output
#       dev.taxon <- data.frame(t(dev.taxon), stringsAsFactors = FALSE)
#       colnames(dev.taxon) <- c("rdev.train", "rdev.test")
#       dev.taxon <- cbind(dev.taxon, taxon)
#     }
#     stopCluster(cl)
#     
#     # Tidy data and output
#     dev.comm$taxon <- as.character(dev.comm$taxon)
#     dev.comm$Fold <- k
#     dev.comm$Variable <- "All variables"
#     dev.comm$Model <- model
#     output.ovat <- rbind(output.ovat, dev.comm)
#   }
#   
#   # Loop through sub-models
#   for (v in 1:length(variables)){
#     sub.model.inputs <- variables[-v]
#     
#     folds <- 1:3
#     # For each sub-model, do k-fold cross-validation
#     for (k in 1:length(folds)){
#       # Get the k-fold calibration/prediction data
#       # Only model taxa that occur in calibration data
#       train.data <- as.data.frame(get(paste("train", k, sep = "")))
#       n.train <- apply(train.data[, !(colnames(train.data) %in% c("SiteId", "SampId"))], 2, sum, na.rm = T)
#       n.train <- sort(n, decreasing = TRUE)
#       n.train <- n.train[n.train > 0]
#       train.data <- train.data[, c("SiteId", names(n.train))]
#       
#       test.data <- as.data.frame(get(paste("test", k, sep = "")))
#       
#       # For each k-fold, fit sub-models to taxa in parallel
#       cl <- makeCluster(cores)
#       registerDoParallel(cl)
#       dev.comm <- foreach(t = 1:length(n.train), .combine='rbind') %dopar% {
#         library(dplyr)
#         # Get taxon
#         taxon <- names(n.train)[t]
#         
#         # Select calibration/prediction observations for taxon
#         train.taxon.obs <- train.data[, c("SiteId", taxon)]
#         test.taxon.obs <- test.data[, c("SiteId", taxon)]
#         # Select variables as needed
#         inputs <- i[, c("SiteId", sub.model.inputs)]
#         
#         # Prepare calibration data
#         tt <- left_join(train.taxon.obs, inputs, by = "SiteId")
#         tt <- na.omit(tt)
#         
#         # Prepare prediction data
#         td <- left_join(test.taxon.obs, inputs, by = "SiteId")
#         td <- na.omit(td)
#         
#         # Run the model
#         f <- paste(taxon, "~", paste(sub.model.inputs, collapse = "+"))
#         f <- as.formula(f)
#         m <- glm(f, family = binomial(link = "logit"), data = tt, na.action = na.omit)
#         
#         # Relative deviance of calibration
#         relative.dev.train <- m$deviance/nrow(tt)
#         
#         # Relative deviance of prediction
#         prob <- predict.glm(m, newdata = td, type = "response") # get model predictions for test sites
#         y <- td[, taxon]
#         dev.resid <- sign(y-0.5)*sqrt(-2*(y*log(prob)+(1-y)*log(1-prob)))
#         r.dev <- sum(dev.resid^2, na.rm = TRUE)
#         relative.dev.test <- r.dev/nrow(td)
#         rm(prob, y, dev.resid, r.dev)
#         
#         # Relative deviance statistics
#         dev.taxon <- c(relative.dev.train, relative.dev.test)
#         
#         # Prepare output
#         dev.taxon <- data.frame(t(dev.taxon), stringsAsFactors = FALSE)
#         colnames(dev.taxon) <- c("rdev.train", "rdev.test")
#         dev.taxon <- cbind(dev.taxon, taxon)
#       }
#       stopCluster(cl)
#       
#       # Tidy data and output
#       dev.comm$taxon <- as.character(dev.comm$taxon)
#       dev.comm$Fold <- k
#       dev.comm$Variable <- variables[v]
#       dev.comm$Model <- model
#       output.ovat <- rbind(output.ovat, dev.comm)
#     }
#   }
#   rm(model, variables, sub.model.inputs, folds, n.train, train.data, test.data, dev.comm, v, k)
#   cat("OVAT k-fold cross-validation completed for model:", m, "\n")
# }
# 
# # For each model and variable dropped, get the mean relative deviance for training/testing
# ovat.mrd.train <- aggregate(rdev.train ~ taxon + Model + Variable, output.ovat, mean)
# ovat.mrd.test <- aggregate(rdev.test ~ taxon + Model + Variable, output.ovat, mean)
# 
# ovat.data <- left_join(ovat.mrd.test, ovat.mrd.train, by = c("taxon", "Model", "Variable"))
# ovat.data$n <- n[ovat.data$taxon]
# ovat.data <- ovat.data[ovat.data$n > 10, ]
# rm(ovat.mrd.train, ovat.mrd.test)
# 
# ### Plot: mRD OVAT k-fold CV ###
# pdf('outputs/GLM_OVAT_mean.pdf', paper = 'a4r')
# for (m in 1:nrow(top.models.to.test)){
#   model <- top.models.to.test[m, "Model"]
#   plot.data <- filter(ovat.data, Model == model)
#   
#   g <- ggplot(plot.data, aes(x = Variable, y = rdev.test, fill = Variable))
#   g <- g + geom_boxplot()
#   # g <- g + facet_grid(.~Fold)
#   g <- g + theme_bw(base_size = 15)
#   # g <- g + theme(strip.background=element_rect(fill="black"))
#   # g <- g + theme(strip.text=element_text(color="white", face="bold"))
#   g <- g + labs(title = "Relative deviance of prediction for iSDM",
#                 subtitle = "Top five predictive models: drop OVAT and 3-fold cross-validation",
#                 x = "Variable dropped from model",
#                 y = "Relative deviance for prediction")
#   g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#                  axis.text.x = element_text(angle = 90))
#   g <- g + scale_fill_brewer(palette = "Set1")
#   g <- g + scale_y_continuous(limits = c(0, 2))
#   print(g)
# }
# dev.off()
# 
# ### Plot: OVAT k-fold CV by fold ###
# output.ovat$n <- n[output.ovat$taxon]
# 
# t <- output.ovat
# f <- output.ovat[output.ovat$Variable == 'All variables', c("rdev.test", "taxon", "Fold", "Variable", "Model")]
# colnames(f)[1] <- c("full.rdev.test")
# t <- left_join(t, f, by = c("taxon", "Fold", "Model"))
# t$difference <- t$full.rdev.test - t$rdev.test
# t <- t[t$Variable.x != 'All variables',]
# pdf('outputs/GLM_OVAT_top_models.pdf', paper = 'a4r')
# for (m in 1:nrow(top.models.to.test)){
#   model <- top.models.to.test[m, "Model"]
#   plot.data <- filter(t, Model == model & n > 10)
#   
#   g <- ggplot(plot.data, aes(x = Variable.x, y = difference))
#   g <- g + geom_boxplot(alpha = 0.5)
#   g <- g + facet_grid(Fold ~ ., scales = "free")
#   g <- g + theme_bw(base_size = 15)
#   
#   g <- g + labs(title = "Deviance of full model minus OVAT drop (n>10)",
#                 subtitle = "Top five predictive models: remove OVAT with 3-fold cross-validation",
#                 x = "Variable dropped from model",
#                 y = "Relative deviance (full model - OVAT)")
#   # g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#   #                axis.text.x = element_text(angle = 90))
#   g <- g + theme(axis.text.x = element_text(angle = 90))
#   g <- g + scale_fill_brewer(palette = "Set1")
#   g <- g + scale_y_continuous(limits = c(-0.25, 0.25))
#   print(g)
# }
# dev.off()
# 
# ### Fit top models to full dataset ###
# top.model.fit <- data.table()
# for (m in 1:nrow(top.models.to.test)){
#   model <- top.models.to.test[m, "Model"]
#   variables <- strsplit(model, " ")
#   variables <- variables[[1]]
#   
#   # Fit the full model
#   dev.comm <- sapply(names(n), function(taxon){
#     # Prepare the data
#     mdata <- bdm.all[, c("SiteId", taxon)]
#     inputs <- i[, c("SiteId", variables)]
#     mdata <- left_join(mdata, inputs, by = "SiteId")
#     
#     # Run the model
#     f <- paste(taxon, "~", paste(variables, collapse = "+"))
#     f <- as.formula(f)
#     m <- glm(f, family = binomial(link = "logit"), data = mdata, na.action = na.omit)
#     D2 <- (m$null.deviance - m$deviance)/m$null.deviance
#     rm(mdata, inputs, f, m)
#     D2
#   })
#   dev.comm <- data.frame(Species = names(n), D2 = dev.comm, Model = model, Variable = "Full model", stringsAsFactors = FALSE, row.names = 1:length(dev.comm))
#   top.model.fit <- rbind(top.model.fit, dev.comm)
#   
#   # OVAT fit with full dataset
#   cl <- makeCluster(cores)
#   registerDoParallel(cl)
#   dev.comm <- foreach(v = 1:length(variables), .combine = 'rbind') %dopar% {
#     sub.model.inputs <- variables[-v]
#     library(dplyr)
#     dev = numeric(0)
#     for (t in 1:length(n)){
#       # Get taxon
#       taxon <- names(n)[t]
#       
#       # Prepare the data
#       mdata <- bdm.all[, c("SiteId", taxon)]
#       inputs <- i[, c("SiteId", sub.model.inputs)]
#       mdata <- left_join(mdata, inputs, by = "SiteId")
#       
#       # Run the model
#       f <- paste(taxon, "~", paste(sub.model.inputs, collapse = "+"))
#       f <- as.formula(f)
#       m <- glm(f, family = binomial(link = "logit"), data = mdata, na.action = na.omit)
#       rm(mdata, inputs, f)
#       dev[t] <- (m$null.deviance - m$deviance)/m$null.deviance
#     }
#     u <- data.frame(Species = names(n), 
#                     D2 = dev, 
#                     Variable = variables[v],
#                     Model = paste(variables, sep = " ", collapse = " "), 
#                     stringsAsFactors = FALSE,
#                     row.names = 1:length(n))
#   }
#   stopCluster(cl)
#   top.model.fit <- rbind(top.model.fit, dev.comm)
# }
# 
# sum(is.infinite(top.model.fit$D2))
# 
# ### Plot: top models with full data ###
# top.model.fit$n <- n[top.model.fit$Species]
# top.model.fit.excl <- top.model.fit[top.model.fit$n > 10, ]
# top.model.fit.excl$D2[is.infinite(top.model.fit.excl$D2)] <- 0
# h <- aggregate(D2 ~ Model + Variable, top.model.fit.excl, mean) # tidy data for horizontal line
# 
# pdf('outputs/GLM_OVAT_excl.pdf', paper = 'a4r')
# for (m in 1:uniqueN(top.model.fit$Model)){
#   model <- unique(top.model.fit$Model)[m]
#   plot.data <- filter(top.model.fit.excl, Model == model)
#   h.data <- h[h$Model == model & h$Variable == 'Full model',]
#   g <- ggplot(plot.data, aes(x = Variable, y = D2, fill = Variable))
#   g <- g + geom_boxplot()
#   g <- g + geom_hline(data=h.data,aes(yintercept=D2))
#   g <- g + theme_bw(base_size = 15)
#   g <- g + labs(title = model,
#                 subtitle = expression(paste("Top five predictive models: drop OVAT and obtain D"^2," for entire BDM")),
#                 x = "Variable dropped from model",
#                 y = expression(paste("D"^2)))
#   g <- g + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
#                  axis.text.x = element_text(angle = 90))
#   g <- g + scale_fill_brewer(palette = "Set1")
#   g <- g + scale_y_continuous(limits = c(0, 1))
#   print(g)
# }
# dev.off()
# 
# ### plot mean relative deviance for testing over all top models
# tmtt <- unique(top.models.to.test[1:20, "Model"])
# pdf('../PredictionError_density.pdf', paper = 'special', width = 11, onefile = TRUE)
# for (m in 1:length(tmtt)){
#   model <- tmtt[m]
#   plot.data <- top.models.excl[top.models.excl$Model %in% tmtt, ]
#   xmin <- min(plot.data$mrd.test)
#   xmax <- max(plot.data$mrd.test)
#   
#   # set up a plot for influence factor
#   plot(numeric(0),numeric(0), 
#        main = paste("jSDM posterior parameter distributions for", kname),
#        xlab = expression(paste(beta[kj]^{taxa})),
#        ylab = "Density"
#   )
#   lines(density(beta.taxa[, k, j]), col=pc)
# }








