# Model-Based Exploratory Analyses of Invertebrate Data
# Inference Model: logistic, hierarchical JSDM with/without site-specific random effects
# Note: dplyr package function left_join() is used instead of base R match() or merge() because:
# - merge() doesn't maintain row order of original datasets; see ?merge:
# "The rows are by default lexicographically sorted on the common columns, but for sort = FALSE are in an unspecified order" ???
# - match() matches based on the first match, cannot be used for duplicates
# - dplyr *_join functions are similar to merge() but do maintain row order!

# dplyr() functions used:
# - left_join(x,y,by="") # merge x and y by column keeping all rows in left dataset
# - select() # select columns

# Global parameters ####
community       <- "community.csv"
predictors      <- "predictors.csv"

comm.corr       <- T
site.effects    <- T
n.latent        <- 2
lat.site        <- F

generate.res    <- F

crossvalid      <- FALSE
dir.modinp      <- "inputs"

sampsize        <- 10000
thin            <- 5
n.chain         <- 1
prob.defpri     <- 0.02
thresh.sig      <- 1
fact.sd         <- 1

# Packages ####
if ( !require("coda") )  { install.packages("coda");  library("coda") }
if ( !require("rstan") ) { install.packages("rstan"); library("rstan") }
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") }
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Read data ####
# read occurrence data:
occur.taxa <- read.table(paste(dir.modinp, community, sep="/"), header=TRUE, sep=",",
                         stringsAsFactors=FALSE)

# read environmental conditions:
env.cond <- read.table(paste(dir.modinp, predictors, sep="/"), header=TRUE, sep=",",
                       stringsAsFactors=FALSE)

inf.fact <- colnames(select(env.cond, -SiteId, -SampId))

# join the environmental conditions to the occurrence data
env.cond <- left_join(occur.taxa[, c("SiteId", "SampId")], env.cond, by = c("SiteId", "SampId"))
 
# Pre-process data ####
# drop rows with incomplete influence factors:
ind <- !apply(is.na(env.cond[,inf.fact]),1,FUN=any)
ind <- ifelse(is.na(ind),FALSE,ind)
occur.taxa <- occur.taxa[ind, ]
env.cond   <- env.cond[ind, c("SiteId", "SampId", inf.fact)]
print(paste(sum(!ind),"sites/samples excluded because of incomplete influence factors"))

sites <- occur.taxa$SiteId
samples <- occur.taxa$SampId

n.sites <- length(unique(sites))
n.samples <- length(samples)

occur.taxa$SiteId <- NULL
occur.taxa$SampId <- NULL

# drop TAXA without observations at the selected sites:
ind <- apply(occur.taxa,2,sum, na.rm = TRUE) > 0
occur.taxa <- occur.taxa[, ind]
n.taxa <- ncol(occur.taxa)

unique.sites <- unique(sites)
siteIND <- match(sites, unique.sites)

# replace NA with -1
occur.taxa.na.encod <- occur.taxa
for(j in 1:ncol(occur.taxa.na.encod)) occur.taxa.na.encod[,j] <- ifelse(is.na(occur.taxa.na.encod[,j]), -1, occur.taxa.na.encod[,j])

# choose data for identification and validation:
ind <- 1:nrow(occur.taxa)
set.ident <- rep(TRUE,length(ind))
set.valid <- set.ident
if (crossvalid) {
  set.ident <- ifelse(floor(ind/2)*2+1==ind,TRUE,FALSE)
  set.valid <- !set.ident
  print("cross validation:")
  print(paste("Number of sites chosen for identification:",sum(set.ident)))
  print(paste("Number of sites chosen for validation:    ",sum(set.valid)))
} else {
  print("no cross validation:")
  print(paste("using all",nrow(occur.taxa),"sites for identification and validation"))
}

# Run GLMs ####
# GLM analysis of occurrence data (needed for initial conditions of Markov chains):

res.glm <- data.frame(matrix(NA,nrow=ncol(occur.taxa),ncol=2+1+length(inf.fact),
                             dimnames=list(colnames(occur.taxa),c("present","absent","alpha",inf.fact))))
for ( j in 1:ncol(occur.taxa) )
{
  taxon <- colnames(occur.taxa)[j]
  n <- sum(occur.taxa[set.ident,j], na.rm = TRUE) # loop through rows and taxa
  res.glm[j,1] <- n
  res.glm[j,2] <- sum(set.ident)-n
  if ( n > 5 )
  {
    occur <- occur.taxa[set.ident,j]
    res.glm.j <- glm(formula(paste("occur ~",paste(inf.fact,collapse="+"))),
                     family="binomial",
                     data=env.cond[set.ident,])
    est <- summary(res.glm.j)$coefficients[,"Estimate"]
    sd  <- summary(res.glm.j)$coefficients[,"Std. Error"]
    est <- ifelse(abs(est)>thresh.sig*sd,est,0)
    res.glm[j,3:ncol(res.glm)] <- as.numeric(est)
    res.glm[j,3] <- as.numeric(summary(res.glm.j)$coefficients[,"Estimate"][1])
  }
  else
  {
    if ( n > 0 )
    {
      occur <- occur.taxa[set.ident,j]
      res.glm.j <- glm(formula(paste("occur ~","1")),
                       family="binomial",
                       data=env.cond[set.ident,])
      res.glm[j,3] <- as.numeric(summary(res.glm.j)$coefficients[,"Estimate"][1])
      res.glm[j,4:ncol(res.glm)] <- 0
    }
    else
    {
      res.glm[j,3] <- -Inf
      res.glm[j,4:ncol(res.glm)] <- 0
    }
  }
}

# Model ####

submodel <- ifelse(comm.corr,"corr","nocorr")
submodel <- paste(submodel,ifelse(site.effects,"site","nosite"),sep="_")
submodel <- paste(submodel,ifelse(n.latent>0,paste(n.latent,"latvar",sep=""),"nolatvar"),sep="_")
if ( n.latent > 0 ) submodel <- paste(submodel,ifelse(lat.site,"site","samp"),sep="")
file.model <- paste("Inv_JSDM_",submodel,".stan",sep="")

### Define priors ####

logit <- function(p) { return(log(p/(1-p))) }
z.range.max <- logit(0.95)-logit(0.05)

mu.alpha.comm.pripar    <- 0
sigma.alpha.comm.pripar <- length(inf.fact)*z.range.max/4
sigma.beta.comm.pripar <- z.range.max/(0.5*apply(env.cond[inf.fact],2,function(k){
  max(k)-min(k)
}))
mu.beta.comm.pripar <- rep(0,length(sigma.beta.comm.pripar))
names(mu.beta.comm.pripar) <- names(sigma.beta.comm.pripar)

# compilation of input data:
data <- list(n_sites                 = n.sites,
             n_samples               = n.samples,
             n_taxa                  = n.taxa,
             n_pred                  = length(inf.fact),
             mu_alpha_comm_pripar    = mu.alpha.comm.pripar,
             sigma_alpha_comm_pripar = sigma.alpha.comm.pripar,
             mu_beta_comm_pripar     = mu.beta.comm.pripar,
             sigma_beta_comm_pripar  = sigma.beta.comm.pripar,
             fact_sd                 = fact.sd,
             x                       = as.matrix(env.cond[set.ident,inf.fact]),
             y                       = as.matrix(occur.taxa.na.encod[set.ident,]))

if ( site.effects ) {
  data$sigma_gamma_pripar <- 2*logit(1-prob.defpri)/10
  data$siteIND <- siteIND
}
if ( n.latent > 0 )
{
  data$sigma_beta_lat <- z.range.max/4
  if ( n.latent > 1 )
  {
    data$n_latent     <- n.latent
  }
}

### Initialize chains ####
# definition of (critical) starting points of Markov chains:
init <- list()
for ( i in 1:n.chain )
{
  init[[i]] <- list()
  init.alpha <- res.glm[,"alpha"]
  init.beta  <- t(res.glm[,inf.fact])
  # cut extreme values:
  min.max <- quantile(init.alpha[init.alpha!=-Inf],probs=c(0.2,0.8))
  init.alpha <- ifelse(init.alpha<min.max[1],min.max[1],ifelse(init.alpha>min.max[2],min.max[2],init.alpha))
  init.alpha <- init.alpha*(1+0.1*rnorm(length(init.alpha))) 
  init[[i]][["alpha_taxa"]] <- init.alpha
  for ( k in 1:length(inf.fact) )
  {
    min.max <- quantile(init.beta[k,],probs=c(0.2,0.8))
    init.beta[k,] <- ifelse(init.beta[k,]<min.max[1],min.max[1],ifelse(init.beta[k,]>min.max[2],min.max[2],init.beta[k,]))
    init.beta[k,] <- init.beta[k,]*(1+0.1*rnorm(ncol(init.beta)))
  }
  init[[i]][["beta_taxa"]] <- init.beta
  init[[i]][["mu_alpha_comm"]] <- mean(init[[i]]$alpha_taxa)
  init[[i]][["mu_beta_comm"]]  <- apply(init[[i]]$beta_taxa,1,mean)
  s.pripar = sqrt(log(1+fact.sd^2))
  init[[i]][["sigma_alpha_comm"]] <- sigma.alpha.comm.pripar
  init[[i]][["sigma_beta_comm"]]  <- sigma.beta.comm.pripar

  if ( site.effects )
  {
    init[[i]][["sigma_gamma"]] <- data$sigma_gamma_pripar
    init[[i]][["gamma_site"]] <- rep(0,n.sites)
  }
}


### Model definition ####

# data
# ----

model.code <- paste("// Model:    Logistic, hierarchical JSDM.\n",
                    "// Submodel: ",submodel,"\n",
                    "// -----------------------------------------------------------------------\n\n",
                    "// data:\n\n",
                    "data {\n",
                    "  int                              n_taxa;\n",
                    "  int                              n_sites;\n",
                    "  int                              n_samples;\n",
                    "  int                              n_pred;\n",
                    "  real                             mu_alpha_comm_pripar;\n",
                    "  real                             sigma_alpha_comm_pripar;\n",
                    "  vector[n_pred]                   mu_beta_comm_pripar;\n",
                    "  vector[n_pred]                   sigma_beta_comm_pripar;\n",
                    "  real                             fact_sd;\n",
                    sep = "")
if ( site.effects )
{
  model.code <- paste(model.code,
                      "  real                             sigma_gamma_pripar;\n",
                      "  int                              siteIND[n_samples];\n",
                      sep="")
}
if ( n.latent > 0 )
{
  model.code <- paste(model.code,
                      "  real                             sigma_beta_lat;\n",
                      sep="")
  if ( n.latent > 1 )
  {
    model.code <- paste(model.code,
                        "  int                            n_latent;\n",
                        sep="")
  }
}                      
model.code <- paste(model.code,
                    "  matrix[n_samples,n_pred]         x;\n",
                    "  int<lower=-1,upper=1>            y[n_samples,n_taxa];\n",
                    "}\n\n", 
                    sep="")

# parameters
# ----------

model.code <- paste(model.code,
                    "// parameters:\n\n",
                    "parameters {\n",
                    "  real                             mu_alpha_comm;\n",
                    "  real<lower=0>                    sigma_alpha_comm;\n",
                    "  vector[n_pred]                   mu_beta_comm;\n",
                    "  vector<lower=0>[n_pred]          sigma_beta_comm;\n",
                    "  vector[n_taxa]                   alpha_taxa;\n",
                    "  matrix[n_pred,n_taxa]            beta_taxa;\n",
                    sep="")
if (comm.corr)
{
  model.code <- paste(model.code,
                      "//corr_matrix[n_pred]              corr_comm;\n",
                      "  cholesky_factor_corr[n_pred]     corrfact_comm;\n",
                      sep="")
}
if ( site.effects )
{
  model.code <- paste(model.code,
                      "  real<lower=0>                    sigma_gamma;\n",
                      "  vector[n_sites]                  gamma_site;\n",
                      sep="")
}
if ( n.latent == 1 )
{
  if ( lat.site )
  {
    model.code <- paste(model.code,
                        "  vector[n_sites]                  x_lat;\n",
                        sep="")
  }
  else
  {
    model.code <- paste(model.code,
                        "  vector[n_samples]                x_lat;\n",
                        sep="")
  }
  model.code <- paste(model.code,
                      "  vector[n_taxa]                   beta_lat;\n",
                      sep="")
}  
if ( n.latent > 1 )
{
  if ( lat.site )
  {
    model.code <- paste(model.code,
                        "  matrix[n_sites,n_latent]         x_lat;\n",
                        sep="")
  }
  else
  {
    model.code <- paste(model.code,
                        "  matrix[n_samples,n_latent]       x_lat;\n",
                        sep="")
  }
  model.code <- paste(model.code,
                      "  matrix[n_latent,n_taxa]          beta_lat;\n",
                      sep="")
}
model.code <- paste(model.code,
                    "}\n\n", 
                    sep="")

# local variables
# ---------------

model.code <- paste(model.code,
                    "// model definition:\n\n",
                    "model {\n\n",
                    "  // local variables:\n\n",
                    "  real z[n_samples,n_taxa];\n",
                    "  real p[n_samples,n_taxa];\n",
                    "  real s_pripar;\n\n",
                    sep="")

# root nodes
# ----------
  
model.code <- paste(model.code,
                    "  // root nodes:\n\n",
                    "  mu_alpha_comm ~ normal(mu_alpha_comm_pripar,sigma_alpha_comm_pripar);\n",
                    "  s_pripar = sqrt(log(1+fact_sd^2));\n",
                    "  sigma_alpha_comm   ~ lognormal(log(sigma_alpha_comm_pripar)-0.5*s_pripar^2,s_pripar);\n",
                    "  for ( k in 1:n_pred ) {\n",
                    "    mu_beta_comm[k] ~ normal(mu_beta_comm_pripar[k],sigma_beta_comm_pripar[k]);\n",
                    "    s_pripar = sqrt(log(1+fact_sd^2));\n",
                    "    sigma_beta_comm[k]   ~ lognormal(log(sigma_beta_comm_pripar[k])-0.5*s_pripar^2,s_pripar);\n",
                    "  }\n",
                    sep="")
if ( site.effects )
{
  model.code <- paste(model.code,
                      "  s_pripar = sqrt(log(1+fact_sd^2));\n",
                      "  sigma_gamma ~ lognormal(log(sigma_gamma_pripar)-0.5*s_pripar^2,s_pripar);\n",
                      sep="")
}
if ( comm.corr )
{
  model.code <- paste(model.code,
                      "  corrfact_comm ~ lkj_corr_cholesky(2);\n",
                      sep="")
}
if ( n.latent == 1 )
{
  if ( lat.site )
  {
    model.code <- paste(model.code,
                        "  for ( i in 1:n_sites ) {\n",
                        sep="")
  }
  else
  {
    model.code <- paste(model.code,
                        "  for ( i in 1:n_samples ) {\n",
                        sep="")
  }
  model.code <- paste(model.code,
                      "    x_lat[i] ~ normal(0,1);\n",
                      "  }\n",
                      "  for ( j in 1:n_taxa ) {\n",
                      "    beta_lat[j] ~ normal(0,sigma_beta_lat);\n",
                      "  }\n",
                      sep="")
}
if ( n.latent > 1 )
{
  model.code <- paste(model.code,
                      "  for ( k in 1:n_latent ) {\n",
                      sep="")
  if ( lat.site )
  {
    model.code <- paste(model.code,
                        "    for ( i in 1:n_sites ) {\n",
                        sep="")
  }
  else
  {
    model.code <- paste(model.code,
                        "    for ( i in 1:n_samples ) {\n",
                        sep="")
  }
  model.code <- paste(model.code,
                      "      if ( k==1 && i==11)\n",
                      "        x_lat[i,k] ~ normal(0,1) T[0,];\n",
                      "      else if (k==2 && i==13)\n",
                      "        x_lat[i,k] ~ normal(0,1) T[0,];\n",
                      "      else\n",
                      "        x_lat[i,k] ~ normal(0,1);\n",
                      "    }\n",
                      "    for ( j in 1:n_taxa ) {\n",
                      "      beta_lat[k,j] ~ normal(0,sigma_beta_lat);\n",
                      "    }\n",
                      "  }\n",
                      sep="")
}
model.code <- paste(model.code,"\n",sep="")

# intermediate nodes
# ------------------

model.code <- paste(model.code,
                    "  // intermediate nodes:\n\n",
                    "  for ( j in 1:n_taxa ) {\n",
                    "    alpha_taxa[j] ~ normal(mu_alpha_comm,sigma_alpha_comm);\n",
                    sep="")
if ( !comm.corr )
{
  model.code <- paste(model.code,
                      "    for ( k in 1:n_pred ) {\n",
                      "      beta_taxa[k,j] ~ normal(mu_beta_comm[k],sigma_beta_comm[k]);\n",
                      "    }\n",
                      sep="")
} else {
  model.code <- paste(model.code,
                      "  //beta_taxa[,j] ~ multi_normal(mu_beta_comm,quad_form_diag(corr_comm,sigma_beta_comm));\n",
                      "    beta_taxa[,j] ~ multi_normal_cholesky(mu_beta_comm,diag_pre_multiply(sigma_beta_comm,corrfact_comm));\n",
                      sep="")
}
model.code <- paste(model.code,
                      "  }\n",
                    sep="")
if ( site.effects )
{
  model.code <- paste(model.code,
                      "  for ( i in 1:n_sites ) {\n",
                      "    gamma_site[i] ~ normal(0,sigma_gamma);\n",
                      "  }\n",
                      sep="")
}
model.code <- paste(model.code,
                    "  \n",
                    sep="")

# end nodes
# ---------

model.code <- paste(model.code,
                    "  // end nodes:\n\n",
                    "  for ( i in 1:n_samples ) {\n",
                    "    for ( j in 1:n_taxa ) {\n",
                    "      z[i,j] = alpha_taxa[j];\n",
                    "      for ( k in 1:n_pred ) {\n",
                    "        z[i,j] = z[i,j] + beta_taxa[k,j]*x[i,k];\n",
                    "      }\n",
                    sep="")
if ( site.effects )
{
  model.code <- paste(model.code,
                      "      z[i,j] = z[i,j] + gamma_site[siteIND[i]];\n",        
                      sep="")
}
if ( n.latent == 1 )
{
  if ( lat.site )
  {
    model.code <- paste(model.code,
                        "      z[i,j] = z[i,j] + x_lat[siteIND[i]]*beta_lat[j];\n",        
                        sep="")
  }
  else
  {
    model.code <- paste(model.code,
                        "      z[i,j] = z[i,j] + x_lat[i]*beta_lat[j];\n",        
                        sep="")
  }
}
if ( n.latent > 1 )
{
  model.code <- paste(model.code,
                      "      for ( k in 1:n_latent ) {\n",
                      "        z[i,j] = z[i,j] + x_lat[siteIND[i],k]*beta_lat[k,j];\n",
                      "      }\n",
                      sep="")
}
model.code <- paste(model.code,
                    "      p[i,j] = 1/(1+exp(-z[i,j]));\n",
                    "      if(y[i,j]>=0) {\n",
                    "        y[i,j] ~ bernoulli(p[i,j]);\n",
                    "      }\n",
                    "    }\n",
                    "  }\n",
                    "}\n\n",
                    sep="")

# generated quantities:
# ---------------------

if ( generate.res )
{
  
  model.code <- paste(model.code,
                      "// generated quantities:\n\n",
                      "generated quantities {\n\n",
                      "  // output variables:\n\n",
                      "  real z_mod[n_samples,n_taxa];\n",
                      "  real p_mod[n_samples,n_taxa];\n",
                      "  real deviance[n_taxa];\n",
                      "\n",
                      sep="")
  
  model.code <- paste(model.code,
                      "  // end nodes:\n\n",
                      "  for ( i in 1:n_samples ) {\n",
                      "    for ( j in 1:n_taxa ) {\n",
                      "      z_mod[i,j] = alpha_taxa[j];\n",
                      "      for ( k in 1:n_pred ) {\n",
                      "        z_mod[i,j] = z_mod[i,j] + beta_taxa[k,j]*x[i,k];\n",
                      "      }\n",
                      sep="")
  if ( site.effects )
  {
    model.code <- paste(model.code,
                        "      z_mod[i,j] = z_mod[i,j] + gamma_site[siteIND[i]];\n",        
                        sep="")
  }
  if ( n.latent == 1 )
  {
    model.code <- paste(model.code,
                        "      z_mod[i,j] = z_mod[i,j] + x_lat[siteIND[i]]*beta_lat[j];\n",        
                        sep="")
  }
  if ( n.latent > 1 )
  {
    model.code <- paste(model.code,
                        "      for ( k in 1:n_latent ) {\n",
                        "        z_mod[i,j] = z_mod[i,j] + x_lat[siteIND[i],k]*beta_lat[k,j];\n",
                        "      }\n",
                        sep="")
  }
  model.code <- paste(model.code,
                      "      p_mod[i,j] = 1/(1+exp(-z_mod[i,j]));\n",
                      "    }\n",
                      "  }\n",
                      "  for ( j in 1:n_taxa ) {\n",
                      "    deviance[j] = 0;\n",
                      "    for ( i in 1:n_samples ) {\n",
                      "      if(y[i,j]>=0) {\n",
                      "        deviance[j] = deviance[j] - 2*(y[i,j]*log(p_mod[i,j])+(1-y[i,j])*log(1-p_mod[i,j]));\n",
                      "      }\n",
                      "    }\n",
                      "  }\n",
                      "}\n\n",
                      sep="")
}

cat(model.code, file=file.model)

# Run model ####
# perform Bayesian inference:
ptm <- proc.time()
res <- stan(file.model,data=data,init=init,iter=sampsize,chains=n.chain,warmup=min(0.5*sampsize,1000),thin=thin)
print(proc.time()-ptm)

# save workspace:
save.image(file=paste("Inv_JSDM_",submodel,".RData",sep=""))

