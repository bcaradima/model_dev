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
community <- "community.csv"
predictors <- "predictors.csv"
inf.fact <- c("A10m", "IAR", "LUD", "Urban", "bFRI", "FRI", "FV", "WV", "Temp", "Temp2")

random.effects <- FALSE

crossvalid    <- FALSE
dir.modinp    <- "."
dir.modout    <- "."

sampsize      <- 3000
n.chain       <- 4
prob.defpri   <- 0.02
thresh.sig    <- 1

# Packages ####
if ( !require("coda") )  { install.packages("coda");  library("coda") }
if ( !require("rstan") ) { install.packages("rstan"); library("rstan") }
if ( !require("dplyr") ) { install.packages("dplyr"); library("dplyr") }
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Read data ####
# read occurrence data:
occur.taxa <- read.table(paste("inputs", community, sep="/"), header=TRUE, sep=",",
                         stringsAsFactors=FALSE)

# read environmental conditions:
env.cond <- read.table(paste("inputs", predictors, sep="/"), header=TRUE, sep=",",
                       stringsAsFactors=FALSE)

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

# Model D1 ####

file.model.D1 <- "Inv_JSDM_D1.stan"

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
             n_taxa                  = ncol(occur.taxa),
             n_pred                  = length(inf.fact),
             mu_alpha_comm_pripar    = mu.alpha.comm.pripar,
             sigma_alpha_comm_pripar = sigma.alpha.comm.pripar,
             mu_beta_comm_pripar     = mu.beta.comm.pripar,
             sigma_beta_comm_pripar  = sigma.beta.comm.pripar,
             fact_sd                 = 1,
             x                       = as.matrix(env.cond[set.ident,inf.fact]),
             y                       = as.matrix(occur.taxa.na.encod[set.ident,]))

if (random.effects) {
  data$sigma_eps_pripar <- 2*logit(1-prob.defpri)/10
  data$siteIND <- siteIND
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
}


### Model definition ####
model.code <- paste(
"// Model D1: logistic, hierarchical JSDM\n",
"// -----------------------------------------------------------------------\n\n",
"// data:\n\n",
"data {\n",
"  int                     n_taxa;\n",
"  int                     n_sites;\n",
"  int                     n_samples;\n",
"  int                     n_pred;\n",
"  real                    mu_alpha_comm_pripar;\n",
"  real                    sigma_alpha_comm_pripar;\n",
"  vector[n_pred]          mu_beta_comm_pripar;\n",
"  vector[n_pred]          sigma_beta_comm_pripar;\n",
"  real                    fact_sd;\n",
sep = "")
if (random.effects)
{
  model.code <- paste(model.code,
  "  real                    sigma_eps_pripar;\n",
  "  int                     siteIND[n_samples];\n",
  sep="")
}
model.code <- paste(model.code,
"  matrix[n_samples,n_pred]  x;\n",
"  int<lower=-1,upper=1>   y[n_samples,n_taxa];\n",
"}\n\n", 
sep="")

model.code <- paste(model.code,
"// parameters:\n\n",
"parameters {\n",
"  real                    mu_alpha_comm;\n",
"  real<lower=0>           sigma_alpha_comm;\n",
"  vector[n_pred]          mu_beta_comm;\n",
"  vector<lower=0>[n_pred] sigma_beta_comm;\n",
sep="")
if (random.effects)
{
  model.code <- paste(model.code,
  "  real<lower=0>           sigma_eps;\n",
  "  vector[n_sites]         eps;\n",
  sep="")
}
model.code <- paste(model.code,
"  vector[n_taxa]          alpha_taxa;\n",
"  matrix[n_pred,n_taxa]   beta_taxa;\n",
"}\n\n", 
sep="")

model.code <- paste(model.code,
"// model definition:\n\n",
"model {\n\n",

"  // local variables:\n\n",
"  real z[n_samples,n_taxa];\n",
"  real p[n_samples,n_taxa];\n",
"  real s_pripar;\n\n",
sep="")
  
if (random.effects)
{
  model.code <- paste(model.code,
  "  // root nodes:\n\n",
  "  s_pripar = sqrt(log(1+fact_sd^2));\n",
  "  sigma_eps ~ lognormal(log(sigma_eps_pripar)-0.5*s_pripar^2,s_pripar);\n",
  "  for ( i in 1:n_sites ) {\n",
  "    eps[i] ~ normal(0,sigma_eps);\n",
  "  }\n\n",
  sep="")
}
  
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

"  // intermediate nodes:\n\n",
"  for ( j in 1:n_taxa ) {\n",
"    alpha_taxa[j] ~ normal(mu_alpha_comm,sigma_alpha_comm);\n",
"    for ( k in 1:n_pred ) {\n",
"      beta_taxa[k,j] ~ normal(mu_beta_comm[k],sigma_beta_comm[k]);\n",
"    }\n",
"  }\n\n",
sep="")

model.code <- paste(model.code,
"  // end nodes:\n\n",
"  for ( i in 1:n_samples ) {\n",
"    for ( j in 1:n_taxa ) {\n",
"      z[i,j] = alpha_taxa[j];\n",
"      for ( k in 1:n_pred ) {\n",
"        z[i,j] = z[i,j] + beta_taxa[k,j]*x[i,k];\n",
"      }\n",
sep="")
if (random.effects)
{
model.code <- paste(model.code,
"      z[i,j] = z[i,j] + eps[siteIND[i]];\n",        
sep="")
}
model.code <- paste(model.code,
"      p[i,j] = 1/(1+exp(-z[i,j]));\n",
"      if(y[i,j]>=0) {\n",
"        y[i,j] ~ bernoulli(p[i,j]);\n",
"      }\n",
"    }\n",
"  }\n",
"}\n",
"\n",
sep="")

cat(model.code, file=file.model.D1)

# Run model ####
# perform Bayesian inference:
ptm <- proc.time()
res.D1 <- stan(file.model.D1,data=data,init=init,iter=sampsize,chains=n.chain,warmup=500)
print(proc.time()-ptm)

# save workspace:
save.image(file="Inv_JSDM_D1.RData")