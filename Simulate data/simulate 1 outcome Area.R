#rm(list=ls(all=TRUE))
library(JM)

n <- 250 # number of subjects
K <- 10 #10  # number of planned repeated measurements per subject, per outcome
t.max <- 19.5 # maximum follow-up time

################################################

# parameters for the linear mixed effects model
betas <- c("Group0" = 3.8992, "Group1" = 4.0115, "Time1" =  1.4031,
           "Time2" = 2.4292, "Time3" = 2.0728)
sigma.y <- 0.6350 # measurement error standard deviation
#sigma.y <- 0.8034 # measurement error standard deviation

betas <- c("Group0" = 3.9423, "Group1" = 4.0603, "Time1" =  1.3749,
           "Time2" = 2.3923, "Time3" = 2.0488)
sigma.y <- 0.6356 # measurement error standard deviation

# parameters for the survival model
gammas <- c("(Intercept)" = -5.0580, "Group" = 0.1492) # coefficients for baseline covariates
#gammas <- c("(Intercept)" = -4.0580, "Group" = -0.64)#-0.6018) # coefficients for baseline covariates
#gammas <- c("(Intercept)" = -4.0580, "Group" = 0.0492) # coefficients for baseline covariates

alpha <- -0.0050228 #-0.50228 # association parameter - area under the curve
#alpha <- -0.02 #-0.0228

phi <- 1.6458 #1.6458 # shape for the Weibull baseline hazard
mean.Cens <- 24#10 #24 #12 # mean of the exponential distribution for the censoring mechanism

D <- diag(c(0.7028, 2.2040, 1.3660, 0.2018)^2)
#D <- diag(c(0.6968, 2.1259, 1.5260, 1.2329)^2)
################################################

Bkn <- c(0, 19.5)
kn <- c(2.1, 5.5)

# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
age <- rnorm(n, 46.39673, 13.69578 )
DF <- data.frame(year = times, drug = factor(rep(group, each = K)), age = rep(age, each = K))
X <- model.matrix(~ 0 + drug + ns(year, knots = kn, Boundary.knots = Bkn), data = DF)
Z <- model.matrix(~ ns(year, knots = kn, Boundary.knots = Bkn), data = DF)


typeOperation <- sample(0:1, 300, replace = T)
# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Group" = group)

################################################

#simulate random effects
library(MASS)
b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
y <- rnorm(n * K, eta.y, sigma.y)

# simulate event times
eta.t <- as.vector(as.matrix(W) %*% gammas)
invS <- function (t, u, i) {
  h <- function (s) {
    group0 <- 1 - group[i]
    group1 <- group[i]
    ages <- age[i] 

    INS <- ins(s, knots = kn, Boundary.knots = Bkn)
    XXi <- cbind(group0*s, group1*s, INS[, 1], INS[, 2], INS[, 3])
    ZZi <- cbind(s, INS)
    f1 <- as.vector(XXi %*% betas + rowSums(ZZi * b[rep(i, nrow(ZZi)), ]))
    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha)
  }
  integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(n)
trueTimes <- numeric(n)
for (i in 1:n) {
  Up <- 5000
  tries <- 5
  Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  while(inherits(Root, "try-error") && tries > 0) {
    tries <- tries - 1
    Up <- Up + 5000
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  }
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}
na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

# simulate censoring times from an exponential distribution,
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- runif(n, 0, 2 * mean.Cens)
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator


################################################


# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
dat.id <- data.frame(IDnr = unique(id), years = Time, status = event, group = W[, 2])
names(dat) <- c("year", "group", "age","IDnr", "serBilir", "years", "status")

data <- dat
data.id <- dat.id

## extract 10 subjects for whom predictions will be be made
#ids <- unique(dat$id)[dat.id$event == 0 & tapply(dat$id, dat$id, length) > 2]
#samp.id <- sample(ids, 10)
#predDat <- dat[dat$id %in% samp.id, ]
#dat <- dat[!dat$id %in% samp.id, ]
#dat.id <- dat.id[-samp.id, ]



