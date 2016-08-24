data <- pbc2[complete.cases(pbc2$serBilir, pbc2$serChol, pbc2$hepatomegaly, pbc2$year, pbc2$years, pbc2$age, pbc2$drug, pbc2$sex, pbc2$edema, pbc2$status2), ]

data$IDnr <- as.numeric(data$id)
data.id <- data[tapply(row.names(data), data$IDnr, tail, 1), ]

data$IDnr <- rep(seq(1,dim(data.id)[1],1),  tapply(data$serBilir, data$IDnr, length))
data$id <- rep(seq(1,dim(data.id)[1],1),  tapply(data$serBilir, data$IDnr, length))
data.id$IDnr <- 1:dim(data.id)[1]

#################################
# fit the longitudinal outcomes using mixed-effetcs models
fm1 <- lme(log(serBilir) ~ ns(year, 3) + sex , data = data,
           na.action = na.omit,
           random = list(IDnr = pdDiag(form = ~ ns(year, 3))))
lmeObject <- fm1

fm2 <- lme(log(serChol) ~ ns(year, 3) + sex , data = data,
           na.action = na.omit,
           random = list(IDnr = pdDiag(form = ~ ns(year, 3))))
lmeObject2 <- fm2

fm3 <- glmer(hepatomegaly ~ year + sex + (year|id) , data = data,
             na.action = na.omit,
             family = binomial)
lmeObject3 <- fm3

#################################
### Set
timeVar <- "year"
lag <- 0
survMod <- "spline-PH"

Time <- data.id$years ### is the survival time

#################################
# for the continuous longitudinal outcome create the design matrices
id <- data$id 
offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))

# 1st longitudinal outcome
formYx <- formula(lmeObject)
TermsX <- lmeObject$terms
mfX <- model.frame(TermsX, data = data)
X <- model.matrix(formYx, mfX)

formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
mfZ <- model.frame(terms(formYz), data = data)
TermsZ <- attr(mfZ, "terms")
Z <- model.matrix(formYz, mfZ)

data.id <- data[!duplicated(id), ]
data.id[[timeVar]] <- pmax(Time - 0, 0)

mfX.id <- model.frame(TermsX, data = data.id) 
mfZ.id <- model.frame(TermsZ, data = data.id)  
Xtime <- model.matrix(formYx, mfX.id)
Ztime <- model.matrix(formYz, mfZ.id)

# 2nd longitudinal outcome
formYx2 <- formula(lmeObject2)
TermsX2 <- lmeObject2$terms
mfX2 <- model.frame(TermsX2, data = data)
X2 <- model.matrix(formYx2, mfX2)

formYz2 <- formula(lmeObject2$modelStruct$reStruct[[1]])
mfZ2 <- model.frame(terms(formYz2), data = data)
TermsZ2 <- attr(mfZ2, "terms")
Z2 <- model.matrix(formYz2, mfZ2)

data.id <- data[!duplicated(id), ]
data.id[[timeVar]] <- pmax(Time - 0, 0)

mfX.id2 <- model.frame(TermsX2, data = data.id) 
mfZ.id2 <- model.frame(TermsZ2, data = data.id) 
Xtime2 <- model.matrix(formYx2, mfX.id2)
Ztime2 <- model.matrix(formYz2, mfZ.id2)

# 3rd longitudinal outcome
formYx3 <- hepatomegaly ~ year + sex
TermsX3 <- terms(lmeObject3)
mfX3 <- model.frame(TermsX3, data = data)
X3 <- model.matrix(formYx3, mfX3)

formYz3 <- ~ year
mfZ3 <- model.frame(terms(formYz3), data = data)
TermsZ3 <- attr(mfZ3, "terms")
Z3 <- model.matrix(formYz3, mfZ3)

data.id <- data[!duplicated(id), ]
data.id[[timeVar]] <- pmax(Time - 0, 0)

mfX.id3 <- model.frame(TermsX3, data = data.id) 
mfZ.id3 <- model.frame(TermsZ3, data = data.id)  
Xtime3 <- model.matrix(formYx3, mfX.id3)
Ztime3 <- model.matrix(formYz3, mfZ.id3)

#################################
# survival submodel
# design matrices for the survival submodel
WD <- model.matrix(~ -1 + age + sex, data.id)[,c(1,3)]

eventD <- data.id$status2
nT <- length(Time)
zeros <- numeric(nT)

x <- list(X = X, Z = Z, WD = if (survMod == "weibull-PH") {
  if (is.null(WD)) cbind(rep(1, nT), rep(0, nT)) else cbind(1,
                                                            WD)
} else {
  if (is.null(WD)) cbind(rep(0, nT), rep(0, nT)) else {
    if (ncol(WD) == 1) cbind(WD, rep(0, nT)) else WD
  }
})

#################################
# set
y.long <- model.response(mfX, "numeric")
y.long2 <- model.response(mfX2, "numeric")
y.long3 <- model.response(mfX3)
y.long3 <- as.numeric(y.long3)-1
y <- list(y = y.long, offset = offset, logT = log(Time), y2 = y.long2, y3 = y.long3, 
          eventD = eventD, zeros = zeros, lag = lag)

###################################
# for the longitudinal outcomes - design matrices for the 15-point Gauss-Kronrod quadrature rule approximation
gaussKronrod <- JMbayes:::gaussKronrod
wk <- gaussKronrod()$wk
sk <- gaussKronrod()$sk

ordsk <- order(sk)
sk <- sk[ordsk]
wk <- wk[ordsk]

K <- length(sk)


P <- Time/2
st <- outer(P, sk + 1)
id.GK <- rep(seq_along(Time), each = K)

data.id2 <- data.id[id.GK, ]
data.id2[[timeVar]] <- c(t(st))

# 1st longitudinal outcome
mfX <- model.frame(TermsX, data = data.id2)  
mfZ <- model.frame(TermsZ, data = data.id2)    
Xs <- model.matrix(formYx, mfX)
Zs <- model.matrix(formYz, mfZ)

# 2nd longitudinal outcome
mfX2 <- model.frame(TermsX2, data = data.id2) 
mfZ2 <- model.frame(TermsZ2, data = data.id2)   
Xs2 <- model.matrix(formYx2, mfX2)
Zs2 <- model.matrix(formYz2, mfZ2)

# 3rd longitudinal outcome
mfX3 <- model.frame(TermsX3, data = data.id2)  
mfZ3 <- model.frame(TermsZ3, data = data.id2)  
Xs3 <- model.matrix(formYx3, mfX3)
Zs3 <- model.matrix(formYz3, mfZ3)


#################################
# set MCMC details
con <- list(program = "JAGS", n.chains = 1, n.iter = 100000,
            n.burnin = 50000, n.thin = 20, n.adapt = 10000, K = 100,
            C = 5000, working.directory = getwd(), bugs.directory = "C:/Program Files/WinBUGS14/",
            openbugs.directory = NULL, clearWD = TRUE, over.relax = TRUE,
            knots = NULL, ObsTimes.knots = TRUE, lng.in.kn = 5, ordSpline = 4,
            bugs.seed = 1, quiet = FALSE)


#################################
# design matrices for the baseline hazard (a B-splines baseline hazard function is asssumed)
kn <- if (is.null(con$knots)) {
  pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
  pp <- tail(head(pp, -1), -1)
  tt <- if (con$ObsTimes.knots) {
    Time
  } else {  Time[event == 1]    }
  quantile(tt, pp, names = FALSE)
} else {
  con$knots
}
kn <- kn[kn < max(Time)]
rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
con$knots <- rr

W2D <- splineDesign(rr, Time, ord = con$ordSpline)
if (any(colSums(W2D) == 0))
  stop("\nsome of the knots of the B-splines basis are set outside the range",
       "\n   of the observed event times for one of the strata; refit the model",
       "\n   setting the control argument 'equal.strata.knots' to FALSE.")

# design matrices for the baseline hazard for the 15-point Gauss-Kronrod quadrature rule approximation
W2sD <- splineDesign(rr, c(t(st)), ord = con$ordSpline)

x <- c(x, list(W2D = W2D, W2sD = W2sD))

#################################
ncX <- ncol(X)
ncZ <- ncol(Z)
ncWD <- ncol(x$WD)
ncW2D <- ncol(x$W2D)
ncX2 <- ncol(X2)
ncZ2 <- ncol(Z2)
ncX3 <- ncol(X3)
ncZ3 <- ncol(Z3)
C <- con$C
nb <- ncZ + ncZ2 + ncZ3

b <- cbind(data.matrix(ranef(lmeObject)),data.matrix(ranef(lmeObject2)),data.matrix(ranef(lmeObject3)$id))

nY <- nrow(b)
sigma2 <- lmeObject$sigma^2
sigma22 <- lmeObject2$sigma^2

#################################
# priors/hyperpriors
mu0 <- rep(0, (ncZ+ncZ2+ncZ3))

betas <- rep(0, ncX)
var.betas <- rep(con$K, ncX)

betas2 <- rep(0, ncX2)
var.betas2 <- rep(con$K, ncX2)

betas3 <- rep(0, ncX3)
var.betas3 <- rep(con$K, ncX3)

alphasD <- DalphasD <- 0
var.alphasD <- var.DalphasD <- con$K

DAalphasD <- 0
var.DAalphasD <- con$K

alphasD2 <- DalphasD2 <- 0
var.alphasD2 <- var.DalphasD2 <- con$K

DAalphasD2 <- 0
var.DAalphasD2 <- con$K

alphasD3 <- DalphasD3 <- 0
var.alphasD3 <- var.DalphasD3 <- con$K

DAalphasD3 <- 0
var.DAalphasD3 <- con$K


REalphasD <- rep(0, nb)
var.REalphasD <- rep(con$K, (nb)) 

gammasD <- rep(0,(ncWD))
var.gammasD <- rep(con$K, (ncWD))


Bs.gammasD <- rep(0, (ncW2D))
var.Bs.gammasD <- rep(con$K/10, (ncW2D))

#################################
# design matrices for the slope of the 1st longitudinal outcome
extraFormY1 <- list(fixed = ~ 0 + dns(year, 3),
                    random = ~ 0 + dns(year, 3),
                    indFixed = 2:4, indRandom = 2:4)


mfX.derivY1 <- model.frame(terms(extraFormY1$fixed), data = data)
TermsX.derivY1 <- attr(mfX.derivY1, "terms")
mfZ.derivY1 <- model.frame(terms(extraFormY1$random), data = data)
TermsZ.derivY1 <- attr(mfZ.derivY1, "terms")
mfX.deriv.idY1 <- model.frame(TermsX.derivY1, data = data.id)
mfZ.deriv.idY1 <- model.frame(TermsZ.derivY1, data = data.id)
Xtime.derivY1 <- model.matrix(extraFormY1$fixed, mfX.deriv.idY1)
Ztime.derivY1 <- model.matrix(extraFormY1$random, mfZ.deriv.idY1)
XderivY1 <- model.matrix(extraFormY1$fixed, mfX.derivY1)
ZderivY1 <- model.matrix(extraFormY1$random, mfZ.derivY1)


mfX.derivY1 <- model.frame(TermsX.derivY1, data = data.id2)
mfZ.derivY1 <- model.frame(TermsZ.derivY1, data = data.id2)
Xs.derivY1 <- model.matrix(extraFormY1$fixed, mfX.derivY1)
Zs.derivY1 <- model.matrix(extraFormY1$random, mfZ.derivY1)

# design matrices for the slope of the 2nd longitudinal outcome
extraFormY2 <- list(fixed = ~ 0 + dns(year, 3),
                    random = ~ 0 + dns(year, 3),
                    indFixed = 2:4, indRandom = 2:4)

mfX.derivY2 <- model.frame(terms(extraFormY2$fixed), data = data)
TermsX.derivY2 <- attr(mfX.derivY2, "terms")
mfZ.derivY2 <- model.frame(terms(extraFormY2$random), data = data)
TermsZ.derivY2 <- attr(mfZ.derivY2, "terms")
mfX.deriv.idY2 <- model.frame(TermsX.derivY2, data = data.id)
mfZ.deriv.idY2 <- model.frame(TermsZ.derivY2, data = data.id)
Xtime.derivY2 <- model.matrix(extraFormY2$fixed, mfX.deriv.idY2)
Ztime.derivY2 <- model.matrix(extraFormY2$random, mfZ.deriv.idY2)
XderivY2 <- model.matrix(extraFormY2$fixed, mfX.derivY2)
ZderivY2 <- model.matrix(extraFormY2$random, mfZ.derivY2)


mfX.derivY2 <- model.frame(TermsX.derivY2, data = data.id2)
mfZ.derivY2 <- model.frame(TermsZ.derivY2, data = data.id2)
Xs.derivY2 <- model.matrix(extraFormY2$fixed, mfX.derivY2)
Zs.derivY2 <- model.matrix(extraFormY2$random, mfZ.derivY2)

# design matrices for the slope of the 3rd longitudinal outcome
extraFormY3 <- list(fixed = ~  1,
                    random = ~ 1,
                    indFixed = 2, indRandom = 2)


mfX.derivY3 <- model.frame(terms(extraFormY3$fixed), data = data)
TermsX.derivY3 <- attr(mfX.derivY3, "terms")
mfZ.derivY3 <- model.frame(terms(extraFormY3$random), data = data)
TermsZ.derivY3 <- attr(mfZ.derivY3, "terms")
mfX.deriv.idY3 <- model.frame(TermsX.derivY3, data = data.id)
mfZ.deriv.idY3 <- model.frame(TermsZ.derivY3, data = data.id)
Xtime.derivY3 <- model.matrix(extraFormY3$fixed, mfX.deriv.idY3)
Ztime.derivY3 <- model.matrix(extraFormY3$random, mfZ.deriv.idY3)
XderivY3 <- model.matrix(extraFormY3$fixed, mfX.derivY3)
ZderivY3 <- model.matrix(extraFormY3$random, mfZ.derivY3)


mfX.derivY3 <- model.frame(TermsX.derivY3, data = data.id2)
mfZ.derivY3 <- model.frame(TermsZ.derivY3, data = data.id2)
Xs.derivY3 <- model.matrix(extraFormY3$fixed, mfX.derivY3)
Zs.derivY3 <- model.matrix(extraFormY3$random, mfZ.derivY3)

# design matrices for the area of the 1st longitudinal outcome
AextraFormY1 <- list(fixed = ~ -1 + year + ins(year, 3) + I(year * (sex == "female")),
                     random = ~ -1 + year + ins(year, 3),
                     indFixed = 1:5, indRandom = 1:4)


mfXA.derivY1 <- model.frame(terms(AextraFormY1$fixed), data = data)
TermsXA.derivY1 <- attr(mfXA.derivY1, "terms")
mfZA.derivY1 <- model.frame(terms(AextraFormY1$random), data = data)
TermsZA.derivY1 <- attr(mfZA.derivY1, "terms")
mfXA.deriv.idY1 <- model.frame(TermsXA.derivY1, data = data.id)
mfZA.deriv.idY1 <- model.frame(TermsZA.derivY1, data = data.id)
XAtime.derivY1 <- model.matrix(AextraFormY1$fixed, mfXA.deriv.idY1)
ZAtime.derivY1 <- model.matrix(AextraFormY1$random, mfZA.deriv.idY1)
XAderivY1 <- model.matrix(AextraFormY1$fixed, mfXA.derivY1)
ZAderivY1 <- model.matrix(AextraFormY1$random, mfZA.derivY1)


mfXA.derivY1 <- model.frame(TermsXA.derivY1, data = data.id2)
mfZA.derivY1 <- model.frame(TermsZA.derivY1, data = data.id2)
XAs.derivY1 <- model.matrix(AextraFormY1$fixed, mfXA.derivY1)
ZAs.derivY1 <- model.matrix(AextraFormY1$random, mfZA.derivY1)

# design matrices for the area of the 2nd longitudinal outcome
AextraFormY2 <- list(fixed = ~ -1 + year + ins(year, 3) + I(year * (sex == "female")),
                     random = ~ -1 + year + ins(year, 3),
                     indFixed = 1:5, indRandom = 1:4)


mfXA.derivY2 <- model.frame(terms(AextraFormY2$fixed), data = data)
TermsXA.derivY2 <- attr(mfXA.derivY2, "terms")
mfZA.derivY2 <- model.frame(terms(AextraFormY2$random), data = data)
TermsZA.derivY2 <- attr(mfZA.derivY2, "terms")
mfXA.deriv.idY2 <- model.frame(TermsXA.derivY2, data = data.id)
mfZA.deriv.idY2 <- model.frame(TermsZA.derivY2, data = data.id)
XAtime.derivY2 <- model.matrix(AextraFormY2$fixed, mfXA.deriv.idY2)
ZAtime.derivY2 <- model.matrix(AextraFormY2$random, mfZA.deriv.idY2)
XAderivY2 <- model.matrix(AextraFormY2$fixed, mfXA.derivY2)
ZAderivY2 <- model.matrix(AextraFormY2$random, mfZA.derivY2)


mfXA.derivY2 <- model.frame(TermsXA.derivY2, data = data.id2)
mfZA.derivY2 <- model.frame(TermsZA.derivY2, data = data.id2)
XAs.derivY2 <- model.matrix(AextraFormY2$fixed, mfXA.derivY2)
ZAs.derivY2 <- model.matrix(AextraFormY2$random, mfZA.derivY2)


x <- c(x, list(Xs.derivY1 = Xs.derivY1, Zs.derivY1 = Zs.derivY1, XAs.derivY1 = XAs.derivY1, ZAs.derivY1 = ZAs.derivY1,
               Xs.derivY2 = Xs.derivY2, Zs.derivY2 = Zs.derivY2, XAs.derivY2 = XAs.derivY2, ZAs.derivY2 = ZAs.derivY2))

#################################
# 2 stage approach - calculation of the mean and sd of the longitudinal outcomes
# 1st longitudinal outcome
muV <- mean(Xtime%*%fixef(fm1) + rowSums(Ztime*ranef(fm1)))
stdV <- sd(Xtime%*%fixef(fm1) + rowSums(Ztime*ranef(fm1)))

muS <- mean(Xtime.derivY1%*%fixef(fm1)[c(2:4)] + rowSums(Ztime.derivY1*ranef(fm1)[c(2:4)]))
stdS <- sd(Xtime.derivY1%*%fixef(fm1)[c(2:4)] + rowSums(Ztime.derivY1*ranef(fm1)[c(2:4)]))

muA <- mean(XAtime.derivY1%*%fixef(fm1) +  rowSums(ZAtime.derivY1*ranef(fm1)))
stdA <- sd(XAtime.derivY1%*%fixef(fm1) +  rowSums(ZAtime.derivY1*ranef(fm1)))

# 2nd longitudinal outcome
muV2 <- mean(Xtime2%*%fixef(fm2) + rowSums(Ztime2*ranef(fm2)))
stdV2 <- sd(Xtime2%*%fixef(fm2) + rowSums(Ztime2*ranef(fm2)))

muS2 <- mean(Xtime.derivY2%*%fixef(fm2)[c(2:4)] + rowSums(Ztime.derivY2*ranef(fm2)[c(2:4)]))
stdS2 <- sd(Xtime.derivY2%*%fixef(fm2)[c(2:4)] + rowSums(Ztime.derivY2*ranef(fm2)[c(2:4)]))

muA2 <- mean(XAtime.derivY2%*%fixef(fm2) +  rowSums(ZAtime.derivY2*ranef(fm2)))
stdA2 <- sd(XAtime.derivY2%*%fixef(fm2) +  rowSums(ZAtime.derivY2*ranef(fm2)))

# 3rd longitudinal outcome
etahat <- Xtime3%*%fixef(fm3) + rowSums(Ztime3*ranef(fm3)$id)

muV3 <- mean(expit(etahat))
stdV3 <- sd(expit(etahat))

muS3 <- mean((Xtime.derivY3%*%fixef(fm3)[c(2)] + rowSums(Ztime.derivY3*ranef(fm3)$id[c(2)]))*
               (exp(etahat)/((1+exp(etahat))^2)  )   )
stdS3 <- sd((Xtime.derivY3%*%fixef(fm3)[c(2)] + rowSums(Ztime.derivY3*ranef(fm3)$id[c(2)]))*
              (exp(etahat)/((1+exp(etahat))^2)  )   )

#################################
# specify parameters of interest
parms <- c("betas","betas2","betas3", "tau","tau2", "inv.D", "gammasD","b", "alphasD", "DalphasD", "DAalphasD", 
           "alphasD2", "DalphasD2", "DAalphasD2","alphasD3", "DalphasD3", 
           "Bs.gammasD")

#################################
Data <- list(N = nY, K = K, offset = offset, X = X, Xtime = Xtime, 
             Xtime.derivY1 = Xtime.derivY1, XAtime.derivY1 = XAtime.derivY1,
             y = y$y, 
             Xs = Xs, Xs.derivY1 = Xs.derivY1, XAs.derivY1 = XAs.derivY1,
             Z = Z, Ztime = Ztime,  Ztime.derivY1 = Ztime.derivY1, ZAtime.derivY1 = ZAtime.derivY1,
             Zs = Zs, Zs.derivY1 = Zs.derivY1,  ZAs.derivY1 = ZAs.derivY1, 
             
             X2 = X2, Xtime2 = Xtime2, 
             Xtime.derivY2 = Xtime.derivY2, XAtime.derivY2 = XAtime.derivY2,
             y2 = y$y2, 
             Xs2 = Xs2, Xs.derivY2 = Xs.derivY2, XAs.derivY2 = XAs.derivY2,
             Z2 = Z2, Ztime2 = Ztime2,  Ztime.derivY2 = Ztime.derivY2, ZAtime.derivY2 = ZAtime.derivY2,
             Zs2 = Zs2, Zs.derivY2 = Zs.derivY2,  ZAs.derivY2 = ZAs.derivY2, 
             
             X3 = X3, Xtime3 = Xtime3, 
             Xtime.derivY3 = Xtime.derivY3, 
             y3 = y$y3, 
             Xs3 = Xs3, Xs.derivY3 = Xs.derivY3, 
             Z3 = Z3, Ztime3 = Ztime3,  Ztime.derivY3 = Ztime.derivY3, 
             Zs3 = Zs3, Zs.derivY3 = Zs.derivY3,  
             
             eventD = eventD, zeros = zeros, 
             WD = x$WD, ncZ = ncol(Z), 
             ncX = ncol(X), 
             ncZ2 = ncol(Z2), 
             ncX2 = ncol(X2),
             ncZ3 = ncol(Z3), 
             ncX3 = ncol(X3),
             
             ncWD = ncol(x$WD),  
             ncX.derivY1 = ncol(XderivY1), ncZ.derivY1 = ncol(ZderivY1),
             ncXA.derivY1 = ncol(XAderivY1), ncZA.derivY1 = ncol(ZAderivY1),
             indFixed = extraFormY1$indFixed, indRandom = extraFormY1$indRandom,
             AindFixed = AextraFormY1$indFixed, AindRandom = AextraFormY1$indRandom,
             
             ncX.derivY2 = ncol(XderivY2), ncZ.derivY2 = ncol(ZderivY2),
             ncXA.derivY2 = ncol(XAderivY2), ncZA.derivY2 = ncol(ZAderivY2),
             indFixed2 = extraFormY2$indFixed, indRandom2 = extraFormY2$indRandom,
             AindFixed2 = AextraFormY2$indFixed, AindRandom2 = AextraFormY2$indRandom,
             
             ncX.derivY3 = ncol(XderivY3), ncZ.derivY3 = ncol(ZderivY3),
             indFixed3 = extraFormY3$indFixed, indRandom3 = extraFormY3$indRandom,
             
             W2D = W2D, 
             W2sD = W2sD, ncW2D = ncol(x$W2D), C = C, P = P,
             wk = wk, nb = nb, 
             mu0 = mu0, 
             priorMean.betas = betas, 
             priorTau.betas = diag(1/var.betas),
             
             priorMean.betas2 = betas2, 
             priorTau.betas2 = diag(1/var.betas2),
             
             priorMean.betas3 = betas3, 
             priorTau.betas3 = diag(1/var.betas3),
             
             priorA.tau = (1/sigma2)^2/10,
             priorB.tau = (1/sigma2)/10, 
             
             priorA.tau2 = (1/sigma22)^2/10,
             priorB.tau2 = (1/sigma22)/10, 
             
             priorMean.gammas = gammasD,
             priorTau.gammas = diag(1/var.gammasD),
             priorMean.alphas = alphasD,
             priorMean.Dalphas = DalphasD, 
             priorMean.DAalphas = DAalphasD, 
                 
             priorMean.alphas2 = alphasD2,
             priorMean.Dalphas2 = DalphasD2, 
             priorMean.DAalphas2 = DAalphasD2, 
             
             priorMean.alphas3 = alphasD3,
             priorMean.Dalphas3 = DalphasD3, 
             
             priorMean.Bs.gammas = Bs.gammasD,
             priorTau.Bs.gammas = diag(1/var.Bs.gammasD),
             priorR.D = diag(1,(ncZ+ncZ2+ncZ3)), priorK.D = (ncZ+ncZ2+ncZ3), 
             
             muV = muV, stdV = stdV, muS = muS, stdS = stdS, muA = muA, stdA = stdA,             
             muV2 = muV2, stdV2 = stdV2, muS2 = muS2, stdS2 = stdS2, muA2 = muA2, stdA2 = stdA2,
             muV3 = muV3, stdV3 = stdV3, muS3 = muS3, stdS3 = stdS3,
             lamA = 0.1, lamB = 0.1)

