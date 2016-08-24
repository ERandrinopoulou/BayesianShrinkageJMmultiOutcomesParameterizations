model <- function ()
{
for (i in 1:N) {
  for (j in offset[i]:(offset[i + 1] - 1)) {
    muy[j] <- inprod(betas[1:ncX], X[j, 1:ncX]) + inprod(b[i, 
                                                           1:ncZ], Z[j, 1:ncZ])
    y[j] ~ dnorm(muy[j], tau)
  }
  for (j in offset[i]:(offset[i + 1] - 1)) {
    muy2[j] <- inprod(betas2[1:ncX2], X2[j, 1:ncX2]) + 
      inprod(b[i, (ncZ + 1):(ncZ + ncZ2)], Z2[j, 1:ncZ])
    y2[j] ~ dnorm(muy2[j], tau2)
  }
  for (j in offset[i]:(offset[i + 1] - 1)) {
    muy3[j] <- inprod(betas3[1:ncX3], X3[j, 1:ncX3]) + 
      inprod(b[i, (ncZ + ncZ2 + 1):(ncZ + ncZ2 + ncZ3)], 
             Z3[j, 1:ncZ3])
    Pr[j] <- max(1.00000E-05, min(0.99999, (exp(muy3[j])/(1 + 
                                                            exp(muy3[j])))))
    y3[j] ~ dbin(Pr[j], 1)
  }
  etaBaselineD[i] <- inprod(gammasD[1:(ncWD)], WD[i, 1:ncWD])
  log.h0.TD[i] <- inprod(Bs.gammasD[1:(ncW2D)], W2D[i, 
                                                    1:ncW2D])
  f.T[i] <- inprod(betas[1:ncX], Xtime[i, 1:ncX]) + inprod(b[i, 
                                                             1:ncZ], Ztime[i, 1:ncZ])
  f.T.derivY1[i] <- inprod(betas[2:4], Xtime.derivY1[i, 
                                                     1:ncX.derivY1]) + inprod(b[i, 2:4], Ztime.derivY1[i, 
                                                                                                       1:ncZ.derivY1])
  fA.T.derivY1[i] <- inprod(betas[1:5], XAtime.derivY1[i, 
                                                       1:ncXA.derivY1]) + inprod(b[i, 1:4], ZAtime.derivY1[i, 
                                                                                                           1:ncZA.derivY1])
  f.T2[i] <- inprod(betas2[1:ncX2], Xtime2[i, 1:ncX2]) + 
    inprod(b[i, (ncZ + 1):(ncZ + ncZ2)], Ztime2[i, 1:ncZ2])
  f.T.derivY2[i] <- inprod(betas2[2:4], Xtime.derivY2[i, 
                                                      1:ncX.derivY2]) + inprod(b[i, 6:8], Ztime.derivY2[i, 
                                                                                                        1:ncZ.derivY2])
  fA.T.derivY2[i] <- inprod(betas2[1:5], XAtime.derivY2[i, 
                                                        1:ncXA.derivY2]) + inprod(b[i, (ncZ + 1):(ncZ + ncZ2)], 
                                                                                  ZAtime.derivY2[i, 1:ncZA.derivY2])
  f.T3[i] <- inprod(betas3[1:ncX3], Xtime3[i, 1:ncX3]) + 
    inprod(b[i, (ncZ + ncZ2 + 1):(ncZ + ncZ2 + ncZ3)], 
           Ztime3[i, 1:ncZ3])
  f.T.derivY3[i] <- inprod(betas3[2], Xtime.derivY3[i, 
                                                    1:ncX.derivY3]) + inprod(b[i, 10], Ztime.derivY3[i, 
                                                                                                     1:ncZ.derivY3])
 log.hazardD[i] <- log.h0.TD[i] + etaBaselineD[i] + alphasD * 
    ((f.T[i] - muV)/stdV) + DalphasD * ((f.T.derivY1[i] - 
                                           muS)/stdS) + DAalphasD * ((fA.T.derivY1[i] - muA)/stdA) + 
    alphasD2 * ((f.T2[i] - muV2)/stdV2) + DalphasD2 * 
    ((f.T.derivY2[i] - muS2)/stdS2) + DAalphasD2 * ((fA.T.derivY2[i] - 
                                                       muA2)/stdA2) + alphasD3 * ((ilogit(f.T3[i]) - muV3)/stdV3) + 
    DalphasD3 * ((f.T.derivY3[i] * (exp(f.T3[i])/(1 + 
                                                    exp(f.T3[i]))^2) - muS3)/stdS3)
  for (k in 1:K) {
    log.h0.sD[i, k] <- inprod(Bs.gammasD[1:(ncW2D)], 
                              W2sD[K * (i - 1) + k, 1:ncW2D])
    f.s[i, k] <- inprod(betas[1:ncX], Xs[K * (i - 1) + 
                                           k, 1:ncX]) + inprod(b[i, 1:ncZ], Zs[K * (i - 
                                                                                      1) + k, 1:ncZ])
    f.s.derivY1[i, k] <- inprod(betas[2:4], Xs.derivY1[K * 
                                                         (i - 1) + k, 1:ncX.derivY1]) + inprod(b[i, 2:4], 
                                                                                               Zs.derivY1[K * (i - 1) + k, 1:ncZ.derivY1])
    fA.s.derivY1[i, k] <- inprod(betas[1:5], XAs.derivY1[K * 
                                                           (i - 1) + k, 1:ncXA.derivY1]) + inprod(b[i, 1:4], 
                                                                                                  ZAs.derivY1[K * (i - 1) + k, 1:ncZA.derivY1])
    f.s2[i, k] <- inprod(betas2[1:ncX2], Xs2[K * (i - 
                                                    1) + k, 1:ncX2]) + inprod(b[i, (ncZ + 1):(ncZ + 
                                                                                                ncZ2)], Zs2[K * (i - 1) + k, 1:ncZ2])
    f.s.derivY2[i, k] <- inprod(betas2[2:4], Xs.derivY2[K * 
                                                          (i - 1) + k, 1:ncX.derivY2]) + inprod(b[i, 6:8], 
                                                                                                Zs.derivY2[K * (i - 1) + k, 1:ncZ.derivY2])
    fA.s.derivY2[i, k] <- inprod(betas2[1:5], XAs.derivY2[K * 
                                                            (i - 1) + k, 1:ncXA.derivY2]) + inprod(b[i, (ncZ + 
                                                                                                           1):(ncZ + ncZ2)], ZAs.derivY2[K * (i - 1) + k, 
                                                                                                                                         1:ncZA.derivY2])
    f.s3[i, k] <- inprod(betas3[1:ncX3], Xs3[K * (i - 
                                                    1) + k, 1:ncX3]) + inprod(b[i, (ncZ + ncZ2 + 
                                                                                      1):(ncZ + ncZ2 + ncZ3)], Zs3[K * (i - 1) + k, 
                                                                                                                   1:ncZ3])
    f.s.derivY3[i, k] <- inprod(betas3[2], Xs.derivY3[K * 
                                                        (i - 1) + k, 1:ncX.derivY3]) + inprod(b[i, 10], 
                                                                                              Zs.derivY3[K * (i - 1) + k, 1:ncZ.derivY3])
    SurvLongD[i, k] <- wk[k] * exp(log.h0.sD[i, k] + 
                                     alphasD * ((f.s[i, k] - muV)/stdV) + DalphasD * 
                                     ((f.s.derivY1[i, k] - muS)/stdS) + DAalphasD * 
                                     ((fA.s.derivY1[i, k] - muA)/stdA) + alphasD2 * 
                                     ((f.s2[i, k] - muV2)/stdV2) + DalphasD2 * ((f.s.derivY2[i, 
                                                                                             k] - muS2)/stdS2) + DAalphasD2 * ((fA.s.derivY2[i, 
                                                                                                                                             k] - muA2)/stdA2) + alphasD3 * ((ilogit(f.s3[i, 
                                                                                                                                                                                          k]) - muV3)/stdV3) + DalphasD3 * ((f.s.derivY3[i, 
                                                                                                                                                                                                                                         k] * (exp(f.s3[i, k])/(1 + exp(f.s3[i, k]))^2) - 
                                                                                                                                                                                                                               muS3)/stdS3))
  }
  log.survivalD[i] <- -exp(etaBaselineD[i]) * P[i] * sum(SurvLongD[i, 
                                                                   ])
  phi[i] <- C - ((eventD[i] * log.hazardD[i])) - (log.survivalD[i])
  zeros[i] ~ dpois(phi[i])
  b[i, 1:nb] ~ dmnorm(mu0[], inv.D[, ])
}
betas[1:ncX] ~ dmnorm(priorMean.betas[], priorTau.betas[, 
                                                        ])
betas2[1:ncX2] ~ dmnorm(priorMean.betas2[], priorTau.betas2[, 
                                                            ])
betas3[1:ncX3] ~ dmnorm(priorMean.betas3[], priorTau.betas3[, 
                                                            ])
tau ~ dgamma(priorA.tau, priorB.tau)
tau2 ~ dgamma(priorA.tau2, priorB.tau2)
gammasD[1:(ncWD)] ~ dmnorm(priorMean.gammas[], priorTau.gammas[, 
                                                               ])
alphasD ~ dnorm(priorMean.alphas, priorTau.alphas)
DalphasD ~ dnorm(priorMean.Dalphas, priorTau.Dalphas)
DAalphasD ~ dnorm(priorMean.DAalphas, priorTau.DAalphas)
alphasD2 ~ dnorm(priorMean.alphas2, priorTau.alphas2)
DalphasD2 ~ dnorm(priorMean.Dalphas2, priorTau.Dalphas2)
DAalphasD2 ~ dnorm(priorMean.DAalphas2, priorTau.DAalphas2)
alphasD3 ~ dnorm(priorMean.alphas3, priorTau.alphas3)
DalphasD3 ~ dnorm(priorMean.Dalphas3, priorTau.Dalphas3)
priorTau.alphas <- pow(sigma.alphas, -1)
priorTau.Dalphas <- pow(sigma.Dalphas, -1)
priorTau.DAalphas <- pow(sigma.DAalphas, -1)
priorTau.alphas2 <- pow(sigma.alphas2, -1)
priorTau.Dalphas2 <- pow(sigma.Dalphas2, -1)
priorTau.DAalphas2 <- pow(sigma.DAalphas2, -1)
priorTau.alphas3 <- pow(sigma.alphas3, -1)
priorTau.Dalphas3 <- pow(sigma.Dalphas3, -1)
sigma.alphas ~ dexp(lamda^2/2)
sigma.Dalphas ~ dexp(lamda^2/2)
sigma.DAalphas ~ dexp(lamda^2/2)
sigma.alphas2 ~ dexp(lamda^2/2)
sigma.Dalphas2 ~ dexp(lamda^2/2)
sigma.DAalphas2 ~ dexp(lamda^2/2)
sigma.alphas3 ~ dexp(lamda^2/2)
sigma.Dalphas3 ~ dexp(lamda^2/2)
lamdaPar ~ dgamma(lamA, lamB)
lamda <- sqrt(2 * lamdaPar)
Bs.gammasD[1:(ncW2D)] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, 
                                                                         ])
inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
