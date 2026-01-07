irt2pl = function (Model = Model, data = data, data.simple = data.simple,
                  CountNum = CountNum, n.class = n.class, Prior = Prior,
                  Par.est0 = Par.est0, Par.SE0 = Par.SE0, D = D, np, Tol = Tol,
                  max.ECycle = max.ECycle, max.MCycle = max.MCycle, n.Quadpts = n.Quadpts,
                  n.decimal = n.decimal, Theta.lim = Theta.lim, Missing = Missing,
                  ParConstraint = ParConstraint, BiasSE = BiasSE, I = I,
                  J = J, Time.Begin = Time.Begin)
{
  node.Quadpts = seq(Theta.lim[1], Theta.lim[2], length.out = n.Quadpts)
  weight.Quadpts = dnorm(node.Quadpts, 0, 1, log = F)
  weight.Quadpts = weight.Quadpts/sum(weight.Quadpts)
  InitialValues = Par.est0
  LH = rep(0, max.ECycle)
  IA = rep(0, J)
  IB = rep(0, J)
  IC = rep(0, J)
  IAB = rep(0, J)
  TA = matrix(0, max.ECycle, J)
  TB = matrix(0, max.ECycle, J)
  deltahat.A = rep(0, J)
  deltahat.B = rep(0, J)
  n.ECycle = 1L
  StopNormal = 0L
  E.exit = 0L
  LLinfo = LikelihoodInfo2pl(data.simple, CountNum, Model, Par.est0,
                             n.Quadpts, node.Quadpts, weight.Quadpts, D)
  LH0 = LLinfo$LH
  f = LLinfo$f
  r = LLinfo$r
  fz = LLinfo$fz
  rz = LLinfo$rz
  while (E.exit == 0L && (n.ECycle <= max.ECycle)) {
    for (j in 1:J) {
      at0 = Par.est0$A[j]
      bt0 = Par.est0$B[j]
      n.MCycle = 1
      M.exit = 0L
      if (abs(at0)>30 | abs(bt0)>20 ) {
        M.exit = 1L
      }
      if (abs(at0)<=0.05 & abs(bt0)<=0.05 ) {
        M.exit = 1L
      }

      while (M.exit == 0L && (n.MCycle <= max.MCycle)) {

        Da = D * at0
        D2a2 = Da * Da
        x.bt = node.Quadpts - bt0
        x.bt2 = x.bt * x.bt
        pstar = 1/(1 + exp(-Da * x.bt))
        psf = pstar * f
        wsf = psf * (1 - pstar)  #w
        fz.psf = fz[j, ] - psf
        la1 = Da * sum(fz.psf * x.bt)
        lb1 = Da * sum(fz.psf)
        laa = D2a2 * sum(wsf * x.bt2)
        lbb = D2a2 * sum(wsf)
        lab = -D2a2 * sum(wsf * x.bt)

        if (Prior$PriorA[j] != -9 && Prior$PriorA[j + J] != -9) {
          la1 = la1 - ((log(at0) - Prior$PriorA[j])/Prior$PriorA[j + J])
          laa = laa + 1/Prior$PriorA[j + J]
        }

        if (Prior$PriorB[j] != -9 && Prior$PriorB[j + J] != -9) {
          lb1 = -lb1 #- ((bt0 - Prior$PriorB[j])/Prior$PriorB[j + J])
          lbb = lbb #+ 1/Prior$PriorB[j + J]
        }
        Weight = laa * lbb - lab * lab #Dm

        if (Weight <=0.000099){
          M.exit = 1L
          at0 = exp(at1)
        }

        Iaa = lbb/Weight
        Ibb = laa/Weight
        Iab = lab/Weight

        at1 = log(at0) + (Iaa * la1 - Iab * lb1)
        bt1 = bt0 + (Ibb * lb1 - Iab * la1)

        if (abs(at1)>30 | abs(bt1)>20 ) {
          M.exit = 1L

        }
        if (abs(at1)<=0.05 & abs(bt1)<=0.05 ) {
          M.exit = 1L

        }


        at0 = exp(at1)
        bt0 = bt1
        n.MCycle = n.MCycle + 1
        # }

      }

      if (is.finite(at0) && is.finite(bt0)
      ) {
        if (ParConstraint) {
          if (at0 >= 0.0001 && at0 <= 6 && bt0 >= -6 &&
              bt0 <= 6) {
            Par.est0$A[j] = at0
            Par.est0$B[j] = bt0
            TA[n.ECycle, j] = at0
            TB[n.ECycle, j] = bt0
            IA[j] = Iaa
            IB[j] = Ibb
            IAB[j] = Iab
          }

        }
        else {
          if (at0 >= 0.0001) {
            Par.est0$A[j] = at0
            TA[n.ECycle, j] = at0
          }
          else {
            if (n.ECycle != 1) {
              TA[n.ECycle, j] = TA[n.ECycle - 1, j]
            } else {
              TA[n.ECycle, j] = at0
            }

          }
          Par.est0$B[j] = bt0
          TB[n.ECycle, j] = bt0
          IA[j] = Iaa
          IB[j] = Ibb
          IAB[j] = Iab
        }
      }
      else {
        if (n.ECycle != 1) {
          TA[n.ECycle, j] = TA[n.ECycle - 1, j]
          TB[n.ECycle, j] = TB[n.ECycle - 1, j]

        }
        else {
          TA[n.ECycle, j] = Par.est0$A[j]
          TB[n.ECycle, j] = Par.est0$B[j]

        }
      }
    }

    LLinfo = LikelihoodInfo2pl(data.simple, CountNum, Model,
                               Par.est0, n.Quadpts, node.Quadpts, weight.Quadpts,
                               D)
    LH[n.ECycle] = LLinfo$LH
    f = LLinfo$f
    r = LLinfo$r
    fz = LLinfo$fz
    rz = LLinfo$rz
    cr = LH[n.ECycle] - LH0
    LH0 = LH[n.ECycle]

    if (abs(cr) < Tol) {
      n.ECycle = n.ECycle + 1
      E.exit = 1
      StopNormal = 1L
    }
    else {
      n.ECycle = n.ECycle + 1
    }
  }

  n.ECycle = n.ECycle - 1
  if (BiasSE == FALSE) {
    start.SEM = 0
    end.SEM = n.ECycle
    delta = rep(0, 5)
    delta0 = rep(0, 5)
    delta1 = rep(0, 5)
    cr.SEM0 = 1
    cr.SEM1 = 1
    cr.SEM2 = 1
    cr.SEM3 = 1
    for (i in 1:(n.ECycle - 1)) {
      deltatemp = exp(-(LH[i + 1] - LH[i]))
      if (deltatemp >= 0.9 && deltatemp <= 0.9999) {
        if (cr.SEM0 == 0) {
          end.SEM = i
        }
        else {
          start.SEM = i
          cr.SEM0 = 0
        }
      }
    }
    Time.Mid = Sys.time()
    message(paste("Estimating SEs via USEM algorithm (Requires about ",
                  as.character(round(difftime(Time.Mid, Time.Begin,
                                              units = "auto"), 2)), " secs).", sep = ""),
            "\n")
    for (j in 1:J) {
      z = start.SEM
      SEM.exit = 0
      cr.SEM1 = 1
      cr.SEM2 = 1
      cr.SEM3 = 1
      while (SEM.exit == 0 && z <= end.SEM) {
        for (ParClass in 1:3) {
          if (ParClass == 1 && z >= 2 && cr.SEM1 < sqrt(Tol)) {
            next
          }
          if (ParClass == 2 && z >= 2 && cr.SEM2 < sqrt(Tol)) {
            next
          }
          deltahat = Par.est0

          if (z>0) {
            if (ParClass == 1) {
              deltahat$A[j] = TA[z, j]
            }
            if (ParClass == 2) {
              deltahat$B[j] = TB[z, j]
            }
          }

          LLinfo = LikelihoodInfo2pl(data.simple, CountNum,
                                     Model, deltahat, n.Quadpts, node.Quadpts,
                                     weight.Quadpts, D)
          f = LLinfo$f
          r = LLinfo$r
          fz = LLinfo$fz
          rz = LLinfo$rz
          if (ParClass == 1 || ParClass == 2) {
            at0 = deltahat$A[j]
            bt0 = deltahat$B[j]
            n.MCycle = 1
            M.exit = 0
            if (abs(at0)>30 | abs(bt0)>20 ) {#
              M.exit = 1L
            }
            if (abs(at0)<=0.05 & abs(bt0)<=0.05 ) {#
              M.exit = 1L
            }

            while (M.exit == 0 && (n.MCycle <= max.MCycle)) {

              Da = D * at0
              D2a2 = Da * Da
              x.bt = node.Quadpts - bt0
              x.bt2 = x.bt * x.bt
              pstar = 1/(1 + exp(-Da * x.bt))
              psf = pstar * f
              wsf = psf * (1 - pstar)
              #wsf[wsf < 0.0000009] = 0
              fz.psf = fz[j, ] - psf
              la1 = Da * sum(fz.psf * x.bt)
              lb1 = Da * sum(fz.psf)
              laa = D2a2 * sum(wsf * x.bt2)
              lbb = D2a2 * sum(wsf)
              lab = -D2a2 * sum(wsf * x.bt)

              if (Prior$PriorA[j] != -9 && Prior$PriorA[j + J] != -9) {
                la1 = la1 - ((log(at0) - Prior$PriorA[j])/Prior$PriorA[j + J])
                laa = laa + 1/Prior$PriorA[j + J]
              }

              if (Prior$PriorB[j] != -9 && Prior$PriorB[j + J] != -9) {
                lb1 = -lb1 #- ((bt0 - Prior$PriorB[j])/Prior$PriorB[j + J])
                lbb = lbb #+ 1/Prior$PriorB[j + J]
              }
              Weight = laa * lbb - lab * lab #Dm

              if (Weight <=0.000099){
                M.exit = 1L
              }

              Iaa = lbb/Weight
              Ibb = laa/Weight
              Iab = lab/Weight

              at1 = log(at0) + (Iaa * la1 - Iab * lb1)
              bt1 = bt0 + (Ibb * lb1 - Iab * la1)


              if (abs(at1)>30 | abs(bt1)>20 ) {
                M.exit = 1L
              }
              if (abs(at1)<=0.05 & abs(bt1)<=0.05 ) {
                M.exit = 1L
              }

              #     if (sqrt((at1 - log(at0))^2 + (bt1 - bt0)^2) < 0.01) {
              #      M.exit = 1L
              #       at0 = exp(at1)
              #       bt0 = bt1
              #     }
              #     else {
              at0 = exp(at1)
              bt0 = bt1
              n.MCycle = n.MCycle + 1
              #    }
            }
            if (is.finite(at0) && is.finite(bt0)) {
              if (ParConstraint) {
                if (at0 >= 0.0001 && at0 <= 6 && bt0 >= -6 && bt0 <= 6)
                {
                  deltahat$A[j] = at0
                  deltahat$B[j] = bt0
                }
              }
              else {
                if (at0 >= 0.0001) {
                  deltahat$A[j] = at0
                }
                else {
                }
                deltahat$B[j] = bt0
              }
            }
          }

          if (z>0) {
            if (ParClass == 1) {
              delta1[1] = (log(deltahat$A[j]) -
                             log(Par.est0$A[j]))/(log(TA[z,  j]) -
                                                    log(Par.est0$A[j]) + 1e-05)
              delta1[2] = (deltahat$B[j] -
                             Par.est0$B[j])/(log(TA[z,  j]) -
                                               log(Par.est0$A[j]) + 1e-05)
              break
            }
            if (ParClass == 2) {
              delta1[3] = (log(deltahat$A[j]) -
                             log(Par.est0$A[j]))/(TB[z, j] -
                                                    Par.est0$B[j] + 1e-05)
              delta1[4] = (deltahat$B[j] - Par.est0$B[j])/(TB[z, j] -
                                                             Par.est0$B[j] + 1e-05)
              break
            }
          }

        }
        cr.SEM1 = sqrt((delta1[1] - delta0[1])^2 +
                         (delta1[2] - delta0[2])^2)
        cr.SEM2 = sqrt((delta1[3] - delta0[3])^2 + (delta1[4] -
                                                      delta0[4])^2)
        cr.SEM3 = abs(delta1[5] - delta0[5])
        if (cr.SEM1 < sqrt(Tol) && cr.SEM2 < sqrt(Tol) &&
            cr.SEM3 < sqrt(Tol) && z >= 2) {
          SEM.exit = 1
        }
        else {
          z = z + 1
        }
        delta0[is.finite(delta1)] = delta1[is.finite(delta1)]
      }
      delta[1] = 1 - delta0[1]
      delta[2] = -delta0[2]
      delta[3] = -delta0[3]
      delta[4] = 1 - delta0[4]
      delta[5] = 1 - delta0[5]
      Weight = delta[1] * delta[4] - delta[2] * delta[3]
      delta1[1] = delta[4]/Weight
      delta1[2] = -delta[2]/Weight
      delta1[3] = -delta[3]/Weight
      delta1[4] = delta[1]/Weight
      delta1[5] = 1/delta[5]
      if (is.finite(delta1[1]) == F || delta1[1] <= 0) {
        delta1[1] = 1
      }
      if (is.finite(delta1[2]) == F || delta1[2] <= 0) {
        delta1[2] = 0
      }
      if (is.finite(delta1[3]) == F || delta1[3] <= 0) {
        delta1[3] = 0
      }
      if (is.finite(delta1[4]) == F || delta1[4] <= 0) {
        delta1[4] = 1
      }
      if (is.finite(delta1[5]) == F || delta1[5] <= 0) {
        delta1[5] = 1
      }
      Par.SE0$SEA[j] = sqrt(Par.est0$A[j] * Par.est0$A[j] *
                              IA[j] * delta1[1] + IAB[j] * delta1[3])
      Par.SE0$SEB[j] = sqrt(IB[j] * delta1[4] + IAB[j] *
                              delta1[2])
      #  Par.SE0$SEC[j] = sqrt(IC[j] * delta1[5])
      if (is.finite(Par.SE0$SEA[j])) {
        if (Par.SE0$SEA[j] > 1) {
          Par.SE0$SEA[j] = sqrt(Par.est0$A[j] * Par.est0$A[j] *
                                  IA[j])
        }
      }
      if (is.finite(Par.SE0$SEB[j])) {
        if (Par.SE0$SEB[j] > 1) {
          Par.SE0$SEB[j] = sqrt(IB[j])
        }
      }

    }
  }
  else {
    message("Directly estimating SEs from inversed Hession matrix.",
            "\n")
    for (j in 1:J) {
      Par.SE0$SEA[j] = sqrt(Par.est0$A[j] * Par.est0$A[j] *
                              IA[j])
      Par.SE0$SEB[j] = sqrt(IB[j])
      # Par.SE0$SEC[j] = sqrt(IC[j])
    }
  }

  Par.est0$A = round(Par.est0$A, n.decimal)
  Par.est0$B = round(Par.est0$B, n.decimal)
  Par.SE0$SEA = round(Par.SE0$SEA, n.decimal)
  Par.SE0$SEB = round(Par.SE0$SEB, n.decimal)
  EM.Map = list(Map.A = TA[1:n.ECycle, ], Map.B = TB[1:n.ECycle,]
  )
  Est.ItemPars = as.data.frame(list(est.a = Par.est0$A, est.b = Par.est0$B,
                                    se.a = Par.SE0$SEA, se.b = Par.SE0$SEB
  ))
  P.Quadpts = lapply(as.list(node.Quadpts), Prob.model, Model = Model,
                     Par.est0 = Par.est0, D = D)
  Joint.prob = mapply("*", lapply(P.Quadpts, function(P, data) {
    apply(data * P + (1 - data) * (1 - P), 2, prod, na.rm = T)
  }, data = t(data)), as.list(weight.Quadpts), SIMPLIFY = FALSE)
  Whole.prob = Reduce("+", Joint.prob)
  LogL = sum(log(Whole.prob))
  Posterior.prob = lapply(Joint.prob, "/", Whole.prob)
  EAP.JP = simplify2array(Joint.prob)
  EAP.Theta = rowSums(matrix(1, I, 1) %*% node.Quadpts * EAP.JP)/rowSums(EAP.JP)
  EAP.WP = EAP.JP * simplify2array(lapply(as.list(node.Quadpts),
                                          function(node.Quadpts, Est.Theta) {
                                            (node.Quadpts - Est.Theta)^2
                                          }, Est.Theta = EAP.Theta))
  hauteur = node.Quadpts[2:n.Quadpts] -
    node.Quadpts[1:(n.Quadpts - 1)]
  base.JP = colSums(t((EAP.JP[, 1:(n.Quadpts - 1)] +
                         EAP.JP[, 2:n.Quadpts])/2) * hauteur)
  base.WP = colSums(t((EAP.WP[, 1:(n.Quadpts - 1)] +
                         EAP.WP[,  2:n.Quadpts])/2) * hauteur)
  EAP.Theta.SE = sqrt(base.WP/base.JP)
  Est.Theta = as.data.frame(list(Theta = EAP.Theta, Theta.SE = EAP.Theta.SE))
  N2loglike = -2 * LogL
  AIC = 2 * np + N2loglike
  BIC = N2loglike + log(I) * np
  Theta.uni = sort(unique(EAP.Theta))
  Theta.uni.len = length(Theta.uni)
  G2 = NA
  df = NA
  G2.P = NA
  G2.ratio = NA
  G2.size = NA
  RMSEA = NA
  if (Theta.uni.len >= 11) {
    n.group = 10
    cutpoint = rep(NA, n.group)
    cutpoint[1] = min(Theta.uni) - 0.0001
    cutpoint[11] = max(Theta.uni) + 0.0001
    for (i in 2:n.group) {
      cutpoint[i] = Theta.uni[(i - 1) * Theta.uni.len/n.group]
    }
    Index = cut(EAP.Theta, cutpoint, labels = FALSE)
  }
  if (Theta.uni.len >= 3 & Theta.uni.len < 11) {
    n.group = Theta.uni.len - 1
    cutpoint = rep(NA, n.group)
    cutpoint[1] = min(Theta.uni) - 0.0001
    cutpoint[n.group] = max(Theta.uni) + 0.0001
    if (Theta.uni.len >= 4) {
      for (i in 2:(Theta.uni.len - 2)) {
        cutpoint[i] = Theta.uni[i]
      }
    }
    Index = cut(EAP.Theta, cutpoint, labels = FALSE)
  }
  if (Theta.uni.len >= 3) {
    X2.item = matrix(NA, n.group, J)
    G2.item = matrix(NA, n.group, J)
    Index.Uni = unique(Index)
    for (k in 1:n.group) {
      data.group = data[Index == Index.Uni[k], ]
      Theta.group = EAP.Theta[Index == Index.Uni[k]]
      Obs.P = colMeans(data.group)
      Exp.P = Reduce("+", lapply(Theta.group, Prob.model,
                                 Model = Model, Par.est0 = Par.est0,
                                 D = D))/nrow(data.group)
      Obs.P[Obs.P >= 1] = 0.99999
      Obs.P[Obs.P <= 0] = 1e-05
      Exp.P[Exp.P >= 1] = 0.99999
      Exp.P[Exp.P <= 0] = 1e-05
      X2.item[k, ] = nrow(data.group) * (Obs.P - Exp.P)^2/(Exp.P *
                                                             (1 - Exp.P))
      Odds1 = log(Obs.P/Exp.P)
      Odds2 = log((1 - Obs.P)/(1 - Exp.P))
      G2.item[k, ] = nrow(data.group) * (Obs.P * Odds1 +
                                           (1 - Obs.P) * Odds2)
    }
    X2 = sum(colSums(X2.item, na.rm = T))
    G2 = sum(2 * colSums(G2.item, na.rm = T))
    df = J * (n.group - 3)
    G2.P = 1 - pchisq(G2, df)
    G2.ratio = G2/df
    RMSEA = sqrt(((X2 - df)/(nrow(data) - 1))/X2)
  }
  else {
    warning("The frequence table is too small to do fit tests.")
  }
  fits.test = list(G2 = G2, G2.df = df, G2.P = G2.P, G2.ratio = G2.ratio,
                   RMSEA = RMSEA, AIC = AIC, BIC = BIC)
  Time.End = Sys.time()
  Elapsed.time = paste("Elapsed time:",
                       as.character(round(difftime(Time.End, Time.Begin,
                                                   units = "auto"),
                                          digits = 4)), "secs")


  return(list(Est.ItemPars = Est.ItemPars, Est.Theta = Est.Theta,
                    Loglikelihood = LogL, Iteration = n.ECycle, EM.Map = EM.Map,
                    fits.test = fits.test, Elapsed.time = Elapsed.time,
                    StopNormal = StopNormal, InitialValues = InitialValues,
                    cr = cr))

}


