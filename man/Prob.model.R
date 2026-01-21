Prob.model<-function (X, Model, Par.est0, D = 1)
{
  if (Model == "Rasch" || Model == "2PL" || Model == "3PL" ||
      Model == "4PL" || Model == "1PLAG" || Model == "1PLG") {
    if (Model == "Rasch") {
      Prob = 1/(1 + exp(-(X - Par.est0$B)))
    }
    if (Model == "2PL") {
      Prob = 1/(1 + exp(-D * Par.est0$A * (X - Par.est0$B)))
    }
    if (Model == "3PL") {
      Prob = Par.est0$C + (1 - Par.est0$C)/(1 + exp(-D *
                                                      Par.est0$A * (X - Par.est0$B)))
    }
    if (Model == "4PL") {
      Prob = Par.est0$C + (1 - Par.est0$S - Par.est0$C)/(1 +
                                                           exp(-D * Par.est0$A * (X - Par.est0$B)))
    }
    if (Model == "1PLAG") {
      P.1pl = 1/(1 + exp(-(X - Par.est0$Beta)))
      P.ag = 1/(1 + exp(-(Par.est0$Alpha * X + Par.est0$Gamma)))
      Prob = P.1pl + (1 - P.1pl) * P.ag
    }
    if (Model == "1PLG") {
      P.1pl = 1/(1 + exp(-(X - Par.est0$Beta)))
      P.g = 1/(1 + exp(-(Par.est0$Gamma)))
      Prob = P.1pl + (1 - P.1pl) * P.g
    }
    Prob[Prob >= 0.9999] = 0.9999
    Prob[Prob < 1e-04] = 1e-04

    return(Prob)
  }
  else {
    stop("The Model user specified does not exist!")
  }
}
