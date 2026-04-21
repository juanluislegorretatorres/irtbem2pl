Prob.model<-function (X, Model, Par.est0, D = 1)
{
  if (Model == "Rasch" || Model == "2PL" || Model == "3PL" ||
      Model == "4PL" || Model == "1PLAG" || Model == "1PLG") {
    if (Model == "Rasch") {
      Prob = 1/(1 + exp(-(X - Par.est0$B)))
    }
    if (Model == "2PL") {
      exponent <- -D * Par.est0$A * (X - Par.est0$B)
      # Limiting the exponent between -50 and 50 avoids errors of Inf or absolute 0
      exponent <- pmin(pmax(exponent, -50), 50)
      Prob = 1 / (1 + exp(exponent))
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
    Prob <- pmin(pmax(Prob, 1e-12), 1 - 1e-12)
    return(Prob)
  }
  else {
    stop("The Model user specified does not exist!")
  }
}


