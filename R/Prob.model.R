#' Compute Item Response Probabilities for Various IRT Models
#'
#' Calculates the probability of a correct response for a given latent trait
#' level (Theta) under specified Item Response Theory (IRT) models.
#'
#' @param X A numeric value or vector representing the latent trait level
#'     (Theta) or quadrature node.
#' @param Model A character string specifying the IRT model type. Options
#'     include "Rasch", "2PL", "3PL", "4PL", "1PLAG", or "1PLG".
#' @param Par.est0 A list containing current item parameter estimates
#'     (A, B, C, S, Beta, Alpha, Gamma, depending on the model).
#' @param D A numeric scaling constant (typically 1 or 1.702). Defaults to 1.
#'
#' @return A numeric vector or matrix containing the calculated correct response
#'     probabilities, bounded between 1e-12 and (1 - 1e-12) for numerical stability.
#' @export
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

