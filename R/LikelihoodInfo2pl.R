#' @title Compute Likelihood Information for the 2PL IRT Model
#'
#' @description  ' Calculates the log-likelihood, deviance, and expected frequencies (f, fz, rz)
#' required for the E-step of the Expectation-Maximization (EM) algorithm
#' under a two-parameter logistic (2PL) model.
#'
#' @param data.simple A numeric matrix representing the simplified unique
#'     response patterns.
#' @param CountNum A numeric vector containing the frequencies of each unique
#'     response pattern.
#' @param Model A character string declaring the type of IRT model to be evaluated
#'     (e.g., "2PL").
#' @param Par.est0 A list containing current item parameter estimates (A and B).
#' @param n.Quadpts A single integer indicating the number of quadrature points.
#' @param node.Quadpts A numeric vector containing the nodes of the quadrature grid.
#' @param weight.Quadpts A numeric vector containing the weights of the quadrature grid.
#' @param D A numeric scaling constant (typically 1 or 1.702).
#'
#' @return A list containing the following statistical components:
#' \itemize{
#'   \item \code{LH} The log-likelihood value of the model given the data.
#'   \item \code{Deviance} The model deviance statistic calculated as -2 * LH.
#'   \item \code{f} A numeric matrix of expected frequencies per quadrature point.
#'   \item \code{fz} A numeric matrix of expected correct response frequencies based on P*.
#'   \item \code{rz} A numeric matrix of expected correct joint frequencies.
#' }
#' @export
#'
LikelihoodInfo2pl<-function (data.simple, CountNum, Model, Par.est0, n.Quadpts,
          node.Quadpts, weight.Quadpts, D)
{
  P.Quadpts = lapply(node.Quadpts, Prob.model, Model = Model,
                     Par.est0 = Par.est0, D = D)
  Joint.prob = mapply("*", lapply(P.Quadpts, function(P, data) {
    apply(data * P + (1 - data) * (1 - P), 2, prod, na.rm = T)
  }, data = data.simple), weight.Quadpts, SIMPLIFY = FALSE)
  Whole.prob = Reduce("+", Joint.prob)
  LH = sum(log(Whole.prob) * CountNum)
  Posterior.prob = lapply(Joint.prob, "/", Whole.prob = Whole.prob)
 # Par.est0$C[Par.est0$C >= 0.5] = 0.4
  f = simplify2array(lapply(lapply(Posterior.prob, "*", CountNum),
                            sum, na.rm = T))
  r = simplify2array(lapply(lapply(lapply(Posterior.prob,
                                          "*", t(data.simple)), "*", CountNum),
                            colSums, na.rm = T))
  if (Model == "2PL" || Model == "3PL"
      ) {
    Pstar = lapply(node.Quadpts, Prob.model, Model = "2PL",
                   Par.est0 = Par.est0, D = D)
    EZ = lapply(mapply("/", Pstar, P.Quadpts, SIMPLIFY = FALSE),
                "*", data.simple)
  }

  if (Model == "1PLG" || Model == "1PLAG") {
    Pstar = lapply(node.Quadpts, Prob.model, Model = "Rasch",
                   Par.est0 = Par.est0, D = D)
    EZ = lapply(mapply("/", Pstar, P.Quadpts, SIMPLIFY = FALSE),
                "*", data.simple)
  }
  EZ.core = mapply("*", lapply(EZ, "t"), Posterior.prob, SIMPLIFY = FALSE)
  fz = simplify2array(lapply(lapply(EZ.core, "*", CountNum),
                             colSums, na.rm = T))
  rz = simplify2array(lapply(lapply(lapply(EZ.core, "*", t(data.simple)),
                                    "*", CountNum), colSums, na.rm = T))
  return(list(LH = LH, Deviance=-2 * LH, f = f, fz = fz, rz = rz))
}
