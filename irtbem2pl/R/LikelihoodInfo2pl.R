#' @title Maximum Likelihood Estimation (MLE)
#' @description  "Likelihood" is a statistical concept that measures how likely
#' it is to observe a data outcome given a set of model parameters.
#' Maximum-likelihood fitting of univariate distributions,
#' allowing parameters to be held fixed if desired
#' @param data.siple Data frame con la base
#' @param CountNum
#' @param Model Two parameter logistic (2PL) model
#' @param Par.est0 A list that consists of item parameters for each item
#' based on the given model. For 2PL model, list(A, B)
#' @param n.Quadpts  A single integer, refers to number of quadrature points
#' per dimension
#' @param node.Quadpts
#' @param weight.Quadpts
#' @param D the scaling constant.  By default 1 normal or 1.702 log.
#'
#' @seealso \code{\link{Input.Checking}},\code{\link{BEMM.3PL}},
#' \code{\link{BEMM.1PLG}},\code{\link{Prob.model}}
#

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
  return(list(LH = LH, f = f, r = r, fz = fz, rz = rz))
}
