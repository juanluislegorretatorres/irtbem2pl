#' @title Calculate the probabilites based on a given model and parameters.
#' @description  Based on the given model, return the correct probabilities
#' of a single examinne with ability X answering each item.
#' @param x A numeric with length=1 consists of an examinee's ability theta.
#' @param Par.est0 A list that consists of item parameters for each item
#' based on the given model. For 2PL model, list(A, B)  are
#' numeric refer to item discrimination  and difficulty
#' parameters for each item, respectively.
#' @param Model caracter declaro para el tipo de modelo a usar.
#' '2PL' - Three parameter logistic (2PL) model proposed by Birnbaum
#' (1968) Adjustment with Marginal Bayesian Modal Item Parameter Estimation.
#' #'
#' P(x=1∣θ,a,b)=(1+exp(−D∗a∗(θ−b)))
#'
#' where x=1 is the correct response, theta is examinne's ability;
#' a, b are the item discrimination, and difficulty parameter,
#' respectively; D is the scaling constant 1 normal.'
#'
#' @param D the scaling constant.  By default 1 normal or 1.702 log.
#'
#' @returns A numeric consists of the correct probabilities of a single
#'examinne with ability X answering each item
#'
#' @references \emph{Item Response Theory Parameter Estimation Techniques}
#' Second Edition, Revised and Expanded Frank B. Baker University ofWisconsin
#' Madison, Wisconsin, U.S.A. Seock-Ho Kim The University ofGeorgia
#' Athens, Georgia, u.S.A.
#'#'
#' Barton, M. A., & Lord, F. M. (1981). An upper asymptote for the
#' three-parameter logistic item response model. ETS Research Report Series,
#' 1981(1), 1-8. doi:10.1002/j.2333-8504.1981.tb01255.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an
#' examinee's ability. In F. M. Lord & M. R. Novick (Eds.), Statistical theories
#' of mental test scores (pp. 395-479). MA: Adison-Wesley.
#'
#' San Martín, E., Del Pino, G., & De Boeck, P. (2006). IRT models for
#' ability-based guessing. Applied Psychological Measurement, 30(3),
#' 183-203. doi:10.1177/0146621605282773
#'@seealso \code{\link{Input.Checking}}, \code{\link{BEMM.3PL}},
#' \code{\link{BEMM.1PLG}}, \code{\link{Prob.model}}
#'
#'
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
