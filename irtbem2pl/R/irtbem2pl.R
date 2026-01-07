#' @title Calibrating 2PL model via (a,b), proposed by Birnbaum
#' (1968). Adjustment with Marginal Bayesian Modal Item Parameter Estimation.
#' @description  This function can estimate the item parameters of the 2PL model
#' via  Bayesian statistical methodology include Cornfield
#' (1969), de Finetti (1974), Edwards, Lindman, and Savage (1963), Lindley
#' (1970a, 1970b, 1971), and Novick and Jackson (1974). Examples of the use
#' of these methods in educational settings can be found in Novick and Jackson
#' (1974), Novick, Jackson, Thayer, and Cole (1972), and Rubin (1980). Lord
#' (1986) compared maximum likelihood and Bayesian estimation methods in
#' IRT.
#'@param data Data frame con la base dicotomica (0, 1)
#' @param PriorA The user specified logarithmic normal distribution prior for item
#' discrimation (a) parameters in the 2PL models.
#' A numeric with two hyperparameters
#' mean and variance of logarithmic normal distribution for all a parameters.
#' By default, PriorA=c(0,0.25), which means a log normal prior of mean=0 and
#' variance=0.25 will be used for all item discrimation parameters.
#' @param PriorB The user specified normal distribution prior for item difficulty (b)
#' parameters in the 2PL models.
#' A numeric with two hyperparameters mean and variance of normal distribution for
#' all b parameters. By default, PriorB=c(0,4), which means a normal prior of
#'  mean=0 and variance=4 will be used for all item difficulty parameters.
#' @param InitialA The user specified starting values for item discrimation (a)
#' parameters in the 2PL model. A 1 (default), single number (numeric), refers to set t
#' his number to be the starting values of a for all items
#' @param InitialB The user specified starting values for item difficulty (b)
#'  parameters in the 2PL model. A 0 (default), single number (numeric),
#'  refers to set this number to be the starting values of b for all items.
#' @param Tol A single number (numeric), refers to convergence threshold
#' for E-step cycles; defaults are 0.01.
#' @param max.ECycle A single integer, refers to maximum number of E-step cycles;
#' defaults are 100L.
#' @param max.MCycle A single integer, refers to maximum number of M-step
#' cycles; defaults are 3L.
#' @param n.Quadpts A single integer, refers to number of quadrature points
#' per dimension (must be larger than 5); defaults are 50L.
#' @param n.decimal A single integer, refers to number of decimal places when o
#' utputs results. defaults are 5L.
#' @param Theta.lim A numeric with two number, refers to the range of
#' integration grid for each dimension; default is c(-4, 4).
#' @param Missing A single number (numeric) to indicate which elements are missing;
#' default is -9. The Missing cannot be 0 or 1.
#' @param ParConstraint A logical value to indicate whether estimates
#' parametes in a reasonable range; default is FALSE. If ParConstraint=TRUE:
#' a in [0.001, 6], b in [-6, 6], c in [0.0001, 0.5].
#' @param BiasSE A logical value to determine whether directly estimating
#' SEs from inversed Hession matrix rather than USEM method, default is FALSE
#' @param D the scaling constant 1 normal or 1.702 log.
#'
#' @details Two parameter logistic (2PL) model proposed by Birnbaum (1968):
#'
#' P(x=1∣θ,a,b)=(1+exp(−D∗a∗(θ−b)))
#'
#' where x=1 is the correct response, theta is examinne's ability;
#' a, b are the item discrimination, and difficulty parameter,
#' respectively; D is the scaling constant 1 normal.
#'
#' The implementation of expected a posteriori (EAP) ability estimation
#' has communalities with the E-step of both marginal maximum likelihood /
#' expectation-maximization (MMLE/EM) and Bayesian modal estimation /
#' expectation-maximization (MBE/EM) estimation of item parameters.
#'
#' @returns This function will return a list includes following:
#'* `Est.ItemPars (a,b):` A dataframe consists of the estimates of a, b and c
#' and corresponding estimated standard errors.
#'* `Est.Theta (ability scale:)` A dataframe consists of the estimates of theta
#'and corresponding estimated standard errors (EAP method).
#'* `Loglikelihood:` The loglikelihood.
#'* `Iteration The:` number of iterations.
#'* `EM.Map:` The parameter estimation history of iterations.
#'* `fits.test:` The model fits information includes G2 test, AIC, BIC and RMSEA.
#'* `Elapsed.time:` The running time of the program.
#'* `InitialValues:` The initial values of item parameters.
#'
#' @author Juan Luis Legorreta Torres \email{jlegorreta2002@yahoo.com.mx}
#'
#' @references \emph{Item Response Theory Parameter Estimation Techniques}
#' Second Edition, Revised and Expanded Frank B. Baker University ofWisconsin
#' Madison, Wisconsin, U.S.A. Seock-Ho Kim The University ofGeorgia
#' Athens, Georgia, u.S.A.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an
#' examinee's ability. In F. M. Lord & M. R. Novick (Eds.), Statistical theories
#' of mental test scores (pp. 395-479). MA: Adison-Wesley.
#'
#' Guo, S., & Zheng, C. (2019). The Bayesian Expectation-Maximization-Maximization
#' for the 3PLM. Frontiers in Psychology, 10(1175), 1-11. doi:10.3389/fpsyg.2019.01175
#'
#' Zheng, C., Meng, X., Guo, S., & Liu, Z. (2018). Expectation-Maximization-Maximization:
#' A feasible MLE algorithm for the three-parameter logistic model based on a mixture
#' modeling reformulation. Frontiers in Psychology, 8(2302), 1-10. doi:10.3389/fpsyg.2017.02302
#'
#' @seealso \code{\link{Input.Checking}},\code{\link{BEMM.3PL}},
#' \code{\link{BEMM.1PLG}},\code{\link{Prob.model}}
#'
#'
#' @export irtbem2pl
#' @examples
#'data(dat01)
#' library(irtbem2pl)
#'  mod_2PL<-irtbem2pl(dat01)
#'
#'  mod_2PL$Est.ItemPars       #show item estimates
#'  mod_2PL$Est.Theta          #show ability estimates
#'
irtbem2pl<-function (data,PriorA = c(0, 0.25), PriorB = c(0, 4), InitialA = 1,
                InitialB = 0, Tol = 0.01, max.ECycle = 100L, max.MCycle = 3L,
                n.decimal = 5L, n.Quadpts = 50L,  Theta.lim = c(-4, 4),
                Missing = -9, ParConstraint = F, BiasSE = F, D = 1)
{
  Time.Begin = Sys.time()
  Model = "2PL"
  Check.results = irtchek2pl(Model = Model, data = data,
             PriorA = PriorA, PriorB = PriorB,
             InitialA = InitialA, InitialB = InitialB,
             Tol = Tol, max.ECycle = max.ECycle,
             max.MCycle = max.MCycle, n.Quadpts = n.Quadpts,
             n.decimal = n.decimal, Theta.lim = Theta.lim,
             Missing = Missing, ParConstraint = ParConstraint,
             BiasSE = BiasSE, D = D)
  data = Check.results$data
  data.simple = Check.results$data.simple
  CountNum = Check.results$CountNum
  I = Check.results$I
  J = Check.results$J
  n.class = Check.results$n.class
  PriorA = Check.results$PriorA
  PriorB = Check.results$PriorB
  Prior = list(PriorA = PriorA, PriorB = PriorB)
  InitialA = Check.results$InitialA
  InitialB = Check.results$InitialB
  max.ECycle = Check.results$max.ECycle
  max.MCycle = Check.results$max.MCycle
  n.Quadpts = Check.results$n.Quadpts
  n.decimal = Check.results$n.decimal
  ParConstraint = Check.results$ParConstraint
  BiasSE = Check.results$BiasSE
  Par.est0 = list(A = InitialA, B = InitialB)
  Par.SE0 = list(SEA = InitialA * 0, SEB = InitialB * 0)
  np = J * 3
  D = Check.results$D

  Est.results = irt2pl(Model = Model, data = data, data.simple = data.simple,
                       CountNum = CountNum, n.class = n.class, Prior = Prior,
                       Par.est0 = Par.est0, Par.SE0 = Par.SE0, D = D, np, Tol = Tol,
                       max.ECycle = max.ECycle, max.MCycle = max.MCycle, n.Quadpts = n.Quadpts,
                       n.decimal = n.decimal, Theta.lim = Theta.lim, Missing = Missing,
                       ParConstraint = ParConstraint, BiasSE = BiasSE, I = I,
                       J = J, Time.Begin = Time.Begin)

  if (Est.results$StopNormal == 1) {
    message("PROCEDURE TERMINATED NORMALLY")
  }
  else {
    message("PROCEDURE TERMINATED WITH ISSUES")
  }
  message("2PL_IRTEMM version: 1.0.1")
  message("Item Parameter Calibration. ",Model, "\n")
  message("Quadrature: ", n.Quadpts, " nodes from ", Theta.lim[1],
          " to ", Theta.lim[2], " were used to approximate Gaussian distribution.")
  message("Method for Items: Ability-based Bayesian
              Expectation-Maximization (BEM) Algorithm.")
  if (BiasSE) {
    message("Method for Item SEs: directly estimating SEs from inversed
                Hession matrix.")
    warning("Warning: The SEs maybe not trustworthy!", sep = "")
  }
  else {
    message("Method for Item SEs: Updated SEM algorithm.")
  }
  message("Method for Theta: Expected A Posteriori (EAP).")
  if (Est.results$StopNormal == 1) {
    message("Converged at LL-Change < ", round(Est.results$cr,
                                               6), " after ",
            Est.results$Iteration, " EMM iterations.",
            sep = "")
  }
  else {
    warning("Warning: Estimation cannot converged under current
                max.ECycle and Tol!",
            sep = "")
    warning("Warning: The reults maybe not trustworthy!",
            sep = "")
    message("Terminated at LL-Change = ", round(Est.results$cr,
                                                6), " after ",
            Est.results$Iteration, " EMM iterations.",
            sep = "")
  }
  message("Running time:", Est.results$Elapsed.time, "\n")
  message("Log-likelihood (LL):", as.character(round(Est.results$Loglikelihood,
                                                     n.decimal)))
  message("Estimated Parameters:", as.character(np))
  message("AIB: ", round(Est.results$fits.test$AIC, n.decimal),
          ", BIC: ", round(Est.results$fits.test$BIC, n.decimal),
          ", RMSEA = ", round(Est.results$fits.test$RMSEA, n.decimal))
  message("G2 (", round(Est.results$fits.test$G2.df, n.decimal),
          ") = ", round(Est.results$fits.test$G2, n.decimal),
          ", p = ", round(Est.results$fits.test$G2.P, n.decimal),
          ", G2/df = ", round(Est.results$fits.test$G2.ratio,
                              n.decimal), sep = "")



  return(Est.results)
}
