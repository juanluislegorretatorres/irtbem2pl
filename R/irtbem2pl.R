#' @title  Calibrate 2PL Model via Marginal Bayesian Modal Estimation
#' @description ' Estimates item parameters of the two-parameter logistic (2PL)
#' model using marginal Bayesian modal estimation via the Expectation-Maximization
#' (EM) algorithm. The results are similar to those of "BILOG-MG",
#' from the book Item Response Theory Parameter Estimation Techniques.
#' @param data A matrix or data.frame consisting of dichotomous data
#'     (1 for correct and 0 for wrong response).
#' @param PriorA A numeric vector with two hyperparameters: the mean and
#'     variance of the log-normal distribution for all a parameters.
#'     Defaults to c(0, 0.25).
#' @param PriorB A numeric vector with two hyperparameters: the mean and
#'     variance of the normal distribution prior for item difficulty (b)
#'     parameters. Defaults to c(0, 4).
#' @param InitialA A single numeric value or vector specifying initial
#'     values for the item discrimination parameters (a). Defaults to 1.
#'     It is recommended to use the point-biserial correlation from
#'     Classical Test Theory.
#' @param InitialB A single numeric value or vector specifying initial
#'     values for the difficulty parameters (b). Defaults to 0.
#'     It is recommended to use the item difficulty index from
#'     Classical Test Theory.
#' @param Tol A single numeric value specifying the convergence threshold
#'     for E-step cycles. Defaults to 0.01.
#' @param max.ECycle A single integer specifying the maximum number of
#'     E-step cycles. Defaults to 50L.
#' @param max.MCycle A single integer specifying the maximum number of
#'     M-step cycles. Defaults to 10L.
#' @param n.Quadpts A single integer specifying the number of quadrature
#'     points per dimension (must be larger than 5). Defaults to 50L.
#' @param n.decimal A single integer specifying the number of decimal
#'     places for the output results. Defaults to 5L.
#' @param Theta.lim A numeric vector of length 2 defining the integration
#'     grid range for each dimension. Defaults to c(-4, 4).
#' @param Missing A single numeric value indicating missing elements.
#'     Defaults to -9. Cannot be 0 or 1.
#' @param ParConstraint A logical value indicating whether to restrict
#'     estimates within a reasonable range. Defaults to FALSE.
#' @param BiasSE A logical value determining whether to estimate SEs
#'     directly from the inverted Hessian matrix. Defaults to FALSE.
#' @param D A numeric scale constant (1 or 1.702). Defaults to 1.
#'
#' @details The two-parameter logistic (2PL) model proposed by Birnbaum (1968)
#' is defined as:
#' \deqn{P(x=1|\Theta,a,b) = \frac{1}{1 + exp(-D \cdot a \cdot (\Theta - b))}}
#' where x=1 is the correct response, \code{Theta} is the examinee's ability,
#' \code{a} is the item discrimination, and \code{b} is the difficulty parameter.
#'
#' @return A list containing the following components:
#' \itemize{
#'   \item \code{Est.ItemPars} A data frame with estimates of parameters a and b, along with their standard errors.
#'   \item \code{Est.Theta} A data frame with EAP ability estimates and their standard errors.
#'   \item \code{Loglikelihood} The final log-likelihood value.
#'   \item \code{Iteration} The number of iterations until convergence.
#'   \item \code{EM.Map} The parameter estimation history across iterations.
#'   \item \code{fits.test} Model fit statistics including G2, AIC, BIC, and RMSEA.
#'   \item \code{Elapsed.time} Total execution time of the program.
#'   \item \code{InitialValues} The initial values used for item calibration.
#' }
#'
#' @author Juan Luis Legorreta Torres \email{jlegorreta2002@@yahoo.com.mx}
#'
#' @references
#' Baker, F. B., & Kim, S.-H. (2004). \emph{Item Response Theory: Parameter
#' Estimation Techniques} (2nd ed.). Marcel Dekker.
#'
#' Birnbaum, A. (1968). Some latent trait models and their use in inferring an
#' examinee's ability. In F. M. Lord & M. R. Novick (Eds.), \emph{Statistical
#' theories of mental test scores} (pp. 395-479). Addison-Wesley.
#'
#' @export
#' @examples
#' data(dat01)
#'
#' # It is recommended to use:
#' # InitialA = point-biserial correlation (Classical Test Theory)
#' # InitialB = item difficulty index (Classical Test Theory)
#'
#' mod_2PL <- irtbem2pl(data = dat01, PriorA = c(0, 0.25), PriorB = c(0, 4),
#'                      InitialA = 1, InitialB = 0, Tol = 0.01, max.ECycle = 50L,
#'                      max.MCycle = 10L, n.decimal = 5L, n.Quadpts = 50L,
#'                      Theta.lim = c(-4, 4), Missing = -9, ParConstraint = FALSE,
#'                      BiasSE = FALSE, D = 1)
#'
#' mod_2PL$Est.ItemPars       # show item estimates
#' mod_2PL$Est.Theta          # show ability estimates

irtbem2pl<-function (data,PriorA = c(0, 0.25), PriorB = c(0, 4), InitialA = 1,
                     InitialB = 0, Tol = 0.01, max.ECycle = 50L, max.MCycle = 10L,
                     n.decimal = 5L, n.Quadpts = 50L,  Theta.lim = c(-4, 4),
                     Missing = -9, ParConstraint = FALSE, BiasSE = FALSE, D = 1)
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
  np = J *2
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
  message("Item Parameter Calibration. ",Model)
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
  message("Running time:", Est.results$Elapsed.time)
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
