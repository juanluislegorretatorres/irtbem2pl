#' @title Checking user speciflied input variables
#'
#' @description  Verifies whether the input variables provided by the
#' user are correct based on the 2PL model specifications. If acceptable,
#' the function formats and returns them as a combined list; otherwise,
#' it throws an error
#'
#' @param Model A character string declaring the type of item model (e.g., "2PL").
#' @param data A matrix or data.frame consisting of dichotomous data
#'     (1 for correct and 0 for wrong response), with missing data coded
#'     as in \code{Missing}.
#' @param PriorA A numeric vector or matrix specifying the log-normal prior
#'     distribution hyperparameters (mean and variance) for item discrimination.
#'     Defaults to c(0, 0.25).
#' @param PriorB A numeric vector or matrix specifying the normal prior
#'     distribution hyperparameters (mean and variance) for item difficulty.
#'     Defaults to c(0, 4).
#' @param InitialA A single numeric value or vector specifying initial
#'     values for item discrimination parameters. Defaults to 1.
#' @param InitialB A single numeric value or vector specifying initial
#'     values for item difficulty parameters. Defaults to 0.
#' @param Tol A single numeric value specifying the convergence threshold
#'     for E-step cycles. Defaults to 0.01.
#' @param max.ECycle A single integer specifying the maximum number of
#'     E-step cycles. Defaults to 50L.
#' @param max.MCycle A single integer specifying the maximum number of
#'     M-step cycles. Defaults to 10L.
#' @param n.decimal A single integer specifying the number of decimal places
#'     for output formatting. Defaults to 5L.
#' @param n.Quadpts A single integer specifying the number of quadrature
#'     points per dimension. Defaults to 50L.
#' @param Theta.lim A numeric vector of length 2 defining the integration
#'     grid boundaries. Defaults to c(-4, 4).
#' @param Missing A single numeric value indicating missing elements.
#'     Defaults to -9. Cannot be 0 or 1.
#' @param ParConstraint A logical value indicating whether parameters are
#'     constrained within a reasonable range. Defaults to FALSE.
#' @param BiasSE A logical value determining whether to estimate SEs
#'     directly from the inverted Hessian matrix. Defaults to FALSE.
#' @param D A numeric scaling constant (1 or 1.702). Defaults to 1.
#'
#' @return A combined list containing formatted data arrays, sample sizes,
#'     prior specifications, initial parameter values, and convergence controls.
#'
#' @author Juan Luis Legorreta Torres \email{jlegorreta2002@yahoo.com.mx}
#'
#'@references \emph{Item Response Theory Parameter Estimation Techniques}
#' Second Edition, Revised and Expanded Frank B. Baker University ofWisconsin
#' Madison, Wisconsin, U.S.A. Seock-Ho Kim The University ofGeorgia
#' Athens, Georgia, u.S.A.
#'
#' @export
#' @examples
#'data(dat01)
#' library(irtbem2pl)
#'  Checking_2PL<-irtchek2pl(Model = "2PL", dat01, D = 1)
#'
 irtchek2pl<-function (Model = "2PL", data, PriorA = c(0, 0.25),
                       PriorB = c(0, 4), InitialA = 1, InitialB = 0,
                       Tol = 0.01, max.ECycle = 50L, max.MCycle = 10L,
                       n.decimal = 5L ,  n.Quadpts = 50L, Theta.lim = c(-4, 4),
                       Missing = -9, ParConstraint = FALSE ,BiasSE = FALSE,D=1)
{
   D = D
  if (is.data.frame(data)) {
    data = data.matrix(data)
  }
  if (is.matrix(data)) {
    I = as.integer(nrow(data))
    J = as.integer(ncol(data))
    if (I == 1 | J == 1) {
      stop("Error: The ncol and nrow of data must bigger than 1!")
    }
    else {
     # if (sum(is.na(data)) != 0) {
      #  stop("Error1: Some elements in data are not 1, 0 or Missing!")
      #}
      if (sum(!data %in% c(0, 1, NA)) != 0) {
        stop("Error2: Some elements in data are not 1, 0 or Missing!")
      }
      else {
        if (sum(data == Missing) != 0) {
          Index.miss = which(data == Missing, arr.ind = T)
          data[data == Missing] = 0
          PI = rowMeans(data)
          PJ = colMeans(data)
          for (i in 1:nrow(Index.miss)) {
            PI.tmp = PI[Index.miss[i, 1]]
            PJ.tmp = PJ[Index.miss[i, 2]]
            P.correct = 0.5 + PI.tmp - PJ.tmp
            if (is.na(P.correct)) {
              data[Index.miss[i, 1], Index.miss[i, 2]] = 0
            }
            else {
              if (P.correct >= 1) {
                data[Index.miss[i, 1], Index.miss[i,
                                                  2]] = 1
              }
              else {
                if (P.correct <= 0) {
                  data[Index.miss[i, 1], Index.miss[i,
                                                    2]] = 0
                }
                else {
                  data[Index.miss[i, 1], Index.miss[i,
                                                    2]] = rbinom(1, 1, P.correct)
                }
              }
            }
          }
        }
        datafull = as.data.frame(data)
        data.class = data.matrix(aggregate(list(Num = rep(1,
                                                          I)), datafull, length))
        data.simple = t(data.class[, 1:J])
        n.class = as.integer(nrow(data.class))
        CountNum = as.integer(data.class[, J + 1])
        data.list = list(data = data, data.simple = data.simple,
                         CountNum = CountNum, n.class = n.class, I = I,
                         J = J,D=D)
      }
    }
  }
  else {
    stop("Error: The type of data must be a matrix or a data.framework!")
  }
  if (Model == "2PL"
      #|#Model == "4PL"
      ) {
    if (is.numeric(PriorA) & length(PriorA) == 2 & is.na(sum(PriorA)) ==
        FALSE) {
      if (PriorA[2] <= 0) {
        stop("Error: PriorA[2] is the variance, and it must bigger than 0!")
      }
      else {
        PriorA = rep(PriorA, each = J)
      }
    }
    else {
      if (length(PriorA) == 1) {
        if (is.na(PriorA)) {
          PriorA = rep(-9, 2 * J)
        }
        else {
          stop("Error: PriorA must have two input values unless PriorA=NA!")
        }
      }
      else {
        if (is.matrix(PriorA) & length(PriorA) == J *
            2) {
          if (ncol(PriorA) == 2 | nrow(PriorA) == J) {
            if (sum(rowSums(is.na(PriorA)) == 1) !=
                0) {
              stop("Error: The type of PriorA[1] and PriorA[2] are different!")
            }
            if (sum(PriorA[, 2] <= 0, na.rm = T) !=
                0) {
              stop("Error: PriorA[2] is the variance, and it must bigger than 0!")
            }
            PriorA[is.na(PriorA)] = -9
            PriorA = as.numeric(PriorA)
          }
          else {
            stop("Error: The dim of matrix PriorA must be c(n.item, 2)!")
          }
        }
        else {
          stop("Error: The class of PriorA must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!")
        }
      }
    }
    if (is.numeric(PriorB) & length(PriorB) == 2 & is.na(sum(PriorB)) ==
        FALSE) {
      if (PriorB[2] <= 0) {
        stop("Error: PriorB[2] is the variance, and it must bigger than 0!")
      }
      else {
        PriorB = rep(PriorB, each = J)
      }
    }
    else {
      if (length(PriorB) == 1) {
        if (is.na(PriorB)) {
          PriorB = rep(-9, 2 * J)
        }
        else {
          stop("Error: PriorB must have two input values unless PriorB=NA!")
        }
      }
      else {
        if (is.matrix(PriorB) & length(PriorB) == J *
            2) {
          if (ncol(PriorB) == 2 | nrow(PriorB) == J) {
            if (sum(rowSums(is.na(PriorB)) == 1) !=
                0) {
              stop("Error: The type of PriorB[1] and PriorB[2] are different!")
            }
            if (sum(PriorB[, 2] <= 0, na.rm = T) !=
                0) {
              stop("Error: PriorB[2] is the variance, and it must bigger than 0!")
            }
            PriorB[is.na(PriorB)] = -9
            PriorB = as.numeric(PriorB)
          }
          else {
            stop("Error: The dim of matrix PriorB must be c(n.item, 2)!")
          }
        }
        else {
          stop("Error: The class of PriorB must be NA or a numeric with length=2 or a matrix with dim=c(n.item, 2)!")
        }
      }
    }

    if (Model =="2PL") {

        Prior.list = list(PriorA = PriorA, PriorB = PriorB
                        )
    }
  }

  total_score = matrix(rowSums(data), ncol = 1)
  corr0 = t(cor(data, total_score, use = "complete.obs"))
  PassRate = matrix(colSums(data)/I, nrow = 1)
  PassRate[PassRate > 1] = 0.9999
  PassRate[PassRate < 1e-05] = 1e-05
  Zscore = -qnorm(PassRate, 0, 1)
  if (Model == "2PL"
      #| #Model == "4PL"
      ) {
    Initial.B = as.numeric(Zscore/corr0)
    Initial.A = as.numeric(corr0/sqrt(1 - corr0^2))
    Initial.B[Initial.B > 2] = 1.8
    Initial.B[Initial.B < -2] = -1.8
    Initial.A[Initial.A > 2] = 1.8
    Initial.A[Initial.A < 0.3] = 0.5
    if (Model == "2PL") {
      InitialValue = list(InitialA = Initial.A, InitialB = Initial.B
                                                   )
    }

  }

  if (Model == "2PL"
      ) {
    if (is.numeric(InitialA) & length(InitialA) == 1) {
      if (InitialA <= 0) {
        stop("Error: InitialA must bigger than 0!")
      }
      else {
        InitialValue$InitialA = rep(InitialA, J)
      }
    }
    else {
      if (is.numeric(InitialA) & length(InitialA) == J) {
        for (j in 1:J) {
          if (is.na(InitialA[j]) == FALSE) {
            InitialValue$InitialA[j] = InitialA[j]
          }
        }
      }
      else {
        if (length(InitialA) == 1) {
          if (is.na(InitialA) == F) {
            stop("Error: The class of InitialA must be NA or numeric with length= 1 or n.item!")
          }
        }
        else {
          stop("Error: The class of InitialA must be NA or numeric with length= 1 or n.item!")
        }
      }
    }
    if (is.numeric(InitialB) & length(InitialB) == 1) {
      InitialValue$InitialB = rep(InitialB, J)
    }
    else {
      if (is.numeric(InitialB) & length(InitialB) == J) {
        for (j in 1:J) {
          if (is.na(InitialB[j]) == FALSE) {
            InitialValue$InitialB[j] = InitialB[j]
          }
        }
      }
      else {
        if (length(InitialB) == 1) {
          if (is.na(InitialB) == F) {
            stop("Error: The class of InitialB must be NA or numeric with length= 1 or n.item!")
          }
        }
        else {
          stop("Error: The class of InitialB must be NA or numeric with length= 1 or n.item!")
        }
      }
    }

  }
  if (is.numeric(Tol)) {
    if (length(Tol) == 1) {
      if (Tol <= 0) {
        stop("Error: The min of Tol must bigger than 0!")
      }
    }
    else {
      stop("Error: The length of Tol must be 1!")
    }
  }
  else {
    stop("Error: The type of Tol must be numeric!")
  }
  if (length(max.ECycle) == 1) {
    if (is.numeric(length(max.ECycle)) | is.integer(length(max.ECycle))) {
      max.ECycle = as.integer(max.ECycle)
      if (max.ECycle < 1) {
        stop("Error: The min of max.ECycle is 1!")
      }
    }
    else {
      stop("Error: The type of max.ECycle must be integer!")
    }
  }
  else {
    stop("Error: The length of max.ECycle must be 1!")
  }
  if (length(max.MCycle) == 1) {
    if (is.numeric(length(max.MCycle)) | is.integer(length(max.MCycle))) {
      max.MCycle = as.integer(max.MCycle)
      if (max.MCycle < 1) {
        stop("Error: The min of max.MCycle is 1!")
      }
    }
    else {
      stop("Error: The type of max.MCycle must be integer!")
    }
  }
  else {
    stop("Error: The length of max.MCycle must be 1!")
  }
  if (length(n.Quadpts) == 1) {
    if (is.numeric(length(n.Quadpts)) | is.integer(length(n.Quadpts))) {
      n.Quadpts = as.integer(n.Quadpts)
      if (n.Quadpts < 5) {
        stop("Error: The min of n.Quadpts is 5!")
      }
      if (n.Quadpts > 256) {
        warning("Too large n.Quadpts will cause produre to
                be extremely time-consuming!")
      }
    }
    else {
      stop("Error: The type of n.Quadpts must be integer!")
    }
  }
  else {
    stop("Error: The length of n.Quadpts must be 1!")
  }
  if (length(n.decimal) == 1) {
    if (is.numeric(length(n.decimal)) | is.integer(length(n.decimal))) {
      n.decimal = as.integer(n.decimal)
      if (n.decimal < 0) {
        stop("Error: The min of n.decimal is 0!")
      }
      if (n.decimal > 17) {
        stop("Error: The max of n.decimal is 17!")
      }
    }
    else {
      stop("Error: The type of n.decimal must be integer!")
    }
  }
  else {
    stop("Error: The length of n.decimal must be 1!")
  }
  Integer.tot = list(n.Quadpts = n.Quadpts, max.ECycle = max.ECycle,
                     max.MCycle = max.MCycle, n.decimal = n.decimal)
  if (is.numeric(Theta.lim)) {
    if (length(Theta.lim) == 2) {
      if (Theta.lim[1] >= Theta.lim[2]) {
        stop("Error: Theta.lim[1] must bigger than Theta.lim[2]!")
      }
    }
    else {
      stop("Error: The length of Theta.lim must be 2!")
    }
  }
  else {
    stop("Error: The type of Theta.lim must be numeric!")
  }
  if (length(Missing) == 1) {
    if (is.numeric(Missing) | is.na(Missing)) {
      if (Missing == 1 | Missing == 0) {
        stop("Error: The value of Missing cannot be 1 or 0!")
      }
    }
    else {
      stop("Error: The type of Missing must be numeric or NA!")
    }
  }
  else {
    stop("Error: The length of Missing must be 1!")
  }
  if (length(ParConstraint) == 1) {
    if (is.logical(ParConstraint) == F) {
      stop("Error: The type of ParConstraint must be logical!")
    }
    else {
      if (ParConstraint) {
        ParConstraint = list(ParConstraint = 1L)
      }
      else {
        ParConstraint = list(ParConstraint = 0L)
      }
    }
  }
  else {
    stop("Error: The length of ParConstraint must be 1!")
  }
  if (length(BiasSE) == 1) {
    if (is.logical(BiasSE) == F) {
      stop("Error: The type of BiasSE must be logical!")
    }
    else {
      if (BiasSE) {
        BiasSE = list(BiasSE = 1L)
      }
      else {
        BiasSE = list(BiasSE = 0L)
      }
    }
  }
  else {
    stop("Error: The length of BiasSE must be 1!")
  }
  return(c(data.list, Prior.list, InitialValue, Integer.tot,
           ParConstraint, BiasSE))
}
