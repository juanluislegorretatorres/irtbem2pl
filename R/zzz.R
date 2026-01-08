.onAttach <- function(libname, irtbem2pl) {
  packageStartupMessage("¡Hola! Gracias por usar mi librería ", irtbem2pl, ".")
  packageStartupMessage("Esta es la versión alfa", packageVersion(irtbem2pl), ".")
  packageStartupMessage("data: A matrix or data.frame consists of dichotomous data
  (1 for correct and 0 for wrong response)")
  packageStartupMessage("Usage:")
  packageStartupMessage("irtbem2pl( data, PriorA = c(0, 0.25),PriorB = c(0, 4),
  InitialA = 1,  InitialB = 0,  Tol = 0.01,  max.ECycle = 100L,  max.MCycle = 3L,
  n.decimal = 5L,  n.Quadpts = 50L,  Theta.lim = c(-4, 4),  Missing = -9,
  ParConstraint = F,   BiasSE = F,  D = 1) ")
}
