#############################################################################
#' Dynamic promoter activity.
#'
#' \code{DIPA} returns a vector of promoter activity over specified
#' times as described in PUBLICATION NAME.
#'
#' This function requires reporter protein data and corresponding time points,
#' and biomass data and corresponding time points stored in R vectors or similar
#' data structure. User can choose the nature of input (fluorescent protein,
#' other protein or mRNA) data via the data.type input. For single-cell data,
#' input biomass as a vector of 1's of corresponding length to gene expression
#' measurements, and consider promoter activity output to be activity per cell.
#'
#' @param x Vector of times over which to calculate promoter activity (in hours).
#' @param t.f Vector of times of reporter protein (or mRNA) measurements (in hours).
#' @param f.data Vector of of reporter protein (or mRNA) measurements.
#' @param t.b Vector of times of biomass measurements (in hours).
#' @param b.data Vector of biomass measurements.
#' @param m Maturation constant (if using fluorescent protein reporter).
#'   Defaults to log(2)/0.45 h^-1.
#' @param d Degradation rate of protein. Defaults to log(2)/0.5 h^-1 (Mateus & Avery
#'   2000).
#' @param beta Translation rate. Defaults to 1.
#' @param d.r mRNA degradation rate. Defaults to 1.
#' @param data.type Indicates type of input data, with numbers corresponding to
#'   number of ODEs necessary to model the reporter system. 1 = mRNA, 2 =
#'   protein other than fluorescent (i.e. luciferase, beta-galactosidase, etc.),
#'   3 = fluorescent protein (i.e. GFP, RFP, etc). For more information, see
#'   references. Defaults to 3.
#'
#' @return Returns a list \code{L} with \code{L$t} the time series input by the
#'   user and \code{L$promact} the calculated promoter activities at these
#'   times.
#'
#' @section Author:
#' Soumya Kannan, Technical University of Denmark, 2017.
#'
#' @references PUBLICATION CITATION.
#'
#' @import pspline
#' @importFrom fda create.bspline.basis fd smooth.monotone eval.monfd
#' @importFrom stats predict smooth.spline

#############################################################################

DIPA <- function(x, t.f, f.data, t.b, b.data, m  = log(2)/0.45,
                              d  = log(2)/0.5, beta = 1, d.r = 1, data.type = 3) {

  # Calculate average step size of data points
  step_size <- mean(diff(t.f))

  # Modify biomass measurements to have a minimum of 1 to avoid computational
  # issues regarding noise
  min_biomass <- min(b.data)
  if (min_biomass < 1) {
    diff_1 <- 1 - min_biomass
    b.data <- b.data + diff_1
  }

  # Calculate appropriate smoothing parameter for fluorescence data given step size
  n_spar <- 0.189 / sqrt(step_size) + 0.07
  if (n_spar > 0.6) {
    n_spar <- 0.6
  }

  if (step_size < 2) {
    # Fit spline model to fluorescence data
    f <- smooth.Pspline(x=t.f, y=f.data, norder=4, spar=n_spar)
    # Construct smoothed fluorescence
    f.t <- predict(object=f, x)
    f.t[which(f.t < 0)] <- 0 # Physical constraints
    # Calculate fluorescence derivatives
    dfdt <- predict(object=f, x, nderiv=1)  ## first derivative of GFP
    dfdt2 <- predict(object=f, x, nderiv=2) ## second derivative of GFP
    dfdt3 <- predict(object=f, x, nderiv=3) ## third derivative of GFP
  }
  else {
    # Fit spline model to fluorescence data
    f <- smooth.spline(x=t.f, y=f.data, spar=n_spar)
    # Construct smoothed fluorescence
    f.t <- predict(object=f, x)$y
    f.t[which(f.t < 0)] <- 0 # Physical constraints
    # Calculate fluorescence derivatives
    dfdt <- predict(object=f, x, nderiv=1)$y ## first derivative of GFP
    dfdt2 <- predict(object=f, x, nderiv=2)$y ## second derivative of GFP
    dfdt3 <- predict(object=f, x, nderiv=3)$y ## third derivative of GFP
  }

  # Calculate nbasis for monotone biomass smoothing
  nbasis <- length(t.b) / 2
  if (nbasis > 15) {
    nbasis <- 15
  }
  else if (nbasis < 7) {
    nbasis <- 7
  }

  # Fit monotonic spline to biomass data
  try_mt <- try({
    basis <- create.bspline.basis(rangeval=c(min(t.b) - 1, max(t.b) + 1), nbasis=nbasis)
    Wfd0 <- fd(coef=as.vector(matrix(0, basis$nbasis, 1)), basisobj=basis)
    g <- smooth.monotone(t.b, b.data, WfdParobj=Wfd0, dbglev=0)
  }, silent=TRUE)
  # If monotonic spline gives an error, use the usual one instead & construct
  # smoothed biomass
  if (class(try_mt) == "try-error") {
    g <- smooth.spline(x=t.b, y=b.data, cv=FALSE)
    if (g$spar > n_spar + 0.1) {
      g <- smooth.spline(x=t.b, y=b.data, spar=n_spar + 0.1)
    }
    g.t <- predict(object=g, x)$y
  }
  else {
    g.t <- g$beta[1] + g$beta[2] * eval.monfd(x, g$Wfdobj)
  }
  g.t[which(g.t < 0.001)] <- 0.001 # Physical constraints


  # Calculate promoter activity
  if (data.type == 1) {
    pAct <- (1 / g.t) * (dfdt + d.r * f.t)
  }
  else if (data.type == 2) {
    pAct <- (1 / (beta * g.t)) * (dfdt2 + (d.r + d) * dfdt + d.r*d * f.t)
  }
  else if (data.type == 3) {
    C_0 <- d.r*d*(d + m)
    C_1 <- (d + d.r)*(d + m) + d*d.r
    C_2 <- (2 * d + m + d.r)

    pAct <- (1 / (m * beta * g.t)) * (dfdt3 + C_2 * dfdt2 + C_1 * dfdt + C_0 * f.t)
  }

  # Return times and calculated promoter activities
  return( list(t=x, promact=pAct) )
}
