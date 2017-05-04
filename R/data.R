#' Constant fluorescence during stationary phase
#'
#' A simulation of biomass and fluorescence data with fluorescence observed 
#' during stationary phase. Variables are as follows:
#'
#' @format A data frame with 30100 rows and 4 variables:
#' \describe{
#'    \item{sim_num}{simulation number (1--100)}
#'    \item{time}{time in hours (0--30)}
#'    \item{biomass}{simulated biomass measurement ()}
#'    \item{fluorescence}{simulated fluorescence measurement()}
#' }
"sim_1"

#' Pulse of fluorescence on entry to stationary phase
#' 
#' A simulation of biomass and fluorescence data with a ulse of fluorescence 
#' on entrance to stationary phase. Variables are as follows:
#' 
#' @format A data frame with 30100 rows and 4 variables:
#' \describe{
#'    \item{sim_num}{simulation number (1--100)}
#'    \item{time}{time in hours (0--30)}
#'    \item{biomass}{simulated biomass measurement ()}
#'    \item{fluorescence}{simulated fluorescence measurement()}
#' }
"sim_2"

#' Constant fluorescence during lag phase
#' 
#' A simulation of biomass and fluorescence with constant fluorescence in lag phase 
#' that decreases exponentially during exponential phase. Variables are as follows:
#' 
#' @format A data frame with 30100 rows and 4 variables:
#' \describe{
#'    \item{sim_num}{simulation number (1--100)}
#'    \item{time}{time in hours (0--30)}
#'    \item{biomass}{simulated biomass measurement ()}
#'    \item{fluorescence}{simulated fluorescence measurement()}
#' }
"sim_3"

#' Step function of promoter activity with different data step sizes
#' 
#' Simulations of biomass and fluorescence data for a step function of promoter activity
#' occurring from 12-16 h. Includes 100 simulations each for 4 different step sizes of 
#' data. Variables are as follows
#' 
#' @format A data frame with 41400 rows and 5 variables:
#' \describe{
#'    \item{step_size}{Step size/distance between data points in time (0.1, 0.5, 1, 1.5)}
#'    \item{sim_num}{simulation number (1--100)}
#'    \item{time}{time in hours (0--30)}
#'    \item{biomass}{simulated biomass measurement ()}
#'    \item{fluorescence}{simulated fluorescence measurement()}
#' }
"step_size_sims"

#' GFP driven by pXylA
#' 
#' A dataset containing biomass and fluorescence data for E. coli strains with plasmid-based
#' expression of GFP driven by a xylose-inducible pXylA promoter. 3 replicates of induced
#' and 3 of uninduced cultures are included. Note: RAW DATA - does not include background
#' correction described in publication.
#' 
#' @format A data frame with 1134 rows and 5 variables:
#' \describe{
#'    \item{sample}{sample number (1--6)}
#'    \item{time}{time in hours ()}
#'    \item{measurement}{biomass/fluorescence measurements}
#'    \item{meas_type}{whether or not the measurement is of biomass or fluorescence ("biomass", "fluorescence")}
#'    \item{induced}{whether or not the sample is induced or control (1 (induced), 0 (uninduced))}
#' } 
"pXylA"