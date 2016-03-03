#' Multi Loci Bayesian nonparametric phylodynamic reconstruction.
#' 
#' @param data \code{multiPhylo} object or list containing vectors of coalescent 
#'   times \code{coal_times}, sampling times \code{samp_times}, and number 
#'   sampled per sampling time \code{n_sampled}.
#' @param lengthout numeric specifying number of grid points.
#' @param pref logical. Should the preferential sampling model be used?
#' @param prec_alpha,prec_beta numerics specifying gamma prior for precision 
#'   \eqn{\tau}.
#' @param beta1_prec numeric specifying precision for normal prior on 
#'   \eqn{\beta_1}.
#' @param simplify logical whether to fully bucket all Poisson points.
#' @param derivative logical whether to calculate estimates of the 
#'   log-derivative.
#' @param forward logical whether to use the finite difference approximations of
#'   the log-derivative as a forward or backward derivative.
#'   
#' @return Phylodynamic reconstruction of effective population size from multiple independent genealogies 
#' at grid points. \code{result} contains the INLA output, \code{data} contains the 
#'   information passed to INLA, \code{grid} contains the grid end points, 
#'   \code{x} contains the grid point centers, \code{effpop} contains a vector 
#'   of the posterior median effective population size estimates, 
#'   \code{effpop025} and \code{effpop975} contain the 2.5th and 97.5th 
#'   posterior percentiles, \code{summary} contains a data.frame of the 
#'   estimates, and \code{derivative} (if \code{derivative = TRUE}) contains a
#'   data.frame summarizing the log-derivative.
#' @export
#' 
#' @examples
#' res = BNPR_Multiple(simulate.tree(n=10,N=2,simulator="ms",args="-T -G 0.1")$out)
#' phylo::plot_BNPR(res)
BNPR_Multiple <- function(data, lengthout = 100, pref=FALSE, prec_alpha=0.01,
                 prec_beta=0.01, beta1_prec = 0.001, simplify = TRUE,
                 derivative = FALSE, forward = TRUE)
{
  if (class(data) != "multiPhylo")
  {
    return(phylodyn::BNPR(data,lengthout,pref,prec_alpha,prec_beta,beta1_prec,simplify,derivative,forward))
  }
  else {
    
    ntrees<-length(data)
    phy <- phylodyn::summarize_phylo(data[[1]])  
    range.upp<-max(phy$coal_times)
    range.low<-min(phy$samp_times)
    
    for (j in 2:ntrees){
    phy <- phylodyn::summarize_phylo(data[[j]]) 
    range.upp<-max(range.upp,max(phy$coal_times))
    range.low<-min(range.low,min(phy$samp_times))
    }
    grid <- seq(range.low, range.upp, length.out = lengthout+1)
    
    #//replace by multiple
    
    result <- infer_coal_multiple(samp_times = phy$samp_times, coal_times = phy$coal_times,
                            n_sampled = phy$n_sampled, grid=grid,
                            prec_alpha = prec_alpha, prec_beta = prec_beta,
                            beta1_prec = beta1_prec, use_samp = pref,
                            simplify = simplify, derivative = derivative)
  
  result$samp_times <- phy$samp_times
  result$n_sampled  <- phy$n_sampled
  result$coal_times <- phy$coal_times
  
  result$effpop    <- exp(-result$result$summary.random$time$`0.5quant`)
  result$effpop975 <- exp(-result$result$summary.random$time$`0.025quant`)
  result$effpop025 <- exp(-result$result$summary.random$time$`0.975quant`)
  
  result$summary <- with(result$result$summary.random$time,
                         data.frame(time = ID, mean = exp(-mean),
                                    sd = sd * exp(-mean),
                                    quant0.025 = exp(-`0.975quant`),
                                    quant0.5 = exp(-`0.5quant`),
                                    quant0.975 = exp(-`0.025quant`)))
  
  if (derivative)
  {
    if (forward)
      ind <- c(1:(lengthout-1), (lengthout-1))
    else
      ind <- c(1, 1:(lengthout-1))
    
    result$derivative <- with(result$result$summary.lincomb,
                              data.frame(time = result$x, mean = -mean[ind], sd = sd[ind],
                                         quant0.025 = -`0.975quant`[ind],
                                         quant0.5   = -`0.5quant`[ind],
                                         quant0.975 = -`0.025quant`[ind]))
  }
  
  
  return(result)
}

  infer_coal_multiple <- function(data=data, grid,
                         prec_alpha = 0.01, prec_beta = 0.01, simplify = FALSE,
                         derivative = FALSE)
  {
    j<-1
    phy <- phylodyn::summarize_phylo(data[[j]]) 
    if (min(phy$coal_times) < min(phy$samp_times))
      stop("First coalescent time occurs before first sampling time")
    
    if (max(phy$samp_times) > max(phy$coal_times))
      stop("Last sampling time occurs after last coalescent time")
    
    if (is.null(phy$n_sampled))
      phy$n_sampled <- rep(1, length(samp_times))
    
    
    coal_data <- phylodyn::coal_stats(grid = grid, samp_times = phy$samp_times, n_sampled = phy$n_sampled,
                                      coal_times = phy$coal_times)
    #check whether having 0 is a good idea
    coal_data<-coal_data[coal_data[,4]!=-100L,]
    coal_data$tree<-rep(1,nrow(coal_data))
    
    for (j in 2:ntrees){
      phy <- phylodyn::summarize_phylo(data[[j]]) 
      if (min(phy$coal_times) < min(phy$samp_times))
        stop("First coalescent time occurs before first sampling time")
      
      if (max(phy$samp_times) > max(phy$coal_times))
        stop("Last sampling time occurs after last coalescent time")
      
      if (is.null(phy$n_sampled))
        phy$n_sampled <- rep(1, length(samp_times))
      
      
    coal_data_temp <- phylodyn::coal_stats(grid = grid, samp_times = phy$samp_times, n_sampled = phy$n_sampled,
                            coal_times = phy$coal_times)
    #check whether having 0 is a good idea
    coal_data_temp<-coal_data_temp[coal_data_temp[,4]!=-100L,]
    coal_data_temp$tree<-rep(j,nrow(coal_data_temp))
    coal_data<-rbind(coal_data,coal_data_temp)
    }
    
    
    if (simplify)
      coal_data <- with(coal_data, phylodyn::condense_stats(time = time, event = event, E=E))
    
    data <- with(coal_data, data.frame(y = event, time = time, E_log = E_log))
    hyper <- list(prec = list(param = c(prec_alpha, prec_beta)))
    formula <- y ~ -1 + f(time, model="rw1", hyper = hyper, constr = FALSE)
    
    if (derivative)
    {
      Imat <- diag(lengthout)
      A <- head(Imat, -1) - tail(Imat, -1)
      field <- grid[-1] - diff(grid)/2
      A <- diag(1/diff(field)) %*% A
      A[A == 0] <- NA
      
      lc_many <- INLA::inla.make.lincombs(time = A)
      
      mod <- INLA::inla(formula, family = "poisson", data = data, lincomb = lc_many,
                        control.predictor = list(compute=TRUE),
                        control.inla = list(lincomb.derived.only=FALSE))
    }
    else
    {
      mod <- INLA::inla(formula, family = "poisson", data = data, offset = data$E_log,
                        control.predictor = list(compute=TRUE))
    }
    
    return(list(result = mod, data = data, grid = grid, x = coal_data$time))
  }
#'@title calculate.moller.hetero
#'@description Approximates the posterior distribution of Ne from a single genealogy at a regular grid of points using INLA package
#'@param coal.factor is a vector with coalescent times in increasing order
#'@param s sampling times and coalescent times in increasing order
#'@param event and indicator vector with 1 for coalescent times and 0 for sampling times that correspond to the s vector
#'@param lengthout number of grid points
#'@param prec_alpha alpha gamma hyperparameter of precision parameter for GP prior
#'@param prec_beta  beta gamma hyperparameter of precision parameter for GP prior
#'@param E.log.zero internal
#'@param alpha TODO
#'@param beta TODO
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

calculate.moller.hetero<-function (coal.factor, s, event, lengthout, prec_alpha = 0.01,
          prec_beta = 0.01, E.log.zero = -100, alpha = NULL, beta = NULL)
{

  if (prec_alpha == 0.01 & prec_beta == 0.01 & !is.null(alpha) &
        !is.null(beta)) {
    prec_alpha = alpha
    prec_beta = beta
  }
  grid <- seq(0, max(s), length.out = lengthout + 1)
  u <- diff(grid)
  field <- grid[-1] - u/2
  sgrid <- grid
  event_new <- 0
  time <- 0
  where <- 1
  E.factor <- 0
  for (j in 1:lengthout) {
    count <- sum(s > sgrid[j] & s <= sgrid[j + 1])
    if (count > 1) {
      points <- s[s > sgrid[j] & s <= sgrid[j + 1]]
      u <- diff(c(sgrid[j], points))
      event_new <- c(event_new, event[(where):(where +
                                                 count - 1)])
      time <- c(time, rep(field[j], count))
      E.factor <- c(E.factor, coal.factor[where:(where +
                                                   count - 1)] * u)
      where <- where + count
      if (max(points) < sgrid[j + 1]) {
        event_new <- c(event_new, 0)
        time <- c(time, field[j])
        E.factor <- c(E.factor, coal.factor[where] *
                        (sgrid[j + 1] - max(points)))
      }
    }
    if (count == 1) {
      event_new <- c(event_new, event[where])
      points <- s[s > sgrid[j] & s <= sgrid[j + 1]]
      if (points == sgrid[j + 1]) {
        E.factor <- c(E.factor, coal.factor[where] *
                        (sgrid[j + 1] - sgrid[j]))
        time <- c(time, field[j])
        where <- where + 1
      }
      else {
        event_new <- c(event_new, 0)
        E.factor <- c(E.factor, coal.factor[where] *
                        (points - sgrid[j]))
        E.factor <- c(E.factor, coal.factor[where + 1] *
                        (sgrid[j + 1] - points))
        time <- c(time, rep(field[j], 2))
        where <- where + 1
      }
    }
    if (count == 0) {
      event_new <- c(event_new, 0)
      E.factor <- c(E.factor, coal.factor[where] * (sgrid[j +
                                                            1] - sgrid[j]))
      time <- c(time, field[j])
    }
  }
  time2 <- time
  event_new2 <- event_new
  E.factor2 <- E.factor
  for (j in 1:lengthout) {
    count <- sum(time2 == field[j])
    if (count > 1) {
      indic <- seq(1:length(event_new2))[time2 == field[j]]
      if (sum(event_new2[indic]) == 0) {
        event_new2 <- event_new2[-indic[-1]]
        time2 <- time2[-indic[-1]]
        temp <- sum(E.factor2[indic])
        E.factor2[indic[1]] <- temp
        E.factor2 <- E.factor2[-indic[-1]]
      }
    }
  }
  E.factor2.log = log(E.factor2)
  E.factor2.log[E.factor2 == 0] = E.log.zero
  data <- list(y = event_new2[-1], event = event_new2[-1],
               time = time2[-1], E = E.factor2.log[-1])
  formula <- y ~ -1 + f(time, model = "rw1", hyper = list(prec = list(param = c(prec_alpha,
                                                                                prec_beta))), constr = FALSE)
  mod4 <- inla(formula, family = "poisson", data = data, offset = E,
               control.predictor = list(compute = TRUE))
  return(list(result = mod4, grid = grid, data = data, E = E.factor2.log))
}

#'@title .calculate.moller1
#'@description Approximates the posterior distribution of Ne from a single genealogy at a regular grid of points using INLA package
#'@param tree is a phylo object with a single tree
#'@param lengthout number of grid points
#'@param L the length for the definition of the grid
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

.calculate.moller1<-function(tree,lengthout,L=1){
  ci<-coalescent.intervals(tree)
  data1<-cbind(ci$interval.length,ci$lineages)
  s<-cumsum(data1[,1])
  coal.factor<-data1[,2]*(data1[,2]-1)*.5
  length.min<-min(data1[,1])
  #maxval<-sum(data1[,1])
  grid<-seq(0,L,length.out=lengthout+1)
  #grid<-c(grid[grid<=maxval],maxval)
  u<-diff(grid)
  field<-grid[-1]-u/2
  sgrid<-grid
  event<-0
  time<-0
  E.factor<-0
  where<-1
  for (j in 1:lengthout){
    count<-sum(s>sgrid[j] & s<=sgrid[j+1])
    if (count>1){
      points<-s[s>sgrid[j] & s<=sgrid[j+1]]
      u<-diff(c(sgrid[j],points))
      event<-c(event,rep(1,count))
      time<-c(time,rep(field[j],count))
      E.factor<-c(E.factor,coal.factor[where:(where+count-1)]*u)
      where<-where+count
      if (max(points)<sgrid[j+1]){
        event<-c(event,0)
        time<-c(time,field[j])
        E.factor<-c(E.factor,coal.factor[where]*(sgrid[j+1]-max(points)))
      }
    }
    if (count==1){
      event<-c(event,1)
      points<-s[s>sgrid[j] & s<=sgrid[j+1]]
      if (points==sgrid[j+1]){
        E.factor<-c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
        time<-c(time,field[j])
        where<-where+1
      }else {
        event<-c(event,0)
        E.factor<-c(E.factor,coal.factor[where]*(points-sgrid[j]))
        E.factor<-c(E.factor,coal.factor[where+1]*(sgrid[j+1]-points))
        time<-c(time,rep(field[j],2))
        where<-where+1}
    }
    if (count==0){
      event<-c(event,0)
      E.factor<-c(E.factor,coal.factor[where]*(sgrid[j+1]-sgrid[j]))
      time<-c(time,field[j])
    }

  }
  ##Fixing for when there are no observations
  combine<-unique(sort(c(grid[-1],s)))
  find2<-max(seq(1,length(combine))[combine<=max(s)])
  event<-event[-1]
  time<-time[-1]
  E.factor<-E.factor[-1]
  grid<-c(grid[grid<=max(s)],max(s))
  data<-list(y=event[1:find2],event=event[1:find2],time=time[1:find2],E=log(E.factor[1:find2]))

  #  data<-list(y=event[-1],event=event[-1],time=time[-1],E=log(E.factor[-1]))
  formula<-y~-1+f(time,model="rw1",hyper=list(prec = list(param = c(.001, .001))),constr=FALSE)
  mod.moller.constant<-inla(formula,family="poisson",data=data,offset=E,control.predictor=list(compute=TRUE))

  return(stan_output(list(result=mod.moller.constant,grid=grid)))
}

#'@title calculate.moller
#'@description Approximates the posterior distribution of Ne from a single genealogy at a regular grid of points using INLA package
#'@param tree is a phylo object with a single tree
#'@param lengthout number of grid points
#'@param L the length for the definition of the grid
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

calculate.moller<-function(tree,lengthout,L=1){
 if (class(tree)=="phylo") {return(.calculate.moller1(tree,lengthout,L))}
 else {
   L<-coalescent.intervals(tree[[1]])$total.depth
   for (j in 2:length(tree)){L<-min(L,coalescent.intervals(tree[[j]])$total.depth)}
   result<-matrix(NA,nrow=lengthout,ncol=length(tree)+1)
   result[,1:2]<-.calculate.moller1(tree[[1]],lengthout,L)[,1:2]
   for (j in (2:length(tree))){
     result[,j+1]<-.calculate.moller1(tree[[j]],lengthout,L)[,2]
   }
   return(result)
 }
}

#'@title plot_INLA
#'@description Plots the output from the inla functions for Ne
#'@param INLA_out Otput from the inla functions
#'@param traj the true trajectory
#'@param xlim the limits of the x-axis, as standard in plot()
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

plot_INLA = function(INLA_out, traj=NULL, xlim=NULL, ...)
{
  mod = INLA_out$result$summary.random$time

  grid = mod$"ID"
  if (is.null(xlim))
  {
    xlim=c(max(grid),0)
  }
  plot(grid,exp(-mod$"0.5quant"),type="l",lwd=2.5,col="blue",log="y",
       xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
       xlim=xlim,ylim=c(min(exp(-mod$"0.975quant"[grid > min(xlim) & grid < max(xlim)])),
                        max(exp(-mod$"0.025quant"[grid > min(xlim) & grid < max(xlim)]))), ...)
  lines(grid,exp(-mod$"0.975quant"),lwd=2.5,col="blue",lty=2)
  lines(grid,exp(-mod$"0.025quant"),lwd=2.5,col="blue",lty=2)
  if (!is.null(traj))
    lines(grid, traj(grid))
}


#'@title plot_INLA2
#'@description Plots a general output with 4 columns: time, median, lower and upper bound
#'@param result Matrix with four columns
#'@param traj the true trajectory
#'@param xlim (optional)
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

plot_INLA2<-function(result, traj=NULL, xlim=NULL, ...){
  grid = result[,1]
  if (is.null(xlim)) {
    xlim = c(max(grid), 0)
  }

  plot(grid,result[,2],type="l",lwd=2.5,col="blue",log="y",
       xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
       xlim=xlim, ylim=c(min(result[grid > min(xlim) & grid < max(xlim),4]),
                         max(result[grid > min(xlim) & grid < max(xlim),3])))
  lines(grid,result[,3],lwd=2.5,col="blue",lty=2)
  lines(grid,result[,4],lwd=2.5,col="blue",lty=2)
  if (!is.null(traj)) lines(grid, traj(grid))
}


#'@title plot_INLA3
#'@description Plots bunch of trees
#'@param result Matrix with many columns
#'@param traj the true trajectory
#'@param xlim (optional)
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}
plot_INLA3<-function(result, traj=NULL, xlim=NULL, ...){
  grid = result[,1]
  if (is.null(xlim)) {
    xlim = c(max(grid), 0)
  }
  result2<-matrix(NA,nrow=nrow(result),ncol=2)
  result2[,2]<-apply(result[,2:ncol(result)],1,mean)
  plot(grid,result2[,2],type="l",lwd=2.5,col="blue",log="y",
       xlab="Time (past to present)",ylab="Scaled Effective Pop. Size",
       xlim=xlim, ylim=c(min(result[,2:ncol(result)]), max(result[,2:ncol(result)])))
  for (j in 2:(ncol(result))){
  lines(grid,result[,j],lwd=1.0,col="gray")}
  if (!is.null(traj)) lines(grid, traj(grid))
}

stan_output<-function(INLA_out){
  mod = INLA_out$result$summary.random$time
  grid = mod$"ID"

  results<-matrix(NA,nrow=length(grid),ncol=4)
  results[,1]<-grid
  results[,4]<-exp(-mod$"0.975quant")
  results[,3]<-exp(-mod$"0.025quant")
  results[,2]<-exp(-mod$"0.5quant")
  return(results)
}
