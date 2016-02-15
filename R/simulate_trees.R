#'@title .generate_newick
#'@description Generates a newick file for heterochronous trees. 
#'@param args is a list with coal_factor, sampling and coalescent times, events and indicator. This is the output of gen_INLA_args().
#'@param sample is a matrix or vector. The first column indicates the number of samples and the second column indicates the sampling times
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}
#'

.generate_newick<-function (args, sample) 
{
    ##This function is repeted from phylodyn (delete and add dependency to phylodyn)
  n <- sum(sample[, 1])
  labels <- paste(rep("t", n), seq(1, n, 1), rep("_", n), rep(sample[, 
                                                                     2], sample[, 1]), sep = "")
  tb <- sample[1, 1]
  s <- 0
  temp_labels <- labels[1:tb]
  temp_times <- rep(sample[1, 2], sample[1, 1])
  initial.row <- 2
  for (j in 2:length(args$event)) {
    if (args$event[j] == 1) {
      s <- args$s[j]
      ra <- sample(tb, 1)
      if (ra < tb) {
        new_label <- paste("(", temp_labels[ra], ":", 
                           s - temp_times[ra], ",", temp_labels[ra + 1], 
                           ":", s - temp_times[ra + 1], ")", sep = "")
        temp_labels[ra] <- new_label
        temp_labels <- temp_labels[-(ra + 1)]
        temp_times[ra] <- s
        temp_times <- temp_times[-(ra + 1)]
      }
      else {
        new_label <- paste("(", temp_labels[ra], ":", 
                           s - temp_times[ra], ",", temp_labels[1], ":", 
                           s - temp_times[1], ")", sep = "")
        temp_labels[1] <- new_label
        temp_labels <- temp_labels[-(ra)]
        temp_times[1] <- s
        temp_times <- temp_times[-(ra)]
      }
      tb <- tb - 1
    }
    else {
      s <- args$s[j]
      if (sample[initial.row, 1] == 1) {
        temp_labels <- c(temp_labels, labels[cumsum(sample[, 
                                                           1])[initial.row]])
        initial.row <- initial.row + 1
        tb <- tb + 1
        temp_times <- c(temp_times, s)
      }else{
          end<-cumsum(sample[,1])[initial.row]
          ini<-cumsum(sample[,1])[initial.row-1]+1
          for (k in ini:end){
              temp_labels<-c(temp_labels,labels[k])
              tb<-tb+1
              temp_times<-c(temp_times,s)
          }
          initial.row<-initial.row+1
      }
    }
  }
  out.tree <- ape::read.tree(text = paste(temp_labels, ";", sep = ""))
  return(list(newick = out.tree, labels = labels))
}

#'@title .coalgen_thinning_iso
#'@description Simulates coalescent times using thinning for isochronous sampling
#'@param sample is a matrix or vector. The first column indicates the number of samples and the second column indicates the sampling times
#'@param trajectory is a function with the population size function
#'@param upper is an upper value of inverse trajectory (optional)
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

.coalgen_thinning_iso<-function(sample, trajectory, upper = 25, ...) {
  s = sample[2]
  n <- sample[1]
  out <- rep(0, n - 1)
  time <- 0
  j <- n
  while (j > 1) {
    time_p <- time + rexp(1, upper * j * (j - 1) * 0.5)
    if (upper< trajectory(time_p, ...) ){
      upper<-trajectory(time_p, ...)*1.2
      time_p <- time + rexp(1, upper * j * (j - 1) * 0.5)
    }
    time<-time_p
    if (runif(1) <= trajectory(time, ...)/upper) {
      out[n - j + 1] <- time
      j <- j - 1
    }
  }
  return(list(intercoal_times = c(out[1], diff(out)), lineages = seq(n, 2, -1),upper=upper))
}


#'@title .coalgen_thinning_hetero
#'@description Simulates coalescent times using thinning for isochronous sampling
#'@param sample is a matrix or vector. The first column indicates the number of samples and the second column indicates the sampling times
#'@param trajectory is a function with the population size function
#'@param upper is an upper value of inverse trajectory (optional)
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}

.coalgen_thinning_hetero<-function (sample, trajectory, upper = 25, ...) {
  n_sampled<-sample[,1]
  s_times<-sample[,2]
  s = sample[1, 2]
  b <- sample[1, 1]
  n <- sum(sample[, 1]) - 1
  m <- n
  nsample <- nrow(sample)
  sample <- rbind(sample, c(0, 10 * max(sample, 2)))
  out <- rep(0, n)
  branches <- rep(0, n)
  i <- 1
  while (i < (nsample)) {
    if (b < 2) {
      b <- b + sample[i + 1, 1]
      s <- sample[i + 1, 2]
      i <- i + 1
    }
    E <- rexp(1, upper * b * (b - 1) * 0.5)
    if (upper< trajectory(E+s, ...) ){
      upper<-trajectory(E+s, ...)*1.2
      E <- rexp(1, upper * b * (b - 1) * 0.5)
    }
    
    if (runif(1) <= trajectory(E + s, ...)/upper) {
      if ((s + E) > sample[i + 1, 2]) {
        b <- b + sample[i + 1, 1]
        s <- sample[i + 1, 2]
        i <- i + 1
      }
      else {
        s <- s + E
        out[m - n + 1] <- s
        branches[m - n + 1] <- b
        n <- n - 1
        b <- b - 1
      }
    }
    else {
      s <- s + E
    }
  }
  while (b > 1) {
    E <- rexp(1, upper * b * (b - 1) * 0.5)
    if (runif(1) <= trajectory(E + s, ...)/upper) {
      s <- s + E
      out[m - n + 1] <- s
      branches[m - n + 1] <- b
      n <- n - 1
      b <- b - 1
    }
    else {
      s <- s + E
    }
  }
  intercoal_times = c(out[1], diff(out))
  lineages = branches
  coal_times<-cumsum(intercoal_times)
  args<-gen_INLA_args(coal_times,s_times,n_sampled)
  out<-.generate_newick(args,cbind(n_sampled,s_times))$newick
  return(out)
}

#'@title gen_INLA_args
#'@description Used for heterochronous sampling. It has nothing to do with INLA
#'@param coal_times indicates the coalescent times
#'@param s_times indicates the sampling times
#'@param n_sampled the number of samples 
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}
#'
gen_INLA_args<-function (coal_times, s_times, n_sampled) {
  n = length(coal_times) + 1
  data = matrix(0, nrow = n - 1, ncol = 2)
  data[, 1] = coal_times
  s_times = c(s_times, max(data[, 1]) + 1)
  data[1, 2] = sum(n_sampled[s_times <= data[1, 1]])
  tt = length(s_times[s_times <= data[1, 1]]) + 1
  for (j in 2:nrow(data)) {
    if (data[j, 1] < s_times[tt]) {
      data[j, 2] = data[j - 1, 2] - 1
    }
    else {
      data[j, 2] = data[j - 1, 2] - 1 + sum(n_sampled[s_times > 
                                                        data[j - 1, 1] & s_times <= data[j, 1]])
      tt = length(s_times[s_times <= data[j, 1]]) + 1
    }
  }
  s = unique(sort(c(data[, 1], s_times[1:length(s_times) - 1])))
  event1 = sort(c(data[, 1], s_times[1:length(s_times) - 1]), index.return = TRUE)$ix
  n = nrow(data) + 1
  l = length(s)
  event = rep(0, l)
  event[event1 < n] = 1
  y = diff(s)
  coal.factor = rep(0, l - 1)
  t = rep(0, l - 1)
  indicator = cumsum(n_sampled[s_times < data[1, 1]])
  indicator = c(indicator, indicator[length(indicator)] - 1)
  ini = length(indicator) + 1
  for (k in ini:(l - 1)) {
    j = data[data[, 1] < s[k + 1] & data[, 1] >= s[k], 2]
    if (length(j) == 0) {
      indicator[k] = indicator[k - 1] + sum(n_sampled[s_times < 
                                                        s[k + 1] & s_times >= s[k]])
    }
    if (length(j) > 0) {
      indicator[k] = j - 1 + sum(n_sampled[s_times < s[k +1] & s_times >= s[k]])
    }
  }
  coal_factor = indicator * (indicator - 1)/2
  return(list(coal_factor = coal_factor, s = s, event = event, indicator = indicator))
}

##This is the main function

#'@title simulate.tree
#'@description Simulates genealogies under the coalescent
#'@param n Number of samples
#'@param N Number of genealogies or trees
#'@param sampling It can take two values: iso or hetero. If the sampling scheme is (iso) it generates an isochronous genealogy (all samples at time 0). Otherwise, it assumes it is heterochronous sampling.If hetero, you need to specify samples
#'@param args String of characters for ms
#'@param Ne The standard is set to 1. It takes any function as input for the demographic model
#'@param max (optional) for thinning simulator. It is an upper bound on 1/Ne
#'@param simulator It can either be null, ms or thinning. If null, it uses ape rcoal.
#'@param sample It takes a matrix of two columns. The first column inciates the number of samples and the second column the sampling times
#'@author Julia Palacios \email{julia.pal.r@@gmail.com}
#'
simulate.tree<-function(n=10,N=1,sampling="iso",args="-T -G 0.1",Ne=1,max=1,simulator=NULL,sample=NULL, ...){
  if (is.null(simulator)) {
    #Generates samples using rcoal with Ne=1. If this is a mistake, specify your simulator (ms,thinning,standard)
    out<-replicate(N,rcoal(n),simplify=FALSE)
    class(out)<-"multiPhylo"
    return(list(out=out,description="Simulation from a constant population size Ne=1"))
  }
  if (simulator=="ms"){
    library("phyclust")
    if (is.null(args)){args<-"-T"}
    out<-read.tree(text=paste(rep(ms(nsam = n, opts = args)[3],N),sep="\n"))  
    return(list(out=out,description=paste("Simulation using ms with args:",args,sep=" ")))
  }
  if (simulator=="thinning"){
    if (is.function(Ne)){
      ##A function is needed for thinning
      fun_inv<-function(t,...){
        return(1/Ne(t,...))
      }
    }
    if (is.null(Ne)){
      fun_inv<-function(t){
        return(1)
      }
    }
     
    #upper is an upper bound on fun_inv
    if (sampling=="iso"){
      
      if (is.null(sample)){
        sample<-c(n,0)
      }
      out<-replicate(N,rcoal(n,br=.coalgen_thinning_iso(sample,fun_inv,max)$intercoal_times),simplify=FALSE)
      class(out)<-"multiPhylo"
      
      return(list(out=out,description="isochronous simulation using thinning with trajectory provided"))
      
      
    }
    
    else{
      #        #sampling is heterochronous
      if (is.null(sample)){ #it is actually isochronous
        sample<-c(n,0)
        out<-replicate(N,rcoal(n,br=.coalgen_thinning_iso(sample,fun_inv,max)$intercoal_times),simplify=FALSE)
        class(out)<-"multiPhylo"
        
        return(list(out=out,description="isochronous simulation using thinning with trajectory provided"))
      }else{
        #It is indeed heterochronous
        out<-replicate(N,.coalgen_thinning_hetero(sample, fun_inv,max),simplify=FALSE)     
        class(out)<-"multiPhylo" 
        return(list(out=out,description="Simulation with thinning for heterochronous coalescent with a specific function"))
        #   
      }
    }
  }
  
}



#'@title .branching.sampling.times
#'@description Used for heterochronous genealogies
#'@param phy takes a phylo object
#'@author Julia Palacios and Vladimir Minin \email{julia.pal.r@@gmail.com}
.branching.sampling.times <- function(phy){
  phy = new2old.phylo(phy)
  if (class(phy) != "phylo")
    stop("object \"phy\" is not of class \"phylo\"")
  tmp <- as.numeric(phy$edge)
  nb.tip <- max(tmp)
  nb.node <- -min(tmp)
  xx <- as.numeric(rep(NA, nb.tip + nb.node))
  names(xx) <- as.character(c(-(1:nb.node), 1:nb.tip))
  xx["-1"] <- 0
  for (i in 2:length(xx)) {
    nod <- names(xx[i])
    ind <- which(phy$edge[, 2] == nod)
    base <- phy$edge[ind, 1]
    xx[i] <- xx[base] + phy$edge.length[ind]
  }
  depth <- max(xx)
  branching.sampling.times <- depth - xx
  return(branching.sampling.times)
}

#'@title heterochronous.gp.stat
#'@description Generates sufficient statistcs 
#'@param phy takes a phylo object
#'@author Julia Palacios and Vladimir Minin \email{julia.pal.r@@gmail.com}
heterochronous.gp.stat <- function(phy){
  b.s.times = .branching.sampling.times(phy)
  int.ind = which(as.numeric(names(b.s.times)) < 0)
  tip.ind = which(as.numeric(names(b.s.times)) > 0)
  num.tips = length(tip.ind)
  num.coal.events = length(int.ind)
  sampl.suf.stat = rep(NA, num.coal.events)
  coal.interval = rep(NA, num.coal.events)
  coal.lineages = rep(NA, num.coal.events)
  sorted.coal.times = sort(b.s.times[int.ind])
  names(sorted.coal.times) = NULL
  #unique.sampling.times = sort(unique(b.s.times[tip.ind]))
  sampling.times = sort((b.s.times[tip.ind]))
  for (i in 2:length(sampling.times)){
    if ((sampling.times[i]-sampling.times[i-1])<0.1){
      sampling.times[i]<-sampling.times[i-1]}
  }
  unique.sampling.times<-unique(sampling.times)
  sampled.lineages = NULL
  for (sample.time in unique.sampling.times){
    sampled.lineages = c(sampled.lineages,
                         sum(sampling.times == sample.time))  
  }
  return(list(coal.times=sorted.coal.times, sample.times = unique.sampling.times, sampled.lineages=sampled.lineages))  
}
