#' @title .George
#' @description For Thibaut
#'
#' @examples .George(1)
#'

.George <- function(x)
{
  print("This function doesn't do anything")
}

#' @title getNodeAges
#' @description Returns a vector with the depths, or heights of the nodes in the tree
#' @description Depends on a hidden function nodeHeight, from phangorn.
#' @description Values by default are backwards from the youngest tip which is zero. They can be forward relative to the root if from_past=TRUE.
#'
#' @author George Shirreff <georgeshirreff@@gmail.com>
#'
#' @return For a tree with n tips, the vector is <tip1, tip2...tipn, root, internal_node2, internal_node3>
#'
#' @param x A phylo object
#' @param from_past If FALSE (default) the node ages increase backwards in time, with the most recent tip being zero. If TRUE the root node is zero and the node ages increase towards the tips of the tree.
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples Phy2Sky(trees)

getNodeAges<-function(x,from_past=F)
{
  # require(phangorn)
  nodeAges<-phangorn:::nodeHeight(x)
  names(nodeAges)<-c(x$tip.label,"root",rep(NA,length(nodeAges)-length(x$tip.label)-1))

  if(from_past)
  {
    return(max(nodeAges)-nodeAges)
  }
  return(nodeAges)
}

#' @title Phy2Sky
#' @description Converts a (multi)Phylo object to a series of skylines
#' @description For a set of unrooted trees, remove the burnin trees and root the rest at
#' @description the root.node. Limit the total number of output trees to be less than max.trees
#' @description Robust to trees with tips sampled at different times
#' @description Can be directly applied to BEAST or MrBayes output, or files can be directly parsed by argument
#'
#' @param trees Either a character string specifying the path to the tree file, or a multiPhylo object
#' @param output_type list (a list of skyline objects)
#' @param output_type matrices (two matrices, each with a column for each tree and a row for each coalescent event in the tree. One matrix specifies the population size and one specifies the time at the end of the interval)
#' @param output_type master (a master table with the rows being every unique time point across all trees, with time specified in the row name, and each column is a different tree)
#' @param output_type conf.int (a matrix with each row being each unique event time and the columns being the 50th, 2.5th and 97.5th percentiles respectively)
#' @param output_type conf.int.plot (outputs the conf.int table but also plots them)
#' @param file_type: the format of the tree file
#' @param package: what package has been used to create the param file (BEAST, MrBayes)
#' @param outgroup: specify if the tree should be re-rooted at this node and then this node removed.
#' @param param.file.name: if there is a parameter output file, specify it here so the clock rate can be fetched
#' @param burninfrac: what proportion of trees should be discarded from the multiPhylo object
#' @param max.trees: the maximum number of trees which will be sampled. If more exist, they are taken from the total at regular intervals.
#' @param fixToDate: if true, the times in the skyline will be linked directly to dates given in the tree name.
#' @param DateFUN: a function for recovering dates from the tip labels in the tree e.g. function(label) gsub(".*_([.0-9]+)","\\1",label)
#' @param timeForward: if TRUE the times of coalescent events wil be specified forwards, otherwise backwards.
#' @param scaling: a scaling factor for the branch lengths, if the branch lengths are not on the same scale as time.
#' @param plot_type: choose between "step", where adjacent points on the skyline are joined by a vertical and a horizontal line, and "linear", where adjacent points on the skyline are joined by a single sloping straight line.
#'
#' @author Lucy Mengqi Li <mengqi.li09@@imperial.ac.uk>, adapted by George Shirreff <georgeshirreff@@gmail.com>
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples Phy2Sky(trees)
#'

#skyline_list<-Phy2Sky(trees=phylo,output_type="list",burninfrac=burninfrac,max.trees=max.trees,youngestTip=youngestTip,timeForwards=T)
#conf.int_ES<-Phy2Sky(trees=trees_ES,epsilon=epsilon,output_type = "conf.int.plot",timeForwards=T,fixToDate=T,Date_FUN=Date_FUN,max.trees = 100)
#skyline_conf.int<-Phy2Sky(trees=trees,burninfrac=0,max.trees=1000,Date_FUN=Date_FUN,output_type = "conf.int.plot",fixToDate = T,timeForwards=T)

Phy2Sky <- function (trees, output_type=c("list","matrices","master","conf.int","conf.int.plot") ,file_type="nex", package=NULL, outgroup=NULL,
                             param.file.name=NULL, burninfrac=0, max.trees=1000,fixToDate=F,Date_FUN=Date_FUN,timeForwards=F,scaling=1,plot_type="step",epsilon=0,x_lab="Time",...) {

  require(ape)
  if(fixToDate){
    if(missing(Date_FUN)) stop("fixing to a date requires a date function")
    if(!missing(timeForwards)) print("timeForwards is ignored if fixToDate=TRUE")
  }

  if(is.character(trees))
  {
    if (file_type=="nex")
      trs <- read.nexus(file=trees)
    else if (file_type=="nwk")
      trs <- read.tree(trees)
  } else {
    trs<-trees
  }

  if(class(trs)=="phylo")
  {
    tmp<-list()
    tmp[[1]]=trs
    class(tmp)<-"multiPhylo"
    trs<-tmp
  }

  index.range <- c(floor(burninfrac*length(trs)+1), length(trs))
  if (diff(index.range) >= max.trees) {
    index.seq <- round(seq(floor(burninfrac*length(trs)+1), length(trs),
                     length.out=max.trees))
  } else {
    index.seq <- floor(burninfrac*length(trs)+1):length(trs)
  }
  if(is.null(outgroup)) {
    rooted.trs <- trs[index.seq]
  } else {
    rooted.trs <- lapply(index.seq, function (i){
      rooted <- root(trs[[i]], outgroup=outgroup, resolve.root=TRUE)  # MrBayes trees are already rooted
      drop.tip(rooted, outgroup, FALSE)
      #drop.tip(trs[[i]], root.node, FALSE)
    })
  }

#   if(!is.null(param.file.name))
#   {
#     pars <- read.table(param.file.name, header=TRUE, skip=1)
#     if(package=="BEAST")
#     {
#       clock.rates <- pars$clock.rate[index.seq]
#     } else if (package=="MrBayes"){
#       clock.rates <- pars$Clockrate[index.seq]
#     } else {
#       stop("parameter file package should be specified")
#     }
#   } else {
#     clock.rates=rep(NA,length(rooted.trs))
#   }
#

  #i=1
  #plot(x)
  #getNodeAges(x,from_past=F)
  #plot(sort(getNodeAges(x,from_past=T)))
# i=1
  skyline_list<-lapply(1:length(rooted.trs), function (i) {
    x <- rooted.trs[[i]]
    if(class(x)!="phylo") stop(paste("object",i,"is not a phylo object"))
    if(!is.binary.tree(x)) stop(paste("tree",i,"is not binary, with unknown consequences"))
    x_dated=x
    class(x_dated) <- "datedPhylo"
    #class(x) <- "phylo"
    sk <- skyline(x_dated)


#     getNodeAges(rooted.trs[[i]])
    if(fixToDate){
      n=length(x$tip.label)
      node_ages=getNodeAges(x)

      tip_time <- rev(range(node_ages[1:n]))

#       tree_time <- c(max(sk$time), #the tree age of the root
#         0) #the tree age of the youngest tip

      real_time <- range(Date_FUN(x$tip.label),na.rm=T) #the dates of the root and tip, respectively

#       m=(real_time[2]-real_time[1])/(tree_time[2]-tree_time[1]) #the gradient between tree time and the time specified on the tips
      m=(real_time[2]-real_time[1])/(tip_time[2]-tip_time[1]) #the gradient between tree time and the time specified on the tips

      sk$time<-real_time[2]+m*sk$time
#       sk$time<-real_time[1]
#         m*(sk$time-tip_time[1])
    } else {
      if(timeForwards)
      {
        sk$time = scaling*(max(sk$time) - sk$time)
      } else {
        sk$time = scaling*(sk$time)
      }
    }




    #list(Skyline=sk, MolClockRate=clock.rates[i])
    list(Skyline=sk)
  })

  if(output_type=="list")
  {
    return(skyline_list)
  }

  time_mat<-sapply(1:length(skyline_list),function(x) skyline_list[[x]]$Skyline$time)
  pop_mat<-sapply(1:length(skyline_list),function(x) skyline_list[[x]]$Skyline$population.size)

  if (output_type=="matrices") {
    return(list(time_mat=time_mat,pop_mat=pop_mat))
  }


  ###Make master time and population objects

  treenum_mat<-matrix(1:length(skyline_list),nrow=nrow(time_mat),ncol=length(skyline_list),byrow=T)
  posnum_mat<-matrix(1:nrow(time_mat),nrow=nrow(time_mat),ncol=length(skyline_list),byrow=F)

  treenum_vec<-treenum_mat[order(time_mat)]
  posnum_vec<-posnum_mat[order(time_mat)]

  master_time_vec=sort(time_mat)
  master_time_vec=c(master_time_vec[1]-1e-300,master_time_vec)
  master_pop_mat<-matrix(NA,nrow=length(master_time_vec),ncol=length(skyline_list))
  master_pop_mat[1,]<-0

  s=1
  for(s in 1:length(time_mat))
  {
    master_pop_mat[s+1,]<-master_pop_mat[s,]
    master_pop_mat[s+1,treenum_vec[s]]<-pop_mat[posnum_vec[s],treenum_vec[s]]
  }

  # rownames(master_pop_mat)<-master_time_vec
  if(output_type=="master")
  {
    out_temp<-cbind(master_time_vec,master_pop_mat)
    colnames(out_temp)<-c("Time",paste("Tree",1:ncol(master_pop_mat),sep=""))
    return(out_temp)
  }

  conf.int_mat<-cbind(master_time_vec,t(apply(master_pop_mat,1, function(x) quantile(x,probs=c(0.5,0.025,0.975)))))
  rownames(conf.int_mat)<-NULL
  colnames(conf.int_mat)[1]<-x_lab

  if(output_type=="conf.int")
  {
    return(conf.int_mat)
  }

  if(output_type=="conf.int.plot")
  {
    conf.int.skyline(conf.int_mat,plot_type=plot_type,epsilon=epsilon,...)
  }

  if(missing(output_type)) stop("no output type specified")
}


#' @title conf.int.skyline
#' @description Makes a skyline plot from a matrix with medians and confidence intervals as produced by Phy2Sky
#'
#' @author George Shirreff <georgeshirreff@@gmail.com>
#'
#' @return A plot with median in black and confidence intervals in red
#'
#' @param conf.int: this type of output from Phy2Sky
#' @param epsilon: a smoothing parameter. Any time intervals smaller than this value will be joined with adjacent bins until all are wider than this.
#' @param plot_type: choose between "step", where adjacent points on the skyline are joined by a vertical and a horizontal line, and "linear", where adjacent points on the skyline are joined by a single sloping straight line.
#'
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples conf.int<-Phy2Sky(trees,output_type="conf.int")


conf.int.skyline<-function(conf.int,epsilon=0,plot_type=c("step"),...)
{
  end_times<-conf.int[,1]

  plot(type="n",x=range(end_times),y=range(conf.int[,-1]),xlab=colnames(conf.int)[1],ylab=expression(N[e]*tau),cex.lab=1,cex.axis=1,...)

  skyline_smooth(end_times,pop_size=conf.int[,2],col=1,epsilon=epsilon,plot_type=plot_type)
  skyline_smooth(end_times,pop_size=conf.int[,3],col=2,epsilon=epsilon,plot_type=plot_type)
  skyline_smooth(end_times,pop_size=conf.int[,4],col=2,epsilon=epsilon,plot_type=plot_type)
}



#' @title skyline_smooth
#' @description Draws a generalized skyline onto a plot
#'
#' @author George Shirreff <georgeshirreff@@gmail.com>
#'
#' @return Draws the lines of a skyline, but with a minimum time interval specified by epsilon
#' @return If draw=F, it returns a list with the end times of the intervals and the population sizes, smoothed as per the value of epsilon.
#'
#' @inheritParams draw_skyline
#' @param epsilon: a smoothing parameter. Any time intervals smaller than this value will be joined with adjacent bins until all are wider than this.
#' @param draw: logical specifying if the skyline should be drawn (TRUE) or simply the points on the graph returned (FALSE)
#'
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples conf.int<-Phy2Sky(trees,output_type="conf.int")
#' @examples end_times<-as.numeric(rownames(conf.int))
#' @examples plot(type="n",x=range(end_times),y=range(conf.int),ylab=expression(N[e]*tau))
#' @examples skyline_smooth(end_times,pop_size=conf.int[,1],col=1,epsilon=0.1)

skyline_smooth<-function(end_times,pop_size,epsilon=0,plot_type=c("step"),col=4)
{
  if(diff(range(end_times))<epsilon) stop("epsilon is bigger than the entire skyline range")
  if(epsilon==0)
  {
    end_times_out<-end_times
    pop_size_out<-pop_size
  } else {
    harmonic<-function(vec)  1/mean(1/vec)

    len<-length(pop_size)

    if(length(end_times)!=len) stop("vectors are different lengths")

    end_times_out=numeric(0)
    pop_size_out=numeric(0)
    intervals=c(0,diff(end_times)) #interval sizes

    #cbind(1:length(end_times), end_times,intervals,pop_size)

    i=len
    while(i>=1)
    {
      combine_vec=i #the start of the intervals to be combined, even if only a single one
      while(sum(intervals[combine_vec])<epsilon) #if the combined interval is less than epsilon
      #shouldn't need to include an i-limit here, because the combine_interval should always be greater than one (provided epsilon isn't bigger than the total time!!)
      {
        if(i>1)
        {
          i=i-1
          combine_vec=c(combine_vec,i) #add another interval to the combined interval
        } else {
          combine_vec=c(last_combine_vec,combine_vec) #if there aren't any left, add the previous run of intervals
          end_times_out<-end_times_out[-1] #if this is the case, you have to take out the last entry of each vector so a combined value can be put in
          pop_size_out<-pop_size_out[-1]
        }
      }

      end_times_out<-c(end_times[max(combine_vec)],end_times_out)
#       these_pop_sizes=pop_size[combine_vec]
#       these_pop_sizes=these_pop_sizes[these_pop_sizes!=0]
      #pop_size_out<-c(harmonic(these_pop_sizes),pop_size_out)
      pop_size_out<-c(harmonic(pop_size[combine_vec]),pop_size_out)

      last_combine_vec=combine_vec #in case the last run of intervals don't add up to >epsilon and you need to combine with the previous run of intervals

      i=i-1
    }
  }
  draw_skyline(end_times_out,pop_size_out,plot_type=plot_type,col=col)
}



#' @title draw_skyline
#' @description Draws a classic skyline onto a plot
#'
#' @author George Shirreff <georgeshirreff@@gmail.com>
#'
#' @return Draws the lines of a skyline (step) or INLA (piece-wise linear)
#'
#' @param end_times: the times at the end of every interval in the skyline.
#' @param pop_size: the effective population size at that time point.
#' @param plot_type: choose between "step", where adjacent points on the skyline are joined by a vertical and a horizontal line, and "linear", where adjacent points on the skyline are joined by a single sloping straight line.
#'
#' @examples library(ape)
#' @examples trees <- rmtree(N=5,n=20)
#' @examples conf.int<-Phy2Sky(trees,output_type="conf.int")
#' @examples end_times<-as.numeric(rownames(conf.int))
#' @examples plot(type="n",x=range(end_times),y=range(conf.int),ylab=expression(N[e]*tau))
#' @examples draw_skyline(end_times,pop_size=conf.int[,1],col=1)

draw_skyline<-function(end_times,pop_size,plot_type=c("linear","step"),col=1)
{
  len<-length(pop_size)
  if(length(end_times)!=len) stop("vectors are different lengths")

  if(plot_type=="step")
  {
    arrows(x0=end_times[-len],x1=end_times[-1],y0=pop_size[-len],code=0,col=col)
    arrows(x0=end_times[-1],y0=pop_size[-len],y1=pop_size[-1],code=0,col=col)
  } else if (plot_type=="linear") {
    arrows(x0=end_times[-len],x1=end_times[-1],y0=pop_size[-len],y1=pop_size[-1],code=0,col=col)
  }
}

