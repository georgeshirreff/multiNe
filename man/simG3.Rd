\name{simG3}
\alias{simG3}
\docType{data}
\title{
%%   simG3
}
\description{
%%  Genind object. Microsatellite data set: two loci, 50 individuals sampled each in three consecutive generations. Simulated with rmetasim: sampled after 1000 generations. Mu = 0.01
}
\usage{data("simG3")}
\format{
  The format is:
Formal class 'genind' [package "adegenet"] with 14 slots
  ..@ strata   : NULL
  ..@ hierarchy: NULL
  ..@ tab      : int [1:150, 1:10] 1 0 1 1 1 0 0 1 1 1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:150] "331" "247" "865" "256" ...
  .. .. ..$ : chr [1:10] "L1.0" "L1.1" "L1.4" "L1.2" ...
  ..@ loc.names: chr [1:2] "L1" "L2"
  ..@ loc.fac  : Factor w/ 2 levels "L1","L2": 1 1 1 1 1 2 2 2 2 2
  ..@ loc.nall : Named int [1:2] 5 5
  .. ..- attr(*, "names")= chr [1:2] "L1" "L2"
  ..@ all.names:List of 2
  .. ..$ L1: chr [1:5] "0" "1" "4" "2" ...
  .. ..$ L2: chr [1:5] "2" "3" "1" "0" ...
  ..@ call     : language repool(g0S, g1S, g2S)
  ..@ ind.names: chr [1:150] "331" "247" "865" "256" ...
  ..@ pop      : Factor w/ 3 levels "0","1","2": 1 1 1 1 1 1 1 1 1 1 ...
  ..@ pop.names: chr [1:3] "0" "1" "2"
  ..@ ploidy   : int [1:150] 2 2 2 2 2 2 2 2 2 2 ...
  ..@ type     : chr "codom"
  ..@ other    : NULL
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(simG3)
## maybe str(simG3) ; plot(simG3) ...
}
\keyword{datasets}
