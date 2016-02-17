\name{simG}
\alias{simG}
\docType{data}
\title{simG}
\description{Genind object. Microsatellite data set for 50 diploid individuals: two loci, sampled in the same generation. Simulated with rmetasim: sampled after 1000 generations. Mu = 0.01}
\usage{data("simG")}
\format{
  The format is:
Formal class 'genind' [package "adegenet"] with 14 slots
  ..@ strata   : NULL
  ..@ hierarchy: NULL
  ..@ tab      : int [1:50, 1:12] 1 0 1 1 1 0 0 1 1 1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:50] "331" "247" "865" "256" ...
  .. .. ..$ : chr [1:12] "L1.0" "L1.1" "L1.2" "L1.3" ...
  ..@ loc.names: chr [1:2] "L1" "L2"
  ..@ loc.fac  : Factor w/ 2 levels "L1","L2": 1 1 1 1 1 1 2 2 2 2 ...
  ..@ loc.nall : Named int [1:2] 6 6
  .. ..- attr(*, "names")= chr [1:2] "L1" "L2"
  ..@ all.names:List of 2
  .. ..$ L1: chr [1:6] "0" "1" "2" "3" ...
  .. ..$ L2: chr [1:6] "0" "1" "2" "3" ...
  ..@ call     : language .local(x = x, i = i, j = j, drop = drop)
  ..@ ind.names: chr [1:50] "331" "247" "865" "256" ...
  ..@ pop      : Factor w/ 1 level "0": 1 1 1 1 1 1 1 1 1 1 ...
  ..@ pop.names: chr "0"
  ..@ ploidy   : int [1:50] 2 2 2 2 2 2 2 2 2 2 ...
  ..@ type     : chr "codom"
  ..@ other    : list()
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
data(simG)
## maybe str(simG) ; plot(simG) ...
}
\keyword{datasets}
