\name{mixcure.penal.mi}
\alias{mixcure.penal.mi}
\title{
Parameter estimation of mixture cure model for multiply imputed data
}
\description{
Parameter estimation of MC model for mids object, via a nonlinear minimization function ('nlm'), which utilizes a newton-type algorithm for optimization of all parameters. Optimization is based on either the usual likelihood or the Firth-type penalized likelihood (FT-PL).
}
\usage{
mixcure.penal.mi(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A mids object from mice().
}
}
\details{
Refers to mixcure.penal.est() in 'mixcuref'.
}
\value{
Refers to mixcure.penal.est() in 'mixcuref'.
}
\references{
Xu C, Bull SB. Penalized maximum likelihood inference under the mixture cure model in sparse data. Statistics in Medicine 2023; 42: 2134–2161.
}
\author{
Changchang Xu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(ANNbcBMdat5_miss)
library(mice)

  imp20 <- mice(ANNbcBMdat5_miss, m=20, maxit = 20, seed = 24415, print=FALSE) #20 times imputation
  data20 <- complete(imp20, action = "long", include = TRUE)
  #note labER 0-positive
  data20$TN<-factor(with(data20,
                         ifelse((Her2=="1" & ER == 2 & PR == 2),"2","1")))
  data20$LuminalA<-factor(with(data20,
                               ifelse((Her2=="1" & (ER == 1 | PR == 1) & KI67 == 1),"2","1")))
  data20$LuminalB<-factor(with(data20,
                               ifelse((Her2=="1" & (ER == 1 | PR == 1) & KI67 == 2),"2","1")))
  imp20sub <- as.mids(data20)

 # CLIP inference for maximum likelihood estimation

 mixcure.all20<- with(data=imp20sub, mixcure.penal.mi(Surv(Time,CENS==1)~ Her2 + LuminalA + TN + MENS0 + TUMCT, init=c(2, rep(0.1,5),-10,rep(0.1,5),1), pl=F))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
