\name{pool.mixcure}
\alias{pool.mixcure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Combine estimates by pooling via Rubin's Rule
}
\description{
The pool.mixcure() function combines the estimates from m repeated complete data analyses. The typical sequence of steps to perform a multiple imputation analysis is:
}
\usage{
pool.mixcure(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class mira (produced by with.mids() or as.mira()), or a list with model fits.
}
}
\details{
The pool.mixcure() function averages the estimates of the complete data model, computes the total variance over the repeated analyses by Rubin's rules (Rubin, 1987, p. 76), and computes the following diagnostic statistics per estimate:

1. Relative increase in variance due to nonresponse r;

2. Residual degrees of freedom for hypothesis testing df;

3. Proportion of total variance due to missingness lambda;

4. Fraction of missing information fmi.

The details please refer to pool() under 'mice' package.
}
\value{
An object of class mipo, which stands for 'multiple imputation pooled outcome'.
}
\references{
Rubin, D.B. (1987). Multiple Imputation for Nonresponse in Surveys. New York: John Wiley and Sons.

(submitted) Xu C, Briollais L, Andrulis I, Bull SB.Inference for multiple imputation with penalized maximum likelihood under the mixture cure model. Statistical Methods in Medical Research.
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

out.RR.mod <- pool.mixcure(mixcure.all20)

 mixcure.all20.penal<- with(data=imp20sub, mixcure.penal.mi(Surv(Time,CENS==1)~ Her2 + LuminalA + TN + MENS0 + TUMCT, init=c(2, rep(0.1,5),-10,rep(0.1,5),1), pl=T))

out.RR.mod.penal <- pool.mixcure(obj=mixcure.all20.penal)

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
