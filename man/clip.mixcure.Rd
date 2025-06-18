\name{clip.mixcure}
\alias{clip.mixcure}
\title{
combined likelihood profile inference for mixture cure models under multiple imputation
}
\description{
This function implements the combination of likelihood profiles (CLIP) method described in Xu et al (2024) for mixture cure model. This method is useful for computing confidence intervals for parameters after multiple imputation of data sets, if the normality assumption on parameter estimates and consequently the validity of applying Rubin’s rules (pooling of variances) is in doubt. The assumption on approximating endpoints of confidence intervals and performing hypothesis tests follow the same idea as paper by Heinze (2013). The default confidence level of the combined/posterior profile distribution is equal to the values specified in the argument ci.level (usually 0.025 and 0.975). The function also accomodates two alternative endpoint estimations by Firth type penalized likelihood (Xu & Bull (2023)) and maximum likelihood.
}
\usage{
clip.mixcure(obj = NULL, variable = NULL, pl = F, ci.level = c(0.025, 0.975), pvalue = TRUE, bound.lo = NULL, bound.up = NULL, iternum = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
The result of analysis of a mice (multiply imputed) object via with(data=, mixcure.penal.mi())
}
  \item{variable}{
The variables for which confidence intervals need be computed. If missing, confidence intervals for all variables will be computed.
}
  \item{pl}{
Whether Firth-type penalized likelihood is applied for point estimation.
}
  \item{ci.level}{
The two confidence levels for the endpoints of the confidence interval.
}
  \item{pvalue}{
If TRUE, will also compute a P-value from the posterior.
}
  \item{bound.lo}{
Bounds (vector of length 2) for the lower limit. Can be left blank. Use only if
problems are encountered.
}
  \item{bound.up}{
Bounds (vector of length 2) for the upper limit. Can be left blank. Use only if
problems are encountered.
}
  \item{iternum}{
Limit for iteration number for calculating confidence interval.
}
}
\details{
Require 'mice' package.
For each confidence limit, this function performs a binary search to evaluate the combined posterior, which is obtained by first transforming the imputed-data likelihood profiles into cumulative distribution functions (CDFs), and then averaging the CDFs to obtain the CDF of the posterior, see Xu(2024) and Heinze(2013).
}
\value{

}
\references{
Heinze G, Ploner M, Beyea J. Confidence intervals after multiple imputation: Combining profile likelihood information from logistic regressions. Statistics in Medicine 2013; 32: 5062–5076.

Xu C, Bull SB. Penalized maximum likelihood inference under the mixture cure model in sparse data. Statistics in Medicine 2023; 42: 2134–2161.

(submitted) Xu C, Briollais L, Andrulis I, Bull SB.Inference for multiple imputation with penalized maximum likelihood under the mixture cure model. Statistical Methods in Medical Research.

}
\author{
Changchang Xu
}
\note{
}
%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Low event rate, 5 variable
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

 out.clip.mod <- clip.mixcure(obj=mixcure.all20, pl=F)

 # CLIP inference for Firth type penalized likelihood estimation

 mixcure.all20.penal<- with(data=imp20sub, mixcure.penal.mi(Surv(Time,CENS==1)~ Her2 + LuminalA + TN + MENS0 + TUMCT, init=c(2, rep(0.1,5),-10,rep(0.1,5),1), pl=T))

 out.clip.mod.penal <- clip.mixcure(obj=mixcure.all20, pl=T)
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
