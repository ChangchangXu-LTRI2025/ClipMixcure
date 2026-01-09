\name{mixcure.clip.pdf}
\alias{mixcure.clip.pdf}
\title{
Probability density calculation for mixture cure model for CLIP method
}
\description{
A function that calculates pdf for the fixed input value of the parameter of interest under mixture cure model. This function is automatically called in clip.mixcure() function.
}
\usage{
mixcure.clip.pdf(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
MC model specification in the form of a conventional survival model, i.e., Surv()~. Note the whatever x variables specified in the latency part of the model are assumed to be also regressed on in the incidence part, and vice versa.
}
  \item{data}{
  a data object in the form of dataframe.
  }
  \item{init}{
  a vector of initial values input for the optimization of the vector of parameters specified in the likelihood. Usually, specified in the order of 1) regression parameters for the incidence part, 2) regression parameters for the latency part, and 3) the shape parameter of Weibull hazard function for the latency part.
  }
  \item{pl}{
  if TRUE, uses FT-PL; otherwise uses the usual likelihood for the model.
  }
  \item{iterlim}{
  specifies the number of maximum iterations for model parameter optimization. DEFAULT set to 200.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
This function returns a probability value between 0 ~ 1.
}
\references{
%% ~put references to the literature/web site here ~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x)
{
  }
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
