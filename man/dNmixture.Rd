% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dNmixture.R
\name{dNmixture}
\alias{dNmixture}
\alias{dNmixture_s}
\alias{dNmixture_v}
\alias{rNmixture_s}
\alias{rNmixture_v}
\alias{dNmixture_BNB_v}
\alias{dNmixture_BNB_s}
\alias{dNmixture_BNB_oneObs}
\alias{dNmixture_BBP_v}
\alias{dNmixture_BBP_s}
\alias{dNmixture_BBP_oneObs}
\alias{dNmixture_BBNB_v}
\alias{dNmixture_BBNB_s}
\alias{dNmixture_BBNB_oneObs}
\alias{rNmixture_BNB_v}
\alias{rNmixture_BNB_s}
\alias{rNmixture_BNB_oneObs}
\alias{rNmixture_BBP_v}
\alias{rNmixture_BBP_s}
\alias{rNmixture_BBP_oneObs}
\alias{rNmixture_BBNB_v}
\alias{rNmixture_BBNB_s}
\alias{rNmixture_BBNB_oneObs}
\title{dNmixture distribution for use in \code{nimble} models}
\usage{
dNmixture_v(x, lambda, prob, Nmin = -1, Nmax = -1, len, log = 0)

dNmixture_s(x, lambda, prob, Nmin = -1, Nmax = -1, len, log = 0)

rNmixture_v(n, lambda, prob, Nmin = -1, Nmax = -1, len)

rNmixture_s(n, lambda, prob, Nmin = -1, Nmax = -1, len)

dNmixture_BNB_v(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len, log = 0)

dNmixture_BNB_s(x, lambda, theta, prob, Nmin = -1, Nmax = -1, len, log = 0)

dNmixture_BNB_oneObs(x, lambda, theta, prob, Nmin = -1, Nmax = -1, log = 0)

dNmixture_BBP_v(x, lambda, prob, s, Nmin = -1, Nmax = -1, len, log = 0)

dNmixture_BBP_s(x, lambda, prob, s, Nmin = -1, Nmax = -1, len, log = 0)

dNmixture_BBP_oneObs(x, lambda, prob, s, Nmin = -1, Nmax = -1, log = 0)

dNmixture_BBNB_v(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len, log = 0)

dNmixture_BBNB_s(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len, log = 0)

dNmixture_BBNB_oneObs(x, lambda, theta, prob, s, Nmin = -1, Nmax = -1, log = 0)

rNmixture_BNB_v(n, lambda, theta, prob, Nmin = -1, Nmax = -1, len)

rNmixture_BNB_s(n, lambda, theta, prob, Nmin = -1, Nmax = -1, len)

rNmixture_BNB_oneObs(n, lambda, theta, prob, Nmin = -1, Nmax = -1)

rNmixture_BBP_v(n, lambda, prob, s, Nmin = -1, Nmax = -1, len)

rNmixture_BBP_s(n, lambda, prob, s, Nmin = -1, Nmax = -1, len)

rNmixture_BBP_oneObs(n, lambda, prob, s, Nmin = -1, Nmax = -1)

rNmixture_BBNB_v(n, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)

rNmixture_BBNB_s(n, lambda, theta, prob, s, Nmin = -1, Nmax = -1, len)

rNmixture_BBNB_oneObs(n, lambda, theta, prob, s, Nmin = -1, Nmax = -1)
}
\arguments{
\item{x}{vector of integer counts from a series of sampling occasions.}

\item{lambda}{expected value of the Poisson distribution of true abundance}

\item{prob}{detection probability (scalar for \code{dNmixture_s}, vector for
\code{dNmixture_v}).}

\item{Nmin}{minimum abundance to sum over for the mixture probability. Set to
-1 to select automatically (not available for beta binomial variations; see
Details).}

\item{Nmax}{maximum abundance to sum over for the mixture probability. Set to
-1 to select automatically (not available for beta binomial variations; see
Details).}

\item{len}{The length of the x vector}

\item{log}{TRUE or 1 to return log probability. FALSE or 0 to return
probability.}

\item{n}{number of random draws, each returning a vector of length
\code{len}. Currently only \code{n = 1} is supported, but the argument
exists for standardization of "\code{r}" functions.}

\item{theta}{abundance overdispersion parameter required for negative
binomial (*NB) N-mixture models. The negative binomial is parameterized
such that variance of x is \code{lambda^2 * theta + lambda}}

\item{s}{detection overdispersion parameter required for beta binomial (BB*)
N-mixture models. The beta binomial is parameterized such that variance of
x is \code{V(x) = N * prob * (1-prob) * (N + s) / (s + 1)}}
}
\value{
For \code{dNmixture_s} and \code{dNmixture_v}: the probability (or likelihood) or log
probability of observation vector \code{x}.

For \code{rNmixture_s} and \code{rNmixture_v}: a simulated detection history, \code{x}.
}
\description{
\code{dNmixture_s} and \code{dNmixture_v} provide Poisson-Binomial mixture
distributions of abundance ("N-mixture") for use in \code{nimble} models.
Overdispersion alternatives using the negative binomial distribution (for
the abundance submodel) and the beta binomial distribution (for the detection
submodel) are also provided.
}
\details{
These nimbleFunctions provide distributions that can be
    used directly in R or in \code{nimble} hierarchical models (via
    \code{\link[nimble]{nimbleCode}} and
    \code{\link[nimble]{nimbleModel}}).

An N-mixture model defines a distribution for multiple counts (typically of
animals, typically made at a sequence of visits to the same site).  The
latent number of animals available to be counted, N, follows a Poisson or
negative binomial distribution. Each count, \code{x[i]} for visit \code{i},
follows a binomial or beta-binomial distribution. The N-mixture distributions
calculate the marginal probability of observed counts by summing over the
range of latent abundance values.

The basic N-mixture model uses Poisson latent abundance with mean
\code{lambda} and binomial observed counts  with size (number of trials) N
and probability of success (being counted) \code{prob[i]}. This distribution
is available in two forms, \code{dNmixture_s} and \code{dNmixture_v}. With
\code{dNmixture_s}, detection probability is a scalar, independent of visit,
so \code{prob[i]} should be replaced with \code{prob} above.  With
\code{dNmixture_v}, detection probability is a vector, with one element for
each visit, as written above.

We also provide three important variations on the traditional N-mixture
model: \code{dNmixture_BNB}, \code{dNmixture_BBP}, and \code{dNmixture_BBNB}.
These distributions allow you to replace the Poisson (P) abundance
distribution with the negative binomial (NB) and the binomial (B) detection
distribution with the beta binomial (BB).

Binomial-negative binomial: BNB N-mixture models use a binomial distribution
for detection and a negative binomial distribution for abundance with scalar
overdispersion parameter \code{theta} (0-Inf). We parameterize such that the
variance of the negative binomial is \code{lambda^2 * theta + lambda}, so
large \code{theta} indicates a large amount of overdisperison in abundance.
The BNB is available in three suffixed forms: \code{dNmixture_BNB_v} is used
if \code{prob} varies between observations, \code{dNmixture_BNB_s} is used if
\code{prob} is scalar (constant across observations), and
\code{dNmixture_BNB_oneObs} is used if only one observation is available at
the site (so both x and prob are scalar).

Beta-binomial-Poisson: BBP N-mixture uses a beta binomial distribution for
detection probabilities and a Poisson distribution for abundance. The beta
binomial distribution has scalar overdispersion parameter s (0-Inf). We
parameterize such that the variance of the beta binomial is \code{N * prob
* (1-prob) * (N + s) / (s + 1)}, with greater s indicating less variance
(greater-than-binomial relatedness between observations at the site) and s ->
0 indicating the binomial. The BBP is available in three suffixed forms:
\code{dNmixture_BBP_v} is used if \code{prob} varies between observations,
\code{dNmixture_BBP_s} is used if \code{prob} is scalar (constant across
observations), and \code{dNmixture_BBP_oneObs} is used if only one
observation is available at the site (so both x and prob are scalar).

Beta-binomial-negative-binomial: dNmixture_BBNB is available using a negative
binomial abundance distribution and a beta binomial detection distribution.
\code{dNmixture_BBNB} is available with \code{_s}, \code{_v}, and
\code{_oneObs} suffixes as above and requires both arguments \code{s} and
\code{theta} as parameterized above.

The distribution dNmixture_oneObs is not provided as the probability given
by the traditional N-mixture distribution for \code{length(x) = 1} is
equivalent to \code{dpois(prob * lambda)}.

For more explanation, see package vignette
(\code{vignette("Introduction_to_nimbleEcology")}).

Compared to writing \code{nimble} models with a discrete latent state of
abundance N and a separate scalar datum for each count, use of these
distributions allows one to directly sum (marginalize) over the discrete
latent state N and calculate the probability of all observations for a site
jointly.

If one knows a reasonable range for summation over possible values of N, the
start and end of the range can be provided as \code{Nmin} and \code{Nmax}.
Otherwise one can set both to -1, in which case values for \code{Nmin} and
\code{Nmax} will be chosen based on the 0.0001 and 0.9999 quantiles of the
marginal distributions of each count, using the minimum over counts of the
former and the maximum over counts of the latter.

The summation over N uses the efficient method given by Meehan et al. (2020,
see Appendix B) for the basic Poisson-Binomial case, extended for the
overdispersion cases in Goldstein and de Valpine (2022).

These are \code{nimbleFunction}s written in the format of user-defined
distributions for NIMBLE's extension of the BUGS model language. More
information can be found in the NIMBLE User Manual at
\href{https://r-nimble.org}{https://r-nimble.org}.

When using these distributions in a \code{nimble} model, the left-hand side
will be used as \code{x}, and the user should not provide the \code{log}
argument.

For example, in \code{nimble} model code,

\code{observedCounts[i, 1:T] ~ dNmixture_v(lambda[i],
prob[i, 1:T],
Nmin, Nmax, T)}

declares that the \code{observedCounts[i, 1:T]} (observed counts
for site \code{i}, for example) vector follows an N-mixture
distribution with parameters as indicated, assuming all the
parameters have been declared elsewhere in the model. As above,
\code{lambda[i]} is the mean of the abundance distribution at site
i, \code{prob[i, 1:T]} is a vector of detection probabilities at
site i, and \code{T} is the number of observation occasions. This
will invoke (something like) the following call to
\code{dNmixture_v} when \code{nimble} uses the model such as for
MCMC:

\code{dNmixture_v(observedCounts[i, 1:T], lambda[i],
prob[i, 1:T],
Nmin, Nmax, T, log = TRUE)}

If an algorithm using a \code{nimble} model with this declaration
needs to generate a random draw for \code{observedCounts[1:T]}, it
will make a similar invocation of \code{rNmixture_v}, with \code{n = 1}.

If the observation probabilities are visit-independent, one would use:

\code{observedCounts[1:T] ~ dNmixture_s(observedCounts[i, 1:T], lambda[i],
prob[i],
Nmin, Nmax, T)}
}
\section{Notes for use with automatic differentiation}{


The N-mixture distributions are the only ones in \code{nimbleEcology} for which
one must use different versions when AD support is needed. See
\code{\link{dNmixtureAD}}.
}

\examples{
# Set up constants and initial values for defining the model
len <- 5 # length of dataset
dat <- c(1,2,0,1,5) # A vector of observations
lambda <- 10 # mean abundance
prob <- c(0.2, 0.3, 0.2, 0.1, 0.4) # A vector of detection probabilities

# Define code for a nimbleModel
 nc <- nimbleCode({
   x[1:5] ~ dNmixture_v(lambda, prob = prob[1:5],
                        Nmin = -1, Nmax = -1, len = 5)

   lambda ~ dunif(0, 1000)

   for (i in 1:5) {
     prob[i] ~ dunif(0, 1)
   }
 })

# Build the model
nmix <- nimbleModel(nc,
                    data = list(x = dat),
                    inits = list(lambda = lambda,
                                 prob = prob))
# Calculate log probability of data from the model
nmix$calculate()
# Use the model for a variety of other purposes...
}
\references{
D. Turek, P. de Valpine and C. J. Paciorek. 2016. Efficient Markov chain Monte
Carlo sampling for hierarchical hidden Markov models. Environmental and
Ecological Statistics 23:549–564. DOI 10.1007/s10651-016-0353-z

Meehan, T. D., Michel, N. L., & Rue, H. 2020. Estimating Animal Abundance
with N-Mixture Models Using the R—INLA Package for R. Journal of Statistical
Software, 95(2). https://doi.org/10.18637/jss.v095.i02

Goldstein, B.R., and P. de Valpine. 2022. Comparing N-mixture Models and
GLMMs for Relative Abundance Estimation in a Citizen Science Dataset.
Scientific Reports 12: 12276. DOI:10.1038/s41598-022-16368-z
}
\seealso{
For occupancy models dealing with detection/nondetection data,
see \code{\link{dOcc}}.
For dynamic occupancy, see \code{\link{dDynOcc}}.
}
\author{
Ben Goldstein, Lauren Ponisio, and Perry de Valpine
}
