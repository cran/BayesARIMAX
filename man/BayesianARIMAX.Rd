\name{BayesARIMAX}
\alias{BayesARIMAX}
\title{Bayesian Estimation of ARIMAX Model}
\usage{BayesARIMAX(Y,X,sd,iter,burnIn,nc,p,d,q)}
\description{The Autoregressive Integrated Moving Average (ARIMA) model is very popular univarite time series model. Its application has been widened by the incorporation of exogenous variable(s) (X) in the model and modified as ARIMAX . The details of the model are documented  effectively by Bierens (1987) <doi:10.1016/0304-4076(87)90086-8>. In this package we estimate the ARIMAX model using Bayesian framework. We have assumed priors to follow Normal distribution and the posteriors are obtained using Markov chain Monte Carlo (MCMC) algorithm.
}

%-
\arguments{
  \item{Y}{
a univariate time series data
}
\item{X}{
a univariate or multivariate time series data to be used as covariate or exogenous variable for Y

}
\item{sd}{
sd is the standard deviation of the Normal priors assigned to each parameter. Default value is taken as 10.

}
\item{iter}{
iter is the number of iterations for the Markov chain Monte Carlo (MCMC) chain. Default value is 100.

}
\item{burnIn}{
burnIn is the number of iterations to be excluded from the estimate of the chain. Default value is 40.

}
\item{nc}{
nc is the number of chains used for calculating the diagnostic statistics. Default value is 2.

}
\item{p}{
p is the order of AR parameter of ARIMA model. Default value is 1.

}
\item{d}{
d is the order of differencing used for making the series sationary. Default value is 1.

}
\item{q}{
q is the he order of MA parameter of ARIMA model. Default value is 1.

}
}
\details{The Autoregressive Integrated Moving Average (ARIMA) model is very popular univariate time sereis model. The details of the model along with its implementation using R has been well documented by Shumway and Stoffer (2017) <https://doi.org/10.1007/978-3-319-52452-8> Its application has been widened by the incorpoartion of exogenous variable(s) (X) in the model and modified as ARIMAX.In this package we have estimated this model using Bayesian technique. Metropilis-Hasting algorithm is used to generate the posterior density of the model parameters. Normal distribution is used as priors for each parameter of the model following Fioruci et al., (2014).

}
\value{It returns the Bayesian estimates of the ARIMAX model.}
\references{
%%
Bierens, H.J.(1987)<https://doi.org/10.1016/0304-4076(87)90086-8>

Fioruci et al.(2014)<https://doi.org/10.1080/02664763.2013.839635>

Metropolis et al.(1953)<https://doi.org/10.1063/1.1699114>

Shumway, R.H. and Stoffer, D.S.(2017)<https://doi.org/10.1007/978-3-319-52452-8>

}
\author{
Achal Lama,Kn Singh and Bishal Gurung
}
\note{This package cab be used to analyse ARIMAX model using Normal priors for the parameters. The users need to identify the tentative ARIMA model by themselves.}

\seealso{arimax,arima}
\examples{
set.seed(121)
Y<- arima.sim(list(order = c(1,1,1),ar=0.7,ma=0.4), n = 49)
X=rnorm(50,4,1)
BayesARIMAX(Y,X)
}

\keyword{Bayesian}
\keyword{ARIMAX}

