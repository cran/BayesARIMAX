
BayesARIMAX<-function(Y,X,sd=10,iter=100,burnIn=40,nc=2,p=1,d=1,q=1)
{
    requireNamespace("coda")
  requireNamespace("forecast")
  model_arimax<-arima(Y,order = c(p,d,q),xreg =X, include.mean = TRUE,  method = c( "ML"),kappa = 1e+06)
  par<-model_arimax$coef
  names(par) <- NULL
  para<-0
  paraprior<-0
  prior_total<-0
  l<-length(par)
  llike<-model_arimax$loglik
  prior <- function(param)
  {
    List <- list()
    for(j in 1:l){
      para[j]= param[j]
      paraprior[j] = dnorm(para[j], sd, log = T)
      List[[j]] <- paraprior[j]
      prior_total<-prior_total+List[[j]]
    }

    return(prior_total)
  }
  posterior <- function(param)
  {
    return (llike + prior(param))
  }
  proposalfunction <- function(param){
    return(rnorm(l,mean = param, sd))
  }
  run_metropolis_MCMC <- function(startvalue, iterations){
    chain = array(dim = c(iterations+1,l))
    chain[1,] = startvalue
    for (i in 1:iterations){
      proposal = proposalfunction(chain[i,])
      probab = exp(posterior(proposal) - posterior(chain[i,]))
      if (runif(1) < probab){
        chain[i+1,] = proposal
      }else{chain[i+1,] = chain[i,]
      }
    }
    return(mcmc(chain))
  }

  List <- list()
  for(i in 1:nc){
    chain= run_metropolis_MCMC(par, iter)
    acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
    List[[i]] <- chain
  }
  summary(List[[i]])

  Diagnostic<-gelman.diag(List, confidence = 0.99, transform=FALSE, autoburnin=TRUE,
                          multivariate=TRUE)
  gelman.plot(List)
  return_list<-list(model_arimax,summary(chain),summary(List[[1]]),Diagnostic,gelman.plot(List))
  return(return_list)
}
