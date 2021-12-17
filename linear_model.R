# (Naive) Linear Model
#
# y_i ~ Normal(mu, 1/tau)
# 
# mu = alpha + beta * x
#
# Priors
# alpha ~ Normal(alpha0, 1/tau_alpha)
# beta ~ Normal(beta0, 1/tau_beta)
# tau ~ Gamma(a, b)

# Set "true" values of 
N = 20


# hyperparameters
alpha0 = 0
beta0 = 0
tau_alpha = 0.01
tau_beta = 0.01
a = 0.01
b = 0.01

# Simulate some data
x = runif(N, 0, 15)
x_center = x - mean(x)
mu = 0.5 + 0.1 * x_center
y = rnorm(N, mu, sd=1.2)
plot(x_center, y)

# Run the MCMC for `its` iterations
its=1000


# Create conditional posterior density functions
joint_log_posterior = function(alpha, beta, tau) {
  
  mu = alpha + beta * x_center
  # Calculate the likelihood
  log_likelihood = sum(dnorm(y, mean=mu, sd=sqrt(1/tau), log=T))
  
  # Priors
  log_prior_alpha = dnorm(alpha, mean=alpha0, sd=sqrt(1/tau_alpha), log=T)
  log_prior_beta = dnorm(beta, mean=beta0, sd=sqrt(1/tau_beta), log=T)
  log_prior_tau = dgamma(tau, shape=a, rate=b, log=T)
  
  log_posterior = log_likelihood + 
    log_prior_alpha + 
    log_prior_beta + 
    log_prior_tau
  
  if(is.nan(log_posterior)) return(-Inf)
  log_posterior
}


# Metropolis algorithm
mh_accept = function(log_pi_x_star, log_pi_x_t) {
  met_ratio = log_pi_x_star - log_pi_x_t
  log(runif(1)) < met_ratio
}


# create two vectors to store the samples of mu and tau
alpha.post = vector(length=its)
beta.post = vector(length=its)
tau.post = vector(length=its)

# Samples initial values from the prior
alpha.post[1] = 0 #rnorm(1,alpha0,1/tau_alpha)
beta.post[1] = 0  #  rnorm(1,beta0,1/tau_beta)
tau.post[1] = 0.1 # rgamma(1,a,b)

# MCMC algorithm starts here
accept = list(alpha=0, beta=0, tau=0)

param = list(alpha = alpha.post[1],
             beta = beta.post[1],
             tau = tau.post[1])
for(i in 2:its) {
  
  # Draw alpha
  alpha_star = rnorm(1, mean=param$alpha, sd=1.0) # Proposal
  pi_star = joint_log_posterior(alpha_star, param$beta, param$tau)
  pi_t = joint_log_posterior(param$alpha, param$beta, param$tau)
  
  is_accepted = mh_accept(pi_star, pi_t)
  if(is_accepted) {
    param$alpha = alpha_star
    accept$alpha = accept$alpha + 1
  }
  
  # Draw beta
  beta_star = rnorm(1, mean=param$beta, sd=0.3) # Proposal
  pi_star = joint_log_posterior(param$alpha, beta_star, param$tau)
  pi_t = joint_log_posterior(param$alpha, param$beta, param$tau)
  
  is_accepted = mh_accept(pi_star, pi_t)
  if(is_accepted) {
    param$beta = beta_star
    accept$beta = accept$beta + 1
  }

  # Draw beta
  tau_star = rnorm(1, mean=param$tau, sd=0.4) # Proposal
  pi_star = joint_log_posterior(param$alpha, param$beta, tau_star)
  pi_t = joint_log_posterior(param$alpha, param$beta, param$tau)
  
  is_accepted = mh_accept(pi_star, pi_t)
  if(is_accepted) {
    param$tau = tau_star
    accept$tau = accept$tau + 1
  }
  
  # Record the new state of the chain
  alpha.post[i] = param$alpha
  beta.post[i] = param$beta
  tau.post[i] = param$tau
}

# Acceptance rates
cat("alpha acceptance:", accept$alpha/(its-1))
cat("beta acceptance:", accept$beta/(its-1))
cat("tau acceptance:", accept$tau/(its-1))

# Traceplots
par(mfrow=c(3, 1))
plot(alpha.post, ylab=expression(alpha), type='l')
plot(beta.post, ylab=expression(beta), type='l')
plot(tau.post, ylab=expression(tau), type='l')

# Plot beta against alpha
plot(alpha.post, beta.post, pch='.')

# OPTIONAL
# discard first 100 for burn in
alpha.post = mu.post[-c(1:100)]
beta.post = beta.post[-c(1:100)]
tau.post = tau.post[-c(1:100)]

par(mfrow=c(2, 1))
plot(mu.post, ylab=expression(mu), type='l')
plot(tau.post, ylab=expression(tau), type='l')

# ACF -- autocorrelation function
acf(mu.post)
acf(tau.post)
