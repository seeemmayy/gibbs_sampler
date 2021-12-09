# Model
# y_i ~ Normal(mu, 1/tau)
# 
# Priors
# mu ~ Normal(mu0, 1/tau0)
# tau ~ Gamma(alpha, beta)

# Set "true" values of 
N = 100
mu = 4
tau = 4

# hyperparameters
alpha = 1
beta = 2
mu0 = 2
tau0 = 2

# Simulate some data
y = rnorm(N,mu,sd=sqrt(1/tau))

# Run the MCMC for `its` iterations
its=1000


# Create conditional posterior density functions
pi_mu = function(mu, tau) {
  dnorm(mu,
        (tau*sum(y)+tau0*mu0)/(N*tau+tau0),
        sd=sqrt(1/(N*tau+tau0)))
}

pi_tau = function(tau, mu) {
  dgamma(tau,
         alpha+N/2,
         beta+0.5*sum((y-mu)^2))
}


# Metropolis algorithm
mh_accept = function(pi_x_star, pi_x_t) {
  met_ratio = pi_x_star / pi_x_t
  runif(1) < met_ratio
}


# create two vectors to store the samples of mu and tau
mu.post = vector(length=its)
tau.post = vector(length=its)

# Samples initial values from the prior
mu.post[1] = rnorm(1,mu0,1/tau0)
tau.post[1] = rgamma(1,alpha,beta)

# MCMC algorithm starts here
mu_accept = 0
tau_accept = 0
for(i in 2:its) {
  
  # Draw mu
  mu_star = rnorm(1, mean=mu.post[i-1], sd=0.15) # Proposal
  pi_mu_star = pi_mu(mu_star, tau.post[i-1])
  pi_mu_t = pi_mu(mu.post[i-1], tau.post[i-1])
  
  is_accepted = mh_accept(pi_mu_star, pi_mu_t)
  if(is_accepted) {
    mu.post[i] = mu_star
    mu_accept = mu_accept + 1
  }
  else {
    mu.post[i] = mu.post[i-1]
  }
  
  # Draw tau
  tau_star = rnorm(1, mean=tau.post[i-1], sd=1.3) # Proposal
  pi_tau_star = pi_tau(tau_star, mu.post[i])
  pi_tau_t = pi_tau(tau.post[i-1], mu.post[i])
  
  is_accepted = mh_accept(pi_tau_star, pi_tau_t)
  if(is_accepted) {
    tau.post[i] = tau_star
    tau_accept = tau_accept + 1
  }
  else {
    tau.post[i] = tau.post[i-1]
  }
}

# Acceptance rates
cat("Mu acceptance:", mu_accept/(its-1))
cat("tau acceptance:", tau_accept/(its-1))

# Traceplots
par(mfrow=c(2, 1))
plot(mu.post, ylab=expression(mu), type='l')
plot(tau.post, ylab=expression(tau), type='l')

# OPTIONAL
# discard first 100 for burn in
mu.post <- mu.post[-c(1:100)]
tau.post <- tau.post[-c(1:100)]

par(mfrow=c(2, 1))
plot(mu.post, ylab=expression(mu), type='l')
plot(tau.post, ylab=expression(tau), type='l')

# ACF -- autocorrelation function
acf(mu.post)
acf(tau.post)
