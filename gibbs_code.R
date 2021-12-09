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

# create two vectors to store the samples of mu and tau
mu.post = vector(length=its)
tau.post = vector(length=its)

# Samples initial values from the prior
mu.post[1] = rnorm(1,mu0,1/tau0)
tau.post[1] = rgamma(1,alpha,beta)

# Create conditional posterior samplers
draw_mu = function(tau) {
  rnorm(1,
        (tau*sum(y)+tau0*mu0)/(N*tau+tau0),
        sd=sqrt(1/(N*tau+tau0)))
}

draw_tau = function(mu) {
  rgamma(1,alpha+N/2,beta+0.5*sum((y-mu)^2))
}

# sampling from each conditional
for(i in 2:its){
  mu.post[i] = draw_mu(tau.post[i-1])
  tau.post[i] = draw_tau(mu.post[i])
}

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

# "Pairplot"
par(mfrow=c(1,1))
plot(mu.post, tau.post, xlab=expression(mu), ylab=expression(tau), type='l')

# posterior density plots
plot(density(mu.post))
plot(density(1/tau.post))
plot(density(y))

# ACF -- autocorrelation function
acf(mu.post)
acf(tau.post)
