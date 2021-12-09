import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


n = 100
mu0 = 10
tau0 = 100
alpha = 2
beta = 3
y = np.random.normal(mu, 1/tau)
mupostsamples = np.random.normal(mu0,1/tau0)
taupostsamples = np.random.gamma(alpha,beta)

for i in range(2, n):
    mupostsamples = np.random.normal(tau0*mu0+sum(y)*taupostsamples/tau0+n*taupostsamples, 1/tau0+n*taupostsamples)
    taupostsamples = np.random.gamma(alpha+n/2, beta+0.5*sum((y-mupostsamples)^2))
    i =+ 1

### my problem is the for loop i think ###
