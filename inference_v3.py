import sys
import msprime
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

from migration2 import one_bin
from migration2 import two_bin

delta = 100
data = np.array([313.406833815]) #variance in pi between
total_sims = 1e6

def one_param():
    #inference procedure for 1 parameter
    # uses a uniform prior. Fagundes et al use a log-uniform prior for migration rates though
    param0 = []
    num_sims = 0
    while num_sims < total_sims:
        # draw parameter from uniform dist
        mig = np.random.uniform(0.01, 0.05)
        # simulate
        sim_stat = one_bin(500,250,250,10000, mig)
        # record parameter if it is close enough
        if np.absolute(sim_stat - data) < delta:
            param0.append(mig)
        num_sims = num_sims + 1
        if num_sims % 1000 == 0:
            print("Finished", num_sims, "simulations", file=sys.stderr)
    
    print("Number of accepted simulations:", len(param0))
    print("Probability of model:", len(param0)/total_sims)
    print("Median estimate of param0:", np.median(param0))
    print("95% confidence interval:", st.t.interval(0.95, len(param0)-1, loc=np.mean(param0), scale=st.sem(param0)))
    prob = len(param0)/total_sims
    aic = 2 - 2*np.log(prob)
    print("AIC:", aic)
    fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
    plt.hist(param0)
    plt.xlabel('Parameter estimate ')
    fig.savefig("1param.png")   # save the figure to file
    plt.close(fig)    # close the figure

def two_params():
    #inference procedure for 2 parameters
    # uses a uniform prior. Fagundes et al use a log-uniform prior for migration rates though
    param0 = []
    param1 = []
    num_sims = 0
    while num_sims < total_sims:
        # draw parameter from uniform dist
        mig0 = np.random.uniform(0.01, 0.05)
        mig1 = np.random.uniform(0.01, 0.05)
        # simulate
        sim_stat = two_bin(500,250,250,10000, mig0, mig1)
        # record parameter if it is close enough
        if np.absolute(sim_stat - data) < delta:
            param0.append(mig0)
            param1.append(mig1)
        num_sims = num_sims + 1
        if num_sims % 1000 == 0:
            print("Finished", num_sims, "simulations", file=sys.stderr)
    
    print("Number of accepted simulations:", len(param0))
    print("Probability of model:", len(param0)/total_sims)
    print("Median estimate of param0:", np.median(param0))
    print("95% HPD:", st.t.interval(0.95, len(param0)-1, loc=np.mean(param0), scale=st.sem(param0)))
    print("Median estimate of param1:", np.median(param1))
    print("95% HPD:", st.t.interval(0.95, len(param1)-1, loc=np.mean(param1), scale=st.sem(param1)))
    prob = len(param0)/total_sims
    aic = 2 - 2*np.log(prob)
    print("AIC:", aic)
    fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
    plt.hist(param0)
    plt.hist(param1)
    plt.xlabel('Parameter estimate')
    fig.savefig("2params.png")   # save the figure to file
    plt.close(fig)    # close the figure

if __name__ == "__main__":
    one_param()
    two_params()

