# try inference based on # segregating sites and pi and variances
import msprime
import math
import numpy as np
from migration2 import one_bin
from migration2 import two_bin

# data comes from migration2.py simulations of low_high model
data = np.array([324.467523549])

# try to match with one time bin
deltas = []
#print("Time bins (k) = 1")
#for i in np.arange(0.01,0.05, 0.01):
#    print("Testing parameter:",i)
#    result = one_bin(500,250,250,10000,i)
#    deltas.append([i,np.sum((np.power((result-data),2)/result))])
#print("mig delta")
#for d in deltas:
#    print(d[0],d[1])

#----
# try to match with two time bins
# 8,000,000 simulations
deltas = []
print("Time bins (k) = 2")
for i in np.arange(0.01,0.05, 0.01):
    for j in np.arange(0.01,0.05, 0.01):
        print("Testing parameters:",i,j)
        result = two_bin(500,250,250,10000,i, j)
        deltas.append([i, j, np.sum((np.power((result-data),2)/result))])
print("mig1 mig2 delta")
for d in deltas:
    print(d[0],d[1],d[2])

#----
# try to match with three time bins
# 32,000,000 simulations
#deltas = []
#for i in np.arange(0.01,0.05, 0.01):
#    for j in np.arange(0.01,0.05, 0.01):
#        for k in np.arange(0.01,0.05, 0.01):
#            print(i,j,k)
#            result = three_bins(500,500,500,10000,i,j,k)
#            deltas.append([i, j, k, np.sum(np.power((result-data),2))])
#print("Three time bins (migration parameter (k=1), migration parameter (k=2), migration parameter (k=3), delta_variance)")
#for d in deltas:
#    print(d)
