# try inference based on # segregating sites and pi and variances
import msprime
import math
import numpy as np

def one_bin(NA,N1,N2,Ts,M):
    NA=NA
    N1=N1
    N2=N2
    Ts=Ts
    M_ave=M
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=0, initial_size=N1),
        msprime.PopulationConfiguration(sample_size=50, initial_size=N2)
    ]
    migration_matrix = [
        [0, M_ave],
        [0, 0]
    ]
    demographic_events = [
        msprime.MassMigration(time=Ts, source=1, destination=0, proportion=1.0)
    ]
    
    #dp = msprime.DemographyDebugger(
    #    Ne=NA,         
    #    population_configurations=population_configurations,
    #    migration_matrix=migration_matrix,
    #    demographic_events=demographic_events)
    #dp.print_history()

    sim = msprime.simulate(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=1e-7, 
        recombination_rate=1e-8,
        length=100000, 
        num_replicates=500000)
    pi = []
    seg = []
    for s in sim:
        pi.append(s.get_pairwise_diversity())
        seg.append(s.get_num_mutations())

    return(np.array([np.mean(pi),np.var(pi),np.mean(seg),np.var(seg)]))

def two_bins(NA,N1,N2,Ts,M1,M2):
    NA=NA
    N1=N1
    N2=N2
    Ts=Ts
    M1=M1
    M2=M2
    
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=0, initial_size=N1),
        msprime.PopulationConfiguration(sample_size=50, initial_size=N2)
    ]
    migration_matrix = [
        [0, M2],
        [0, 0]
    ]
    demographic_events = [
        msprime.MigrationRateChange(time=Ts/2, rate=M1, matrix_index=(0, 1)),
        #msprime.MigrationRateChange(time=Ts/2, rate=M1, matrix_index=(1, 0)),
        msprime.MassMigration(time=Ts, source=1, destination=0, proportion=1.0)
    ]
    
    #dp = msprime.DemographyDebugger(
    #    Ne=NA,         
    #    population_configurations=population_configurations,
    #    migration_matrix=migration_matrix,
    #    demographic_events=demographic_events)
    #dp.print_history()

    sim = msprime.simulate(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=1e-7,
        recombination_rate=1e-8,
        length=100000, 
        num_replicates=500000)
    pi = []
    seg = []
    for s in sim:
        pi.append(s.get_pairwise_diversity())
        seg.append(s.get_num_mutations())

    return(np.array([np.mean(pi),np.var(pi),np.mean(seg),np.var(seg)]))

def three_bins(NA,N1,N2,Ts,M1,M2,M3):
    NA=NA
    N1=N1
    N2=N2
    Ts=Ts
    M1=M1
    M2=M2
    M3=M3
    
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=0, initial_size=N1),
        msprime.PopulationConfiguration(sample_size=50, initial_size=N2)
    ]
    migration_matrix = [
        [0, M2],
        [0, 0]
    ]
    demographic_events = [
        msprime.MigrationRateChange(time=Ts/3, rate=M1, matrix_index=(0, 1)),
        msprime.MigrationRateChange(time=2*(Ts/3), rate=M1, matrix_index=(0, 1)),
        msprime.MassMigration(time=Ts, source=1, destination=0, proportion=1.0)
    ]
    
    #dp = msprime.DemographyDebugger(
    #    Ne=NA,         
    #    population_configurations=population_configurations,
    #    migration_matrix=migration_matrix,
    #    demographic_events=demographic_events)
    #dp.print_history()

    sim = msprime.simulate(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=1e-7,
        recombination_rate=1e-8,
        length=100000, 
        num_replicates=500000)
    pi = []
    seg = []
    for s in sim:
        pi.append(s.get_pairwise_diversity())
        seg.append(s.get_num_mutations())

    return(np.array([np.mean(pi),np.var(pi),np.mean(seg),np.var(seg)]))


# data comes from migration.py simulations of low_high model
data = np.array([19.9988318604,71.4938202807,89.598257,553.867747562])

# try to match with one time bin
#deltas = []
#for i in np.arange(0.01,0.05, 0.01):
#    print(i)
#    result = one_bin(500,500,500,10000,i)
#    deltas.append([i,np.sum(np.power((result-data),2))])
#print("One time bin (migration parameter, delta_variance)")
#for d in deltas:
#    print(d)

#----
# try to match with two time bins
# 8,000,000 simulations
#deltas = []
#for i in np.arange(0.01,0.05, 0.01):
#    for j in np.arange(0.01,0.05, 0.01):
#        print(i,j)
#        result = two_bins(500,500,500,10000,i, j)
#        deltas.append([i, j, np.sum(np.power((result-data),2))])
#print("Two time bins (migration parameter (k=1), migration parameter (k=2), delta_variance)")
#for d in deltas:
#    print(d)

#----
# try to match with three time bins
# 32,000,000 simulations
deltas = []
for i in np.arange(0.01,0.05, 0.01):
    for j in np.arange(0.01,0.05, 0.01):
        for k in np.arange(0.01,0.05, 0.01):
            print(i,j,k)
            result = three_bins(500,500,500,10000,i,j,k)
            deltas.append([i, j, k, np.sum(np.power((result-data),2))])
print("Three time bins (migration parameter (k=1), migration parameter (k=2), migration parameter (k=3), delta_variance)")
for d in deltas:
    print(d)
