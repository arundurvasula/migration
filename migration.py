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
    
    dp = msprime.DemographyDebugger(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    #dp.print_history()

    sim = msprime.simulate(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=1e-7, 
        recombination_rate=1e-8,
        length=100000, 
        num_replicates=1000000)
    pi = []
    seg = []
    for s in sim:
        pi.append(s.get_pairwise_diversity())
        seg.append(s.get_num_mutations())

    print(np.mean(pi))
    print(np.var(pi))
    print(np.mean(seg))
    print(np.var(seg))
 
def two_bin(NA,N1,N2,Ts,M1,M2):
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
    
    dp = msprime.DemographyDebugger(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    #dp.print_history()

    sim = msprime.simulate(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=1e-7,
        recombination_rate=1e-8,
        length=100000, 
        num_replicates=1000000)
    pi = []
    seg = []
    for s in sim:
        pi.append(s.get_pairwise_diversity())
        seg.append(s.get_num_mutations())

    print(np.mean(pi))
    print(np.var(pi))
    print(np.mean(seg))
    print(np.var(seg))  
    
print("constant")
#one_bin(500,500,500,10000,0.02)
print("high low")
#two_bin(500,500,500,10000,0.03,0.01)
print("low high")
two_bin(500,500,500,10000,0.01,0.03)
