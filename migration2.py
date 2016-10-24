import msprime
import math
import numpy as np

def pairwise_diffs(hap0, hap1):
    return sum(1 for a, b in zip(hap0, hap1) if a != b)

def one_bin(NA,N1,N2,Ts,M):
    NA=NA
    N1=N1
    N2=N2
    Ts=Ts
    M_ave=M
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=50, initial_size=N1),
        msprime.PopulationConfiguration(sample_size=50, initial_size=N2)
    ]
    migration_matrix = [
        [0, M_ave],
        [M_ave, 0]
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
    
    
    replicates=1
    length=100000
    sim = msprime.simulate(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=1e-7,
        recombination_rate=1e-8,
        length=length, 
        num_replicates=replicates)
        
    pairwise_diff = []
    for j, s in enumerate(sim):
        s0 = len(s.get_samples(0))
        s1 = len(s.get_samples(1))
        haps = [h for h in s.haplotypes()]
        h0 = haps[0:s0]
        h1 = haps[s0:s0+s1-1]
        
        for hap0 in h0:
            for hap1 in h1:
                pairwise_diff.append(sum(1 for a, b in zip(hap0, hap1) if a != b))
        #pairwise_diff.append([pairwise_diffs(hap0,hap1) for hap0 in h0 for hap1 in h1])
            
    return(np.var(np.array(pairwise_diff)))
    
 
def two_bin(NA,N1,N2,Ts,M1,M2):
    NA=NA
    N1=N1
    N2=N2
    Ts=Ts
    M1=M1
    M2=M2
    
    population_configurations = [
        msprime.PopulationConfiguration(sample_size=50, initial_size=N1),
        msprime.PopulationConfiguration(sample_size=50, initial_size=N2)
    ]
    migration_matrix = [
        [0, M1],
        [M1, 0]
    ]
    demographic_events = [
        msprime.MigrationRateChange(time=Ts/2, rate=M2, matrix_index=(0, 1)),
        msprime.MigrationRateChange(time=Ts/2, rate=M2, matrix_index=(1, 0)),
        msprime.MassMigration(time=Ts, source=1, destination=0, proportion=1.0)
    ]
    
    dp = msprime.DemographyDebugger(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events)
    #dp.print_history()
    replicates=1
    length=100000
    sim = msprime.simulate(
        Ne=NA,         
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        mutation_rate=1e-7,
        recombination_rate=1e-8,
        length=length, 
        num_replicates=replicates)
    pairwise_diff = []
    for j, s in enumerate(sim):
        s0 = len(s.get_samples(0))
        s1 = len(s.get_samples(1))
        haps = [h for h in s.haplotypes()]
        h0 = haps[0:s0]
        h1 = haps[s0:s0+s1-1]
        for hap0 in h0:
            for hap1 in h1:
                pairwise_diff.append(sum(1 for a, b in zip(hap0, hap1) if a != b))
    return(np.var(np.array(pairwise_diff)))
if __name__ == "__main__":    
    print("constant")
    print(one_bin(500,250,250,10000,0.02))
    #print("high low")
    #print(two_bin(500,250,250,10000,0.03,0.01))
    #print("low high")
    #print(two_bin(500,250,250,10000,0.01,0.03))
