# Testing modeling STRs with msprime

```bash
conda activate msprime-env
# conda install numpy pandas scipy
cd /home/amann11/2025-modeling_strs
ipython3
```
	
### 1. Shared parameters across all models

```python
import msprime
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.stats import norm

np.random.seed(548) # for reproducibility

# shared demographic parameters
gen_time = 25
original_n = 1_000 # original pop size

# str parameters
lo = 5 # min number of repeats
hi = 20 # max number of repeats

num_loci = 20 # number of STRs to simulate

# str mutation rate model
mean_log = -7 # meant rate of ~1e-3
sigma_log = 1 # ranges from 10-4 to 10-2 (typical range of mutations for Y-chromosome STRs)
str_mut_rate = np.random.lognormal(mean_log, sigma_log, num_loci) # array of mutation rates
str_mut_rate = np.clip(str_mut_rate, 1e-4, 5e-3) # max rate of 5e-3 
positions = np.arange(num_loci + 1)
rate_map = msprime.RateMap(position=positions, rate=str_mut_rate)
print("RateMap:", rate_map)

# get normalized distribution of repeats
x = np.arange(lo, hi+1)
probs = norm.pdf(x, loc=12, scale=2)  # mean=12, sd=2
root_dist = list(probs/probs.sum())   # normalized
```

### 2. Define each model

```python
# single founder event, slow pop growth
def single_founder():
	demography = msprime.Demography()
	demography.add_population(name="A", initial_size=original_n)

	# founder event 900 years ago (36 generations)
	founder_time = 900/gen_time
	demography.add_population_parameters_change(
		time=founder_time,
		initial_size=50, # size of the founder population
		growth_rate=np.log(original_n/50)/(900/gen_time), # slow recovery after inital founder event
		population="A")

	# simulate ancestry
	ts = msprime.sim_ancestry(
		samples={"A": 100},
		sequence_length=num_loci,
		demography=demography,
		ploidy=1,
		random_seed=245)

	return ts

# two founder lineages with slow growth
def two_founder():
	demography=msprime.Demography()
	demography.add_population(name="ancestral", initial_size=10_000)
	demography.add_population(name="A", initial_size=original_n)
	demography.add_population(name="B", initial_size=original_n)

	# populations split 500 years ago
	split_time = 500/gen_time
	demography.add_population_split(time=split_time, derived=["A", "B"], ancestral="ancestral")

	# independent founder events at the same time
	founder_time = 900/gen_time
	for pop in ["A", "B"]:
		demography.add_population_parameters_change(
			time=founder_time,
			initial_size=30,
			growth_rate=np.log(original_n/30)/(900/gen_time),
			population=pop)

	# sample from both populations
	ts = msprime.sim_ancestry(
		samples={"A": 50, "B": 50},
		sequence_length=num_loci,
		demography=demography,
		ploidy=1,
		random_seed=245)

	return ts

# single founder + later migrations
def migration():
    demography = msprime.Demography()
    demography.add_population(name="A", initial_size=original_n)
    demography.add_population(name="B", initial_size=original_n)  
    
    # founder event in population A
    founder_time = 900/gen_time
    demography.add_population_parameters_change(
        time=founder_time,
        initial_size=50,
        growth_rate=np.log(original_n/50)/(900/gen_time),
        population="A"
    )
    
    # migration from B to A starting 300 years ago
    migration_time = 300 / gen_time
    demography.set_migration_rate(source="B", dest="A", rate=1e-3)  # moderate gene flow between groups
    
    ts = msprime.sim_ancestry(
        samples={"A": 100},  # population B is unobserved
        demography=demography,
        sequence_length=num_loci,
        ploidy=1,
        random_seed=245
    )

    return ts	
```

### 3. Run models

```python
# run models and add mutations
for i, model in enumerate([single_founder, two_founder, migration], 1):
	ts = model()
	ts = msprime.sim_mutations(
		ts,
		rate=rate_map,
		model=msprime.SMM(
			lo=lo,
			hi=hi,
			root_distribution=root_dist),
		random_seed=312 + i)

# save haplotypes in repeat number format
	haplotypes = ts.genotype_matrix()
    df = pd.DataFrame(
        haplotypes.T,
        columns=[f"STR_{i}" for i in range(ts.num_sites)],
        index=[f"Sample_{i}" for i in ts.samples()]
    )
    df.to_csv(f"model_{i}_haplotypes.tsv", sep="\t")
    print(f"Saved Model {i} ({df.shape[0]} samples)")
```

















# add in stepwise mutations
ts_str = msprime.sim_mutations(
    ts,
    rate=rate_map,  # per-generation, per-locus mutation rate
    model=msprime.SMM(
        lo=lo,
        hi=hi,
        root_distribution=root_dist 
    ),
    random_seed=312
)

# extract first 20 str loci for each individual
str_loc = []
num_loc = 19

# create dict to store genotypes
dat = {}
for sample in ts_str.samples():
	sample_data = {}
	for i, variant in enumerate(ts_str.variants()):
		if i >= num_loc:
			break
		pos = variant.position
		sample_data[pos] = int(variant.alleles[variant.genotypes[sample]])
		dat[f"Sample_{sample}"] = sample_data

# convert to dataframe
df = pd.DataFrame.from_dict(dat, orient='index')
# clean up headers
df = df.add_prefix("STR_")

# get min and max of each str loci
df.describe().loc[['min', 'max']]

# print to file
df.to_csv("str_genotypes-simple.txt", index=True, sep="\t")
```

### 2. Two founder lineages with slow growth

```python


```

### 3. Single founder event and later migration from distinct paternal lineage

