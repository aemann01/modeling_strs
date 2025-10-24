# Y-chromosomal STR founder event simultaions

### 1. Set up environment

```bash
conda activate msprime-env
# conda install numpy pandas scipy
cd /home/amann11/2025-modeling_strs
ipython3
```

### 2. Load required libraries

```python
import msprime
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy.stats import norm
```

### 3. Set up STR mutation rates

```python
# str parameters
np.random.seed(548) # make sure mutations are the same each time
num_loci = 20 # number of STRs to simulate
mean_log = -7 # meant rate of ~1e-3
sigma_log = 1 # ranges from 10-4 to 10-2 (typical range of mutations for Y-chromosome STRs)
str_mut_rate = np.random.lognormal(mean_log, sigma_log, num_loci) # array of mutation rates
str_mut_rate = np.clip(str_mut_rate, 1e-4, 5e-3) # max rate of 5e-3 
positions = np.arange(num_loci + 1)

rate_map = msprime.RateMap(position=positions, rate=str_mut_rate)
print("RateMap:", rate_map)

# str parameters
lo = 5 # min number of repeats
hi = 20 # max number of repeats

# get normalized distribution of repeats
x = np.arange(lo, hi+1)
probs = norm.pdf(x, loc=12, scale=2)  # mean=12, sd=2
root_dist = list(probs/probs.sum())   # normalized
```

### 4. General population parameters

```python
# general parameters
gen_time = 25
found_time = 900 / gen_time # time at the founder event
original_n = 1_000 # original pop size
founder_n = 25 # founder pop size
recovery_time = 900 / gen_time # time to recover to original size

# growth rate (exponential recovery from founder to original size)
alpha = np.log(original_n / founder_n) / recovery_time
```

### 4. Simple founder event 

```python
# initialize empty demographic model
demography = msprime.Demography()

# define original population
demography.add_population(name="original_pop", description="Hypothetical mainland population", initial_size=original_n)
# founder event
demography.add_population_parameters_change(
	time=found_time,
	initial_size=founder_n,
	growth_rate=alpha,
	population="original_pop")

# simulate genology of haploid Y-chromosome
ts = msprime.sim_ancestry(
	samples=100,
	sequence_length=num_loci,
	demography=demography,
	ploidy=1,
	random_seed=245)

print(ts)

# check how the population has grown over time
growth = np.linspace(found_time, found_time + recovery_time, 5)
for x in growth:
	nt = founder_n * np.exp(alpha * (x - found_time))
	print(f"At generation={x:.0f}: N={nt:.1f} (target={original_n})")

# add in stepwise mutations
ts_str = msprime.sim_mutations(
    ts,
    rate=rate_map,  # per-generation, per-locus mutation rate
    model=msprime.SMM(
        lo=lo,
        hi=hi,
        root_distribution=root_dist 
    ),
    random_seed=312,
)

# check first locus to make sure mutations added
print(f"Simulated {ts_str.num_sites} STR loci")
print(f"First variant: {next(ts_str.variants())}")

## NOTE! The small population size of this founder simulation means that not all STR loci will have mutations (and therefore not be in the output) -- we can assume an ancestral state for these but they won't be useful in analysis

# extract first 20 str loci for each individual (only get 18 with this model)
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
df.columns = df.columns.str.replace(r'\.0$', '', regex=True)

# add in columns for MDM intput
df.insert(0, 'Pop', 'simple')
df['HG'] = 'R'
# get min and max of each str loci
df.describe().loc[['min', 'max']]

# print to file
df.to_csv("str_genotypes-simple.csv", index=True)
```

### 5. Founder event with two lineages

```python
# initialize empty demographic model
demography = msprime.Demography()

# define original populations
demography.add_population(name="ancestral", description="ancestral pop to founder lineages", initial_size=5_000)
demography.add_population(name="original_pop", description="Hypothetical mainland population", initial_size=original_n)
demography.add_population(name="original_pop2", description="Secondary mainland population", initial_size=original_n)

# split time between founder lineages
split_time = 50/gen_time
demography.add_population_split(time=split_time, derived=["original_pop", "original_pop2"], ancestral="ancestral")

# add in some asymmetric migration between populations 
migration_start_time = split_time + 5/gen_time  # 5 gens after split
migration_end_time = split_time + 10/gen_time   # lasts 10 generations

demography.set_migration_rate(source="original_pop", dest="original_pop2", rate=1e-4) 
demography.set_migration_rate(source="original_pop2", dest="original_pop", rate=1e-4)

demography.add_symmetric_migration_rate_change(
    time=migration_start_time,
    populations=["original_pop", "original_pop2"],
    rate=1e-4 
)

demography.add_symmetric_migration_rate_change(
    time=migration_end_time,
    populations=["original_pop", "original_pop2"],
    rate=0  # stops gene flow
)

# independent founder events at 900 years ago
for pop in ["original_pop", "original_pop2"]:
	demography.add_population_parameters_change(
		time=found_time,
		initial_size=founder_n,
		growth_rate=alpha,
		population=pop)

# simulate genology of haploid Y-chromosome
ts = msprime.sim_ancestry(
	samples={"original_pop": 50, "original_pop2": 50}, # 50 samples from each founder lineage
	sequence_length=num_loci,
	demography=demography,
	ploidy=1,
	random_seed=245)

print(ts)

# check how the population has grown over time
growth = np.linspace(found_time, found_time + recovery_time, 5)
for x in growth:
	nt = founder_n * np.exp(alpha * (x - found_time))
	print(f"At generation={x:.0f}: N={nt:.1f} (target={original_n})")

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

# check first locus to make sure mutations added
print(f"Simulated {ts_str.num_sites} STR loci")
print(f"First variant: {next(ts_str.variants())}")

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
df.columns = df.columns.str.replace(r'\.0$', '', regex=True)

# add in columns for MDM intput
df.insert(0, 'Pop', 'simple')
df['HG'] = 'R'
# get min and max of each str loci
df.describe().loc[['min', 'max']]

# print to file
df.to_csv("str_genotypes-twolineages.csv", index=True)
```

### 6. Single founder event with later migration for distinct lineage

```python
# initialize empty demographic model
demography = msprime.Demography()

# define original populations
demography.add_population(name="ancestral", description="ancestral pop to founder lineages", initial_size=5_000)
demography.add_population(name="original_pop", description="Hypothetical mainland population", initial_size=original_n)
demography.add_population(name="original_pop2", description="Secondary mainland population", initial_size=original_n)

# split time between founder lineages
split_time = 10/gen_time
demography.add_population_split(time=split_time, derived=["original_pop", "original_pop2"], ancestral="ancestral")

# add in some asymmetric migration between populations 
migration_start_time = split_time + 5/gen_time  # 5 gens after split
migration_end_time = split_time + 50/gen_time   # lasts 50 generations

demography.set_migration_rate(source="original_pop", dest="original_pop2", rate=1e-4) 
demography.set_migration_rate(source="original_pop2", dest="original_pop", rate=1e-4)

demography.add_symmetric_migration_rate_change(
    time=migration_start_time,
    populations=["original_pop", "original_pop2"],
    rate=1e-4 
)

demography.add_symmetric_migration_rate_change(
    time=migration_end_time,
    populations=["original_pop", "original_pop2"],
    rate=0  # stops gene flow
)

# initial founder event 900 years ago
demography.add_population_parameters_change(
	time=found_time,
	initial_size=founder_n,
	growth_rate=alpha,
	population="original_pop")

# migration of secondary population into first founder starting 200 years ago
migration_time = 200/gen_time
demography.set_migration_rate(source="original_pop2", dest="original_pop", rate=1e-3)

# simulate genology of haploid Y-chromosome
ts = msprime.sim_ancestry(
	samples={"original_pop": 100}, # sample from original founder population (secondary is unobserved)
	sequence_length=num_loci,
	demography=demography,
	ploidy=1,
	random_seed=245)

print(ts)

# check how the population has grown over time
growth = np.linspace(found_time, found_time + recovery_time, 5)
for x in growth:
	nt = founder_n * np.exp(alpha * (x - found_time))
	print(f"At generation={x:.0f}: N={nt:.1f} (target={original_n})")

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

# check first locus to make sure mutations added
print(f"Simulated {ts_str.num_sites} STR loci")
print(f"First variant: {next(ts_str.variants())}")

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
df.columns = df.columns.str.replace(r'\.0$', '', regex=True)

# add in columns for MDM intput
df.insert(0, 'Pop', 'simple')
df['HG'] = 'R'
# get min and max of each str loci
df.describe().loc[['min', 'max']]

# print to file
df.to_csv("str_genotypes-migration.csv", index=True)
```

