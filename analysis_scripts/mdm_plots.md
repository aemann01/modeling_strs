# Modeling STR demographic history project

### 1. First run replicates of our simple, single founder event model

```bash
# running replicates from command line
cd /home/amann11/2025-modeling_strs/scripts
# running single migration, 10 patrilines founder group
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 01-single_migration-10pat.slim 1> 01-out.txt
# single migration, 50 patrilines founder group
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 02-single_migration-50pat.slim 1> 02-out.txt
# single migration, 100 patrilines founder group
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 03-single_migration-100pat.slim 1> 03-out.txt
```

Move them into the results file -- need to put this in the script at one point

```bash
mv single_migration_sim_* ../results
cd ../results
```

### 2. 10 patrilines

Load required libraries

```R
library(tidyverse)
library(dplyr)
library(ggplot2)
```

Pull all replicate files into R, merge and add replicate identifier

```R
# get list of files
file_list <- list.files(pattern = "single_migration_sim_10pat_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)
```

Calculate MDM for each replicate (modified from Shiny app from Mann et al. 2024)

```R
# get modal haplotype function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# calculate mdm for each replicate and each haplogroup function
calculate_mdm <- function(data, hg_filter) {
  # filter by haplogroup
  hg_data <- data %>% filter(HG == hg_filter)
  
  # unique combinations of replicate + Pop
  groups <- hg_data %>% distinct(replicate, Pop)
  
  mdm_results <- map_dfr(1:nrow(groups), function(i) {
    rep_num <- groups$replicate[i]
    pop_val <- groups$Pop[i]
    
    # Filter to this replicate & population
    rep_data <- hg_data %>%
      filter(replicate == rep_num, Pop == pop_val) %>%
      select(starts_with("locus"))
    
    if(nrow(rep_data) == 0) return(NULL)  # skip empty groups
    
    # calculate modal haplotype
    mode_hap <- sapply(rep_data, getmode)
    
    # calculate mdm for each individual
    mdm_values <- c()
    for(row in 1:nrow(rep_data)) {
      temp <- rbind(rep_data[row,], mode_hap)
      temp <- temp - as.list(temp[1,])
      mdm_values <- c(mdm_values, sum(abs(temp[2,])))
    }
    
    data.frame(
      replicate = rep_num,
      Pop = pop_val,
      mdm = mdm_values
    )
  })
  
  return(mdm_results)
}
```

Calculate for each haplogroup

```R
mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )
```

Now calculate frequency and mean

```R
freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.481    0.313     0.552
# 2 source    0.0662   0.101     0.00237

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.483    0.294       0.543
# 2 source    0.0829   0.123       0
```

And plot

```R
pop_colors <- c("found" = "#d95f02", "source" = "#1b9e77")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_10pat_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_10pat_HG2.pdf")
p
dev.off()
```

### 3. 50 patrilines

```R
# get list of files
file_list <- list.files(pattern = "single_migration_sim_50pat_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

# calculate mdm for each haplogroup
mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.305    0.228     0.317
# 2 source    0.0890   0.131     0.00840

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found      0.275   0.236      0.289
# 2 source     0.103   0.128      0.0261

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_50pat_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_50pat_HG2.pdf")
p
dev.off()
```

### 4. 100 patrilines

```R
# get list of files
file_list <- list.files(pattern = "single_migration_sim_100pat_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

# calculate mdm for each haplogroup
mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found      0.287   0.259       0.286
# 2 source     0.182   0.216       0.120

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found      0.208   0.246      0.123
# 2 source     0.173   0.236      0.0284

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_100pat_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_100pat_HG2.pdf")
p
dev.off()
```

### 5. Running replicates for migration from source population

```bash
# running replicates from command line
cd /home/amann11/2025-modeling_strs/scripts
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 04-source_migration-low.slim 1> 04-out.txt
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 05-source_migration-med.slim 1> 05-out.txt
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 06-source_migration-high.slim 1> 06-out.txt
```

Move them into the results file -- need to put this in the script at one point

```bash
mv source_migration_sim_* ../results
cd ../results
```

### 6. Low migration rate

```R    
library(tidyverse)
library(dplyr)
library(ggplot2)

# get list of files
file_list <- list.files(pattern = "source_migration_sim_low_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

# get modal haplotype function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# calculate mdm for each replicate and each haplogroup function
calculate_mdm <- function(data, hg_filter) {
  # filter by haplogroup
  hg_data <- data %>% filter(HG == hg_filter)
  
  # unique combinations of replicate + Pop
  groups <- hg_data %>% distinct(replicate, Pop)
  
  mdm_results <- map_dfr(1:nrow(groups), function(i) {
    rep_num <- groups$replicate[i]
    pop_val <- groups$Pop[i]
    
    # Filter to this replicate & population
    rep_data <- hg_data %>%
      filter(replicate == rep_num, Pop == pop_val) %>%
      select(starts_with("locus"))
    
    if(nrow(rep_data) == 0) return(NULL)  # skip empty groups
    
    # calculate modal haplotype
    mode_hap <- sapply(rep_data, getmode)
    
    # calculate mdm for each individual
    mdm_values <- c()
    for(row in 1:nrow(rep_data)) {
      temp <- rbind(rep_data[row,], mode_hap)
      temp <- temp - as.list(temp[1,])
      mdm_values <- c(mdm_values, sum(abs(temp[2,])))
    }
    
    data.frame(
      replicate = rep_num,
      Pop = pop_val,
      mdm = mdm_values
    )
  })
  
  return(mdm_results)
}

mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.153    0.173      0.0773
# 2 source    0.0881   0.113      0.0155

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.149    0.163      0.0617
# 2 source    0.0842   0.131      0

pop_colors <- c("found" = "#d95f02", "source" = "#1b9e77")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_lowSource_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_lowSource_HG2.pdf")
p
dev.off()
```

### 7. Medium migration rate

```R    
# get list of files
file_list <- list.files(pattern = "source_migration_sim_med_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.101    0.136     0.00247
# 2 source    0.0847   0.128     0

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.0886   0.109     0.0217
# 2 source    0.0799   0.110     0.00180

pop_colors <- c("found" = "#d95f02", "source" = "#1b9e77")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_medSource_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_medSource_HG2.pdf")
p
dev.off()
```

### 8. High migration rate (reaching total population replacement)

```R    
# get list of files
file_list <- list.files(pattern = "source_migration_sim_high_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.107    0.145     0.00597
# 2 source    0.0995   0.140     0.00447

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.0960   0.125     0.00519
# 2 source    0.0945   0.124     0.00376

pop_colors <- c("found" = "#d95f02", "source" = "#1b9e77")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_highSource_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_highSource_HG2.pdf")
p
dev.off()
```

### 9. Running replicates for migration from divergent population

```bash
# running replicates from command line
cd /home/amann11/2025-modeling_strs/scripts
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 07-divergent_source_migration-low.slim 1> 07-out.txt
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 08-divergent_source_migration-med.slim 1> 08-out.txt
seq 1 100 | parallel -j20 slim -d "replicate={}" -d "seed={#}" 09-divergent_source_migration-high.slim 1> 09-out.txt
```

Move them into the results file -- need to put this in the script at one point

```bash
mv divergent_migration_sim_* ../results
cd ../results
```






### 10. Low migration rate divergent pop

```R    
library(tidyverse)
library(dplyr)
library(ggplot2)

# get list of files
file_list <- list.files(pattern = "divergent_migration_sim_low_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

# get modal haplotype function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# calculate mdm for each replicate and each haplogroup function
calculate_mdm <- function(data, hg_filter) {
  # filter by haplogroup
  hg_data <- data %>% filter(HG == hg_filter)
  
  # unique combinations of replicate + Pop
  groups <- hg_data %>% distinct(replicate, Pop)
  
  mdm_results <- map_dfr(1:nrow(groups), function(i) {
    rep_num <- groups$replicate[i]
    pop_val <- groups$Pop[i]
    
    # Filter to this replicate & population
    rep_data <- hg_data %>%
      filter(replicate == rep_num, Pop == pop_val) %>%
      select(starts_with("locus"))
    
    if(nrow(rep_data) == 0) return(NULL)  # skip empty groups
    
    # calculate modal haplotype
    mode_hap <- sapply(rep_data, getmode)
    
    # calculate mdm for each individual
    mdm_values <- c()
    for(row in 1:nrow(rep_data)) {
      temp <- rbind(rep_data[row,], mode_hap)
      temp <- temp - as.list(temp[1,])
      mdm_values <- c(mdm_values, sum(abs(temp[2,])))
    }
    
    data.frame(
      replicate = rep_num,
      Pop = pop_val,
      mdm = mdm_values
    )
  })
  
  return(mdm_results)
}

mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.153    0.173      0.0773
# 2 source    0.0881   0.113      0.0155

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.149    0.163      0.0617
# 2 source    0.0842   0.131      0

pop_colors <- c("found" = "#d95f02", "source" = "#1b9e77")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_lowSource_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_lowSource_HG2.pdf")
p
dev.off()
```

### 7. Medium migration rate

```R    
# get list of files
file_list <- list.files(pattern = "source_migration_sim_med_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.101    0.136     0.00247
# 2 source    0.0847   0.128     0

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.0886   0.109     0.0217
# 2 source    0.0799   0.110     0.00180

pop_colors <- c("found" = "#d95f02", "source" = "#1b9e77")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_medSource_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_medSource_HG2.pdf")
p
dev.off()
```

### 8. High migration rate (reaching total population replacement)

```R    
# get list of files
file_list <- list.files(pattern = "source_migration_sim_high_rep")
# read in and combine files
combined_data <- map_dfr(file_list, function(file) {
	rep_number <- str_extract(file, "rep(\\d+)", group = 1) # extract rep number from file name
	data <- read_csv(file, show_col_types = FALSE)
  
  # Add replicate column
  data <- data %>% mutate(replicate = as.numeric(rep_number))
  
  return(data)
})
dim(combined_data)

mdm_data_HG1 <- calculate_mdm(combined_data, "HG1")
mdm_data_HG2 <- calculate_mdm(combined_data, "HG2")

# sanity check that the replicates aren't clones
mdm_data_HG1 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

mdm_data_HG2 %>%
  group_by(replicate) %>%
  summarise(
    n_obs = n(),
    mean_mdm = mean(mdm),
    var_mdm  = var(mdm)
  )

freq_df <- mdm_data_HG1 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG1 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.107    0.145     0.00597
# 2 source    0.0995   0.140     0.00447

# calculate median proportion of neighbor haplotypes (mdm <= 1) per replicate
mdm_data_HG2 %>%
  group_by(replicate, Pop) %>%
  summarise(prop_neighbor = mean(mdm <= 1), .groups = "drop") %>%
  group_by(Pop) %>%                          # group by population
  summarise(
    mean_prop = mean(prop_neighbor),
    sd_prop = sd(prop_neighbor),
    median_prop = median(prop_neighbor),     # add median
    .groups = "drop"
  )
#   Pop    mean_prop sd_prop median_prop
#   <chr>      <dbl>   <dbl>       <dbl>
# 1 found     0.0960   0.125     0.00519
# 2 source    0.0945   0.124     0.00376

pop_colors <- c("found" = "#d95f02", "source" = "#1b9e77")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_highSource_HG1.pdf")
p
dev.off()

# do same for HG2
freq_df <- mdm_data_HG2 %>%
  count(replicate, Pop, mdm, name = "n") %>%
  group_by(replicate, Pop) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

mean_df <- freq_df %>%
  group_by(Pop, mdm) %>%
  summarise(mean_freq = mean(freq), .groups = "drop")

p <- ggplot() +
  geom_step(data = freq_df,
            aes(x = mdm, y = freq, group = interaction(Pop, replicate)),
            color = "grey70",
            alpha = 0.5,
            linewidth = 0.6) +
  
  # population mean lines colored
  geom_step(data = mean_df,
            aes(x = mdm, y = mean_freq, color = Pop),
            linewidth = 1.4) +
  
  scale_color_manual(values = pop_colors) +
  
  labs(
    x = "Number of repeat differences from modal haplotype",
    y = "Frequency",
    color = "Population"
  ) +
  ylim(0.0, 1.0) +
  theme_minimal()
pdf("mdm_highSource_HG2.pdf")
p
dev.off()
```

