library(tidyverse)
library(readxl)
library(glue)

#### Get correspondance between individual ids and samples barcodes ####
# Read patients between wave 1
wave1 <-
  read_xlsx("data/raw/COVID PAT min 3 PAXGene 161120.xlsx", sheet = 2) %>% 
  #filter(!str_detect(IPP, "CHC")) %>% # Remove CHC patients due to the current lack of metadata; UPD: already have
  select(PATID, DATE = DPREL)

# Read patients between wave 2
wave2 <- 
  read_xlsx("data/raw/COVID PAT 2nd wave min 3 PAXGene 210521.xlsx", sheet = 1) %>% 
  #filter(!str_detect(Origine, "CHC")) %>% # Remove CHC patients due to the current lack of metadata; UPD: already have
  select(PATID = PatId, DATE = `Date prel.`)

# Combine patients from both waves
df <- 
  full_join(wave1, wave2) %>% 
  group_by(PATID) %>% 
  arrange(DATE, .by_group=T) %>% 
  distinct()

# Get time points for each patient assuming day of the first collection is zero
df <- df %>% 
  mutate(DAY = 
  {diff(DATE)} %>% 
    as.numeric() %>% 
    cumsum() %>% 
    c(0, .)
  ) %>% 
  select(PATID, DAY)

# Arrange by largest sampling date
N <- 50 # set up the threshold by day, i.e. to include samples with >= 3 time points earlier than this threshold
df <- df %>% 
  filter(DAY <= N) %>% 
  mutate(n = n()) %>% 
  filter(n >= 3) %>% 
  mutate(MAX_DAY = max(DAY)) %>% 
  arrange(desc(MAX_DAY)) %>% 
  select(-MAX_DAY) %>% 
  ungroup %>% 
  mutate(PATID = factor(PATID, levels = unique(PATID)))

# genotyped <- read_lines("data/raw/15042021_qc_step_3_pca.samples.to_keep.txt") %>% str_remove_all("COV.|.DUP.1") %>% str_subset("^B00[0-9]*")
genotyped <- read_delim("data/raw/GT/liege/cov_liege_15042021.fam", delim = " ", col_names = FALSE) %>% pull(2) %>% str_remove_all("COV.|.DUP.1")

id2barcode <- read_xlsx("data/raw/BelCovid-Uliege cohort_1st waive.xlsx", sheet = 5, skip = 9, col_names = FALSE) %>% select(PATID = 1, BARCODE = 3) %>% drop_na()
id2barcode$WAVE <- 0
id2barcode$GENOTYPED <- 0
for (i in 1:nrow(id2barcode)) {
  if (id2barcode$PATID[i] %in% wave1$PATID) {id2barcode$WAVE[i] <- 1} # so if patient is in both files, it will be assigned to wave1
  else if (id2barcode$PATID[i] %in% wave2$PATID) {id2barcode$WAVE[i] <- 2}
  
  if (id2barcode$BARCODE[i] %in% genotyped) {id2barcode$GENOTYPED[i] <- 1}
}

# check for how many patietns from wave 1 and/or 2 you do not have gt-data
wave1_nogt <- unique(wave1$PATID)[!unique(wave1$PATID) %in% {id2barcode %>% filter(GENOTYPED == 1) %>% pull(PATID)}]
wave2_nogt <- unique(wave2$PATID)[!unique(wave2$PATID) %in% {id2barcode %>% filter(GENOTYPED == 1) %>% pull(PATID)}] %>% .[!. %in% wave1_nogt]
length(unique(c(wave1_nogt, wave2_nogt)))
# So for 7 people from so-called wave1 and for 36 people from so-called wave2 there is no GT data up to now (2020-06-06), in total for 43.
# Ask Souad about tables of samples which will be genotyped in Finland, take those from *nogt variables who are not in the list - those are needed to be genotyped in addition to PDIAG and vaccinated.
