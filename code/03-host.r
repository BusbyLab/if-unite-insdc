# Develop a host-specific and taxonomy-corrected reference from the UNITE + INSDC release ####

# Determine whether enough arguments were supplied ####
arguments <- commandArgs(T)

if(length(arguments) != 2){
    stop('Error: 1st argument = # of threads for parallel processing, 2nd argument = host taxon')
} else {
    writeLines('Two arguments provided')
}

# Identify each argument ####
threads <- arguments[[1]] |> as.integer()
plant <- arguments[[2]] |> as.character()

# Check arguments for validity ####
if(is.na(threads) == T){
    stop('Error: Only integer arguments accepted for first argument')
}
if(is.na(plant) == T){
    stop('Error: Only character arguments accepted for second argument')
}
if(threads < 1){
    stop('Error: At least one thread is needed')
} else {
    writeLines(c(paste(threads, 'threads requested'),
                 paste(plant, 'set as host')))
}

# Load packages ####
library(Biostrings)
library(tidyr)
library(stringr)
library(readr)
library(dplyr)

# Fix the host taxon argument, if needed ####
plant <- plant |> str_to_sentence()

# Define input/output directories ####
in.path <- '02-rename'
out <- '03-host'
unlink(out, recursive = T)
dir.create(out, recursive = T)

# Identify ranks of interest ####
ranks <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

# Read in the FASTA file ####
fun <- readDNAStringSet(file.path(in.path, 'if-unite-insdc.fa.gz'))
names(fun) <- names(fun) |> str_remove(" [[:digit:]]+$")

# Read in and process the Index Fungorum reference ####
fg <- read_csv(list.files('data', 'IFexportGN-*', full.names = T), skip = 1,
               col_names = c('record', 'taxon', 'author', 'year', 'current', ranks[1:5]),
               num_threads = threads) |>
    filter(is.na(current) == F,
           is.na(kingdom) == F,
           kingdom == 'Fungi') |> 
    mutate(parts = str_count(taxon, "[[:graph:]]+"),
           genus = taxon |> str_extract("^[[:graph:]]+"),
           species = if_else(parts >= 2,
                             taxon |> str_replace_all(' ', '_'), paste0(taxon, '_sp'))
    )
fg$record <- fg$record |> as.character()
fg$current <- fg$current |> as.character()

# Read in the United States National Fungus Collections Fungus-Host Dataset ####
unfiltered <- read_csv(list.files('data', 'Fungus-Host-Data_*', full.names = T), skip = 1,
                       col_names = c('fhcounter', 'taxon', 'host', 'fhlrcounter',
                                     'litnum', 'state', 'country', 'notes', 'citation'),
                       num_threads = threads)

# Check if the specified host is in the unfiltered dataset ####
check <- grepl(plant, unfiltered$host)
if(check[check == T] |> length() < 1){
    stop(paste(plant, 'not found'))
} else {
    writeLines(paste(plant, 'found'))
}

# Continue filtering and processing the dataset if the host is present ####
usda <- unfiltered |>
    filter(grepl(plant, host) == T) |> 
    mutate(taxon = taxon |> str_remove(" sp.$"),
           parts = str_count(taxon, "[[:graph:]]+"),
           genus = taxon |> str_extract("^[[:graph:]]+"),
           species = if_else(parts >= 2,
                             taxon |> str_replace_all(' ', '_'), paste0(taxon, '_sp'))
    )

# Join the host-associated fungal taxa with the Index Fungorum reference ####
join <- inner_join(usda, fg, by = 'taxon', relationship = 'many-to-many') |>
    select(ends_with('.x'),
           record = current,
           fhcounter, host, fhlrcounter, litnum, state, country, notes, citation)

# Rejoin the portion that matched, setting the current record as the new record to match on ####
rejoin <- inner_join(join, fg, by = 'record') |> unique() |> 
    group_by(species.x) |>
    mutate(occurrences = n()) |>
    ungroup() |> 
    mutate(keep.stricto = if_else(occurrences == 1, T, F),
           keep.lato = if_else(occurrences == 1 |
                                   occurrences > 1 & record == current &
                                   genus.x == genus &
                                   species.x == species, T, F)) |> 
    filter(keep.lato == T)

# Fix uncertain ranks ####
fill <- rejoin |> mutate(phylum = if_else(grepl('certae sedis', phylum) == T, paste0(kingdom, '_'), phylum),
                         class = if_else(grepl('certae sedis', class) == T, paste0(phylum, '_'), class),
                         order = if_else(grepl('certae sedis', order) == T, paste0(class, '_'), order),
                         family = if_else(grepl('certae sedis', family) == T, paste0(order, '_'), family)) |>
    mutate(across(c(phylum, class, order, family), \(x) str_replace(x, "_+", '_')),
           phylum = phylum |> str_replace('_', '_phylum_incertae_sedis'),
           class = class |> str_replace('_', '_class_incertae_sedis'),
           order = order |> str_replace('_', '_order_incertae_sedis'),
           family = family |> str_replace('_', '_family_incertae_sedis'),
           header = paste(kingdom, phylum, class, order, family, genus, species, sep = ';'))

# Get a vector of unique headers ####
headers <- fill$header |> unique()

# Write out the subset database ####
database <- fill |>
    select(-ends_with('.x'), -starts_with('keep.'),
           -parts, -occurrences, -header)
database |> write.csv(file.path(out, 'host.csv'), row.names = F)

# Keep sequences for host-associated taxa ####
keep <- fun[names(fun) %in% headers]
keep <- keep[sort(names(keep))]
keep |> writeXStringSet(file.path(out, 'host.fa.gz'),
                        width = 20001,
                        compress = T)
