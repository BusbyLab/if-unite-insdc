# Fix the fungal taxonomy in the UNITE + INSDC release ####

# Load packages ####
library(Biostrings)
library(tidyr)
library(stringr)
library(readr)
library(dplyr)

# Define input/output directories ####
in.path <- '01-prep'
out <- '02-rename'
logs <- file.path(out,  'logs')
unlink(out, recursive = T)
dir.create(logs, recursive = T)

# Identify ranks of interest ####
ranks <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

# Read in the FASTA file
full <- readDNAStringSet(file.path(in.path, 'fun.fa.gz'), )

# Read in the unique fungal taxa ####
fun <- full |> names() |>
    unique() |> data.frame(tax = _) |>
    separate(tax, ranks, ';')
fun$taxon <- fun$species |>
    str_remove("_sp$") |> 
    str_replace_all('_', ' ')

# Read in and process the Index Fungorum reference ####
fg <- read_csv(file.path('data', 'IFexportGN-2023-03-07.csv'), skip = 1,
               col_names = c('record', 'taxon', 'author', 'year', 'current', ranks[1:5])) |>
    filter(is.na(current) == F,
           is.na(kingdom) == F) |> 
    mutate(parts = str_count(taxon, "[[:graph:]]+"),
           genus = taxon |> str_extract("^[[:graph:]]+"),
           species = if_else(parts >= 2,
                             taxon |> str_replace_all(' ', '_'), paste0(taxon, '_sp'))
           )
fg$record <- fg$record |> as.character()
fg$current <- fg$current |> as.character()

# Join the unique fungal taxa with the Index Fungorum reference ####
join <- inner_join(fun, fg, by = 'taxon') |>
    select(ends_with('.x'), record = current)

# Rejoin the portion that matched, setting the current record as the new record to match on ####
rejoin <- inner_join(join, fg, by = 'record') |> unique() |> 
    group_by(species.x) |>
    mutate(occurrences = n()) |>
    ungroup() |> 
    mutate(keep.stricto = if_else(occurrences == 1, T, F),
           keep.lato = if_else(occurrences == 1 |
                                   occurrences > 1 & record == current &
                                   kingdom.x == kingdom &
                                   phylum.x == phylum &
                                   class.x == class &
                                   order.x == order &
                                   family.x == family &
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
           family = family |> str_replace('_', '_family_incertae_sedis'))

# Create a conversion key ####
key <- fill |> mutate(header.in = paste(kingdom.x, phylum.x, class.x, order.x, family.x, genus.x, species.x, sep = ';'),
                      header.out = paste(kingdom, phylum, class, order, family, genus, species, sep = ';')) |> 
    select(header.in, header.out) |> 
    arrange(header.in) |> unique()

replicas <- names(full)[names(full) %in% key$header.in] |>
    sort() |> data.frame(header.in = _) |> count(header.in) %>% 
    lapply(rep, .$n) |> as.data.frame() |>
    left_join(key, by = 'header.in')

# Only retain recognized sequences in a specified order ####
full <- full[replicas$header.in]

# Update the names ####
names(full) <- replicas$header.out
full |> writeXStringSet(file.path(out, 'if-unite-insdc.fa'), width = 11000)

# Dereplicate the sequences ####
system(paste('vsearch --fasta_width 0 --derep_id', file.path(out, 'if-unite-insdc.fa'),
             '--log', file.path(logs, 'derep.txt'),
             '--output - | gzip -5 >', file.path(out, 'if-unite-insdc.fa.gz')))

# Remove the uncompressed file ####
unlink(file.path(out, 'if-unite-insdc.fa'))
