# Extract ITS and 5.8S subregions from full-ITS modified reference sequences ####

# Accept an argument for the number of threads and check it for validity ####
threads <- commandArgs(T) |> as.integer()

if(is.na(threads) == T){
    stop('Error: Only integer arguments accepted for first argument')
}
if(length(threads) < 1){
    stop('Error: Please specify the number of threads to launch')
}
if(length(threads) > 1){
    stop('Error: Too many arguments have been provided')
}
if(threads < 1){
    stop('Error: At least one thread is needed')
} else {
    cat(threads, 'threads requested', '\n')
}

# Install packages from GitHub ####
remotes::install_github('brendanf/inferrnal@v0.99.5', dependencies = F)
remotes::install_github('brendanf/LSUx@v0.99.6', dependencies = F)

# Load packages ####
library(LSUx)

# Define input/output directories ####
in.path <- '02-rename'
out <- '03-alt'
unlink(out, recursive = T)
dir.create(out, recursive = T)

# Extract the fungal ITS1, 5.8S, and ITS2 subregions, specifying the truncated 32S covariance model ####
fun.extract <- lsux(file.path(in.path, 'if-unite-insdc.fa.gz'),
                    cm_32S = system.file(file.path('extdata', 'fungi_32S_LR5.cm'),
                                         package = 'LSUx'),
                    ITS1 = T,
                    cpu = threads)
fun.extract |> saveRDS(file.path(out, 'fun_extract.rds'))

# Extract the eukaryotic 5.8S and ITS2 subregions, specifying a eukaryotic 32S covariance model ####
# euk.extract <- lsux(list.files('data', 'sh_general_release_dynamic_all_*', full.names = T),
#                     cm_32S = file.path('data', 'RF02543.cm'),
#                     ITS1 = F,
#                     cpu = threads,
#                     mxsize = 16384)
# euk.extract |> saveRDS(file.path(out, 'euk_extract.rds'))