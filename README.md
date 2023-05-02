IF + UNITE + INSDC: improved fungal taxonomy for ITS reference sequences
================
Kyle A. Gervers
2023-05-02

## Overview

While the UNITE general FASTA releases represent high-quality,
dynamically clustered references for fungal taxon assignment, most
releases consistently ignore taxa that I care about. Also, I’ve noticed
that the taxonomic ranks applied are not always the most up-to-date when
compared with what Index Fungorum reports, and many of the included
sequences lack detailed taxonomic resolution.

Starting from the most recent (as of 2023-05-02) fungal UNITE+INSD
release, this repo does the following:

- Removes sequences not identified to genus
- De-replicates sequences with the same header taxonomy and sequence
- Updates the taxonomy of the remaining sequences according to Index
  Fungorum
- De-replicates sequences (again) with the same header taxonomy and
  sequence
- Formats sequences headers in a format compatible with the
  `assignTaxonomy()` function from `dada2`

Because UNITE does not include the authority in the applied taxonomy,
this approach ends up shedding more sequences as a conservative measure,
only including sequences with unambiguous taxonomy. Also, because data
dumps of Index Fungorum taxonomy (which are difficult to find anyways!)
appear to only allow indexing at the genus-level, a genus-level ID
requirement is placed on the UNITE+INSD release. This approach also
doesn’t apply the dynamic clustering that UNITE applies for it’s species
hypotheses and operates with the assumption that the IDs in the
UNITE+INSD release (which have apparently already undergone some QAQC by
UNITE) are reliable.

The updated fungal taxonomy reference is found in the file
`if-unite-insdc.fa.gz`, located in the `02-rename` folder.

This whole repo will (hopefully) become obsolete if:

- future releases include species hypotheses of taxa currently only
  found in UNITE+INSD releases
- fungal taxonomy of future releases matches Index Fungorum’s

Planned improvements include:

- inclusion of sequences at higher ranks (dependent on Index Fungorum
  data dumps)
- application of dynamic clustering thresholds

## Reproducibility

All packages were installed and managed with `conda`.

    conda 23.3.1
    name: /home/gerverska/projects/if-unite-insdc/env
    channels:
      - conda-forge
      - bioconda
      - defaults
    dependencies:
      - bioconductor-biostrings
      - r-base=4.2.2
      - r-dplyr
      - r-markdown
      - r-readr
      - r-rmarkdown
      - r-stringr
      - r-tidyr
      - vsearch
    prefix: /home/gerverska/projects/if-unite-insdc/env

Install the above bioinformatic environment from `config.yml` using the
script `00-build.sh`

    # Clone the repo (using the GitHub CLI tool) ####
    gh repo clone gerverska/if-unite-insdc

    # Run the build script ####
    bash code/00-build.sh
