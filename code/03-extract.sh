# Extract ITS and 5.8S subregions from full-ITS reference sequences ####

# Check the validity of the parallelization argument ####
if [[ $# -lt 1 ]]; then
  echo 'Error: Please specify the number of threads to launch'
  exit 1
fi

if [[ $# -gt 1 ]]; then
  echo 'Error: Too many arguments have been provided'
  exit 1
fi

if [[ $1 =~ [^0-9]+ ]]; then
  echo 'Error: Only integer arguments are accepted'
  exit 1
fi

if [[ $1 -lt 1 ]]; then
  echo 'Error: At least one thread is needed'
  exit 1
fi

# Prepare directories for input, output and log files ####
in='02-rename'
host='03-host'
out='04-extract'
logs="${out}/logs"
its1="${out}/its1"
five8s="${out}/5.8s"
its2="${out}/its2"
rm -r $out
mkdir -p scratch $logs $its1 $five8s $its2
touch $out/README.md

# Uncompress the inputs ####
pigz -p 4 -cd $in/if-unite-insdc.fa.gz > scratch/if-unite-insdc.fa
pigz -p 4 -cd data/sh_general_release_dynamic_all_* > scratch/sh_general_release_dynamic_all.fa

# Extract the ITS1 region from each sequence in the modified UNITE fungal release ####
ITSx \
-i scratch/if-unite-insdc.fa \
--preserve T \
--save_regions 'ITS1' \
--complement F -t "fungi" \
--graphical F \
--cpu $1 \
-o $its1/fun

# Extract the 5.8S region from each sequence in the modified UNITE fungal release ####
ITSx \
-i scratch/if-unite-insdc.fa \
--preserve T \
--save_regions '5.8S' \
--complement F \
-t "fungi" \
--only_full T \
--graphical F \
--cpu $1 \
-o $five8s/fun

# Extract the 5.8S region from each sequence in the UNITE eukaryote release (no singletons) ####
ITSx \
-i scratch/sh_general_release_dynamic_all.fa \
--preserve T \
--save_regions '5.8S' \
--complement F \
--only_full T \
--graphical F \
--cpu $1 \
-o $five8s/euk

# Extract the ITS2 region from each sequence in the modified UNITE fungal release ####
ITSx \
-i scratch/if-unite-insdc.fa \
--preserve T \
--save_regions 'ITS2' \
--complement F \
-t "fungi" \
--graphical F \
--cpu $1 \
-o $its2/fun

# Compress files in each output directory ####
pigz -p $1 $its1/*
pigz -p $1 $five8s/*
pigz -p $1 $its2/*

# Remove scratch output ####
rm -r scratch
