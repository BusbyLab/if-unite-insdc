# Fix UNITE + INSD headers and update taxonomy ####

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
in='data'
out='01-prep'
logs="${out}/logs"
rm -r $out
mkdir -p $logs
touch $out/README.md

# Define the databases ####
fun=$(find $in -name "UNITE_public*" | grep -v 'UNITE_public_all*')

# Simplify fungal headers and remove duplicate sequences ####
pigz -p 4 -cd $fun | \
sed '/_gen_Incertae_sedis/{N;d;}' | sed -E "s/(^>)(.+\|)(.+)(\|.+)/\1\3/" | sed -E "s/[a-z]{1}__//g" | \
vsearch \
--fasta_width 0 \
--log $logs/derep-vsearch.txt \
--derep_id - \
--output - | pigz -p $1 > $out/fun.fa.gz