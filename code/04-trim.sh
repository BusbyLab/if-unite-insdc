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
in='03-extract'
its1='its1/fun.ITS1.fasta.gz'
euk='5.8s/euk.5_8S.fasta.gz'
fun='5.8s/fun.5_8S.fasta.gz'
its2='its2/fun.ITS2.fasta.gz'
out='04-trim'
rm -r $out
mkdir -p $out/its1 $out/5.8s $out/its2
touch $out/README.md

# Remove unique numbers from the header of each entry ####
pigz -p 4 -cd $in/$its1 | \
sed 's/ [0-9]+$//' | \
pigz -p $1 > $out/$its1

pigz -p 4 -cd $in/$euk | \
sed 's/ [0-9]+$//' | \
pigz -p $1 > $out/$euk

pigz -p 4 -cd $in/$fun | \
sed 's/ [0-9]+$//' | \
pigz -p $1 > $out/$fun

pigz -p 4 -cd $in/$its2 | \
sed 's/ [0-9]+$//' | \
pigz -p $1 > $out/$its2