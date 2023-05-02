# Fix UNITE + INSD headers and update taxonomy ####

# Prepare directories for input, output and log files ####
in='data'
out='01-prep'
logs="${out}/logs"
rm -r $out
mkdir -p $logs
touch $out/README.md

echo 'Define the databases'
fun=$(find $in -name "UNITE_public*" | grep -v 'UNITE_public_all*')

echo 'Simplify fungal headers and remove duplicate sequences'
gzip -cd $fun | \
sed '/_gen_Incertae_sedis/{N;d;}' | sed -E "s/(^>)(.+\|)(.+)(\|.+)/\1\3/" | sed -E "s/[a-z]{1}__//g" | \
vsearch --fasta_width 0 --log $logs/derep-vsearch.txt \
--derep_id - --output - | gzip -5 > $out/fun.fa.gz