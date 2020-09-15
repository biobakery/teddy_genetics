#!/bin/sh

if [[ $# -eq 0 ]] ; then
      echo "Usage:"
      echo "    subset.sh [-h] <help> [-i] <input> [-p] <prefix> [-n] <size_chunks> [-t] <threads>"
    exit 0
fi

input=""
prefix=""

while getopts ":t:n:p:i:h" opt; do
  case ${opt} in
    h )
      echo "Usage:"
      echo "    subset.sh [-h] <help> [-i] <input> [-p] <prefix> [-n] <size_chunks> [-t] <threads>"
      echo "        -h display help"
      echo "        -i input file (.tsv) to be split"
      echo "        -p prefix for output files"
      echo "        -n size of chunks"
      echo "        -t threads"
      exit 0
      ;;
    i )
      i=${OPTARG}
      ;;
    p )
      p=${OPTARG}
      ;;
    n )
      n=${OPTARG}
      ;;
    t  )
      t=${OPTARG}
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      ;;
    : )
      echo "Invalid option: $OPTARG requires an argument" 1>&2
      ;;
  esac
done
shift $((OPTIND -1))

#

cat "${i}" | parallel --header : --pipe -N"${n}" 'cat > batch_{#}.txt'

ls batch_*.txt | awk '{print "Rscript transpose.R -i "$0" -o BATCH_"$0""}' > run_transpose.sh

cat run_transpose.sh | parallel -j "${t}"

for file in BATCH_*; do mv "${file}" "${file/BATCH_/${p}.}"; done

rm batch_*.txt

#
