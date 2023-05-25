## separate files into different files, but awk on csv is tricky, better with tsv
awk -F"," '
  NR > 1 {
    if($NF == "Fall") {
      print $1, $NF > "fall.txt"
    } else if($NF == "Winter") {
      print $1, $NF > "winter.txt"
    } else if($NF == "Spring") {
      print $1, $NF > "spring.txt"
    } else {
      print $1, $NF > "other.txt"
    }
  }
' data/wgs-chinook-samples.csv


## use the output from awk in bash loops
WINTERS=$(awk -F"," '$NF == "Winter" {print $1}' data/wgs-chinook-samples.csv)

for i in $WINTERS; do
  echo "FASTQS are: $i.R1.fq.gz     $i.R2.fq.gz"
done
