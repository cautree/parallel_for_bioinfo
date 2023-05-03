

process get_minimap_res {
errorStrategy 'ignore'


input:
tuple val(pair_id), val(name_fa_1),path(fa_1),  val(name_fa_2), path(fa_2)

output:
path ("*")

"""

minimap2 -x ava-ont ${fa_1} ${fa_2} | cut -f 1-12 | awk 'BEGIN { OFS = "," } ;    {\$1="$name_fa_2";  \$6="$name_fa_1";print}'  > ${pair_id}.overlaps.txt

"""

}

// the use of OFS in awk to save file in csv file
// the use of \ inside awk