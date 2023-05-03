process get_blast_res {

errorStrategy 'ignore'

//publishDir path: "${params.outdir}/blastn_metric", mode: 'copy'

input:
tuple val(pair_id), val(name_fa_1),path(fa_1),  val(name_fa_2), path(fa_2)


output:
path ("*.csv")

"""

blastn -query ${fa_1} -subject ${fa_2}  > temp0.txt
cat temp0.txt | grep -E "Score|Identities|Strand" > temp.txt
cat temp.txt | grep "Score =" | head -3  | sed "s|Score = ||g"  |  sed "s|Expect = ||g" >  Score_Expect
cat temp.txt | grep Identities | head -3  | sed "s|Identities = ||g"  |  sed "s|Gaps = ||g" >  Identities_Gaps
cat temp.txt | grep Strand | head -3  | sed "s|Strand=||g"   >  Strand
cat temp0.txt | grep Length=    | head -n 1 |  sed "s|Length=||g"  | awk '1;1;1' > s1_length
cat temp0.txt | grep Length=    | tail -n 1 |  sed "s|Length=||g"  | awk '1;1;1' > s2_length

echo -e "${pair_id}\n${pair_id}\n${pair_id}" > pair_id
echo -e "${name_fa_1}\n${name_fa_1}\n${name_fa_1}" > name1
echo -e "${name_fa_2}\n${name_fa_2}\n${name_fa_2}" > name2


paste pair_id name1 name2  s1_length s2_length Score_Expect Identities_Gaps Strand | tr "\t" "," > ${pair_id}.csv

rm temp.txt Score_Expect Identities_Gaps Strand

"""


}

