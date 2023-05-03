
params.bam_files = "s3://s3_bucket/bam/plasmid_*.md.bam"

params.ref = "./ref/Escherichia_coli.fa*"
params.fa= "./ref/Escherichia_coli.fa"
params.fa_index = "./ref/Escherichia_coli.fai"
params.wgs = true
params.dev = false
params.number_of_inputs = 2
params.analysis = "unmapped_to_Escherichia_coli"
params.run = "run1"
params.outfile ="run1_unmapped_mapping_to_Escherichia_coli"
params.outdir = "./nextflow_results"


if(params.dev) { 
   path_s3 = "company-dev/analysis"
} else { 
   path_s3 = "company-analysis"
}



process get_unmapped {

publishDir path: "s3://$path_s3/$params.run/$params.analysis/ecoli/unmapped_bam"

input:
tuple val(pair_id), path(bam_file)

output:
tuple val(pair_id), path("*unmapped.bam") 

"""
samtools view -b -f 4 $bam_file | samtools sort  -  -o ${pair_id}.unmapped.bam

"""

}




process bam_to_fastq {

publishDir path: "s3://$path_s3/$params.run/$params.analysis/ecoli/unmapped_fastq"
//publishDir path: "${params.outdir}/fastq", mode: "copy"


input:
tuple val(pair_id), path(bam_file)

output:
tuple val(pair_id), path("*_r1.fastq"), path("*_r2.fastq")


"""
bedtools bamtofastq -i $bam_file -fq ${pair_id}_r1.fastq -fq2 ${pair_id}_r2.fastq

"""


}



workflow{

ref_files = file( params.ref)
fa_file = file(params.fa)
fa_index_file = file(params.fa_index)

bam_files = Channel
     .fromPath( params.bam_files)
     .map{ item -> tuple(item.getBaseName(2), item) }
     .take( params.dev ? params.number_of_inputs : -1 )
     
unmapped = get_unmapped(bam_files)

fastq = bam_to_fastq(unmapped)

}


