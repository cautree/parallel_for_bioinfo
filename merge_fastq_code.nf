params.run = "20220616_lanes"
params.plate = "220616-FASTQ"


//220616-UDI-PlexA_H02_L004 from file 220616-UDI-PlexA_H02_L004_R2_001.fastq.gz become 220616-UDI-PlexA_H02
def getSampleId( prefix ) {
     n=prefix.size()
     prefix.substring(0, n-5)

}

Channel
     .fromFilePairs("s3://s3_bucket/" + params.run + "/" + params.plate + "/" + "*_R{1,2}_001.fastq.gz", flat: true)
     .map { prefix, file1, file2 -> tuple(getSampleId(prefix), file1, file2) }
     .groupTuple()
     .set { fastq_ch }


process merge_fastq {

publishDir path: 'fastq', pattern: '*.fastq.gz', mode: 'copy'

input:
     tuple val(prefix), path(reads1), path(reads2) from fastq_ch

output:
     tuple val(prefix), path('*_R1_001.fastq.gz'), path('*_R2_001.fastq.gz') 


"""
     cat ${reads1}   | gzip > ${prefix}._R1_001.fastq.gz
     cat ${reads2}   | gzip > ${prefix}._R2_001.fastq.gz
     
"""


}



