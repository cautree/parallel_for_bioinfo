BAM=$1

BIN_VERSION="1.4.0"
BASE="${HOME}/plant2/deepvariant-run"
INPUT_DIR="${BASE}/input"
REF="Cannabis_sativa.fa"
OUTPUT_DIR="${BASE}/output"
DATA_DIR="${INPUT_DIR}/data"
bn=$(basename $BAM .md.bam)
#OUTPUT_VCF="230315-A10-6-F.output.vcf.gz"
#OUTPUT_GVCF="230315-A10-6-F.output.g.vcf.gz
OUTPUT_VCF=${bn}".output.vcf.gz"
OUTPUT_GVCF=${bn}".output.g.vcf.gz"
nproc=8

#mkdir -p "${OUTPUT_DIR}"
#mkdir -p "${INPUT_DIR}"
#mkdir -p "${DATA_DIR}"


sudo docker run \
    -v "${DATA_DIR}":"/input" \
    -v "${OUTPUT_DIR}:/output" \
    gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}"  \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="/input/${REF}" \
    --reads="/input/${BAM}" \
    --output_vcf=/output/${OUTPUT_VCF} \
    --output_gvcf=/output/${OUTPUT_GVCF} \
    --num_shards=$(nproc) \
    --intermediate_results_dir /output/intermediate_results_dir





