#!/bin/bash
#PBS -l select=1:ncpus=4:mem=32gb
#PBS -l walltime=02:00:00
#PBS -N indel_visualisation
#PBS -koed

echo $(date +%T)

module purge
module load tools/prod
module load SAMtools/1.22.1-GCC-14.2.0
module load BCFtools/1.22-GCC-14.2.0

eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
cd $WD

export REF_PATH="https://www.ebi.ac.uk/ena/cram/md5/%s"

# Outputs directory
OUT_DIR=${WD}/PDOs_outs/indel_analysis
mkdir -p ${OUT_DIR}

SAREK=/rds/general/project/spatialtranscriptomics/live/sarek_mutect
VCF1072=${SAREK}/variant_calling/mutect2/PDO_1072_vs_NT_1072/PDO_1072_vs_NT_1072.mutect2.filtered.vcf.gz
VCF1090=${SAREK}/variant_calling/mutect2/PDO_1090_vs_NT_1090/PDO_1090_vs_NT_1090.mutect2.filtered.vcf.gz

echo "Extracting PASS indels..."
bcftools view -f PASS -v indels ${VCF1072} -Oz -o ${OUT_DIR}/PDO_1072_pass_indels.vcf.gz
bcftools index ${OUT_DIR}/PDO_1072_pass_indels.vcf.gz

bcftools view -f PASS -v indels ${VCF1090} -Oz -o ${OUT_DIR}/PDO_1090_pass_indels.vcf.gz
bcftools index ${OUT_DIR}/PDO_1090_pass_indels.vcf.gz

echo "Generating quality metric TSVs..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/TLOD\t%INFO/DP\t%INFO/GERMQ\t%INFO/NALOD\t%INFO/NLOD\t%INFO/MBQ\t%INFO/MPOS\t[%AD\t%AF\t]\n' ${OUT_DIR}/PDO_1072_pass_indels.vcf.gz > ${OUT_DIR}/PDO_1072_indels_metrics.tsv

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/TLOD\t%INFO/DP\t%INFO/GERMQ\t%INFO/NALOD\t%INFO/NLOD\t%INFO/MBQ\t%INFO/MPOS\t[%AD\t%AF\t]\n' ${OUT_DIR}/PDO_1090_pass_indels.vcf.gz > ${OUT_DIR}/PDO_1090_indels_metrics.tsv

echo "Running R analysis and visualization..."
Rscript ${WD}/analysis/cnv/Auto_indel_visualisation.R

echo "Done!"
echo $(date +%T)
