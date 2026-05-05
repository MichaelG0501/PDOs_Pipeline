#!/bin/bash
#PBS -l select=1:ncpus=2:mem=24gb
#PBS -l walltime=24:00:00
#PBS -N Auto_ASGARD_Download
#PBS -koed
echo $(date +%T)
module purge
module load tools/dev
WD=/rds/general/project/tumourheterogeneity1/ephemeral/PDOs_Pipeline
REF_ROOT=/rds/general/project/spatialtranscriptomics/ephemeral/Auto_drug_reversal_refs/asgard_l1000
RAW_DIR=$REF_ROOT/raw
PLAIN_DIR=$REF_ROOT/plain
mkdir -p "$RAW_DIR" "$PLAIN_DIR"
download_one() {
  url=$1
  file=$2
  if [[ ! -s "$RAW_DIR/$file" ]]; then
    curl -L --fail --show-error --continue-at - --output "$RAW_DIR/$file" "$url"
  fi
  plain=${file%.gz}
  if [[ ! -s "$PLAIN_DIR/$plain" ]]; then
    gzip -dc "$RAW_DIR/$file" > "$PLAIN_DIR/$plain"
  fi
}
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz" "GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz" "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz" "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz" "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_cell_info.txt.gz" "GSE92742_Broad_LINCS_cell_info.txt.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz" "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz"
download_one "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_sig_info.txt.gz" "GSE92742_Broad_LINCS_sig_info.txt.gz"
cd $WD
printf "asgard_l1000_reference_root,%s\n" "$REF_ROOT" > PDOs_outs/Auto_drug_reversal/asgard_reference_download_manifest.csv
echo $(date +%T)
