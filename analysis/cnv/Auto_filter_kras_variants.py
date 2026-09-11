import os
import sys
import gzip
import csv
import requests
import json

# List of PDOs to process
pdos = [
    "PDO_1070_vs_NT_1070",
    "PDO_1072_vs_NT_1072",
    "PDO_1090_vs_NT_1090",
    "PDO_1121_vs_NT_1121",
    "PDO_1141_vs_NT_1141",
    "PDO_1181_vs_NT_1181",
    "PDO_629",
    "PDO_727"
]

vcf_base_dir = "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/annotation/mutect2"
output_csv = "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/PDOs_KRAS_oncogenic_variants.csv"

# Load the local OncoKB database provided by the user
oncokb_db_path = "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/annotation/mutect2/KRAS_variants.csv"
oncokb_db = {}
if os.path.exists(oncokb_db_path):
    with open(oncokb_db_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            oncokb_db[row['Alteration']] = row
else:
    print(f"Warning: {oncokb_db_path} not found.")

def is_oncogenic_gof_local(protein_change):
    if not protein_change or protein_change == "":
        return False, "Unknown", "Unknown"
        
    # Remove 'p.' prefix if present
    hgvsp = protein_change.replace("p.", "")
    
    if hgvsp in oncokb_db:
        entry = oncokb_db[hgvsp]
        oncogenic = entry.get("Oncogenic Status", "")
        effect = entry.get("Mutation Effect", "")
        
        is_onco = oncogenic in ["Oncogenic", "Likely Oncogenic"]
        is_gof = effect in ["Gain-of-function", "Likely Gain-of-function"]
        return (is_onco and is_gof), oncogenic, effect
        
    # If not in the local DB, assume False
    return False, "Unknown/Neutral (Not in DB)", "Unknown/Neutral (Not in DB)"

results = []

for pdo in pdos:
    vcf_path = os.path.join(vcf_base_dir, pdo, f"{pdo}.mutect2.filtered_VEP.ann.vcf.gz")
    
    if not os.path.exists(vcf_path):
        print(f"WARNING: VCF not found for {pdo} at {vcf_path}")
        continue
        
    print(f"Processing {pdo}...")
    
    with gzip.open(vcf_path, 'rt') as f:
        csq_format = []
        for line in f:
            if line.startswith("##INFO=<ID=CSQ"):
                # Extract Format from description
                desc = line.split("Format: ")[1].strip().strip('">')
                csq_format = desc.split("|")
                continue
                
            if line.startswith("#"):
                continue
                
            parts = line.strip().split("\t")
            info = parts[7]
            
            csq_str = ""
            for item in info.split(";"):
                if item.startswith("CSQ="):
                    csq_str = item[4:]
                    break
                    
            if not csq_str:
                continue
                
            # A variant can have multiple transcripts
            kras_found = False
            best_protein_change = ""
            best_consequence = ""
            
            for transcript_csq in csq_str.split(","):
                csq_vals = transcript_csq.split("|")
                if len(csq_vals) == len(csq_format):
                    csq_dict = dict(zip(csq_format, csq_vals))
                    symbol = csq_dict.get("SYMBOL", "")
                    
                    if symbol == "KRAS":
                        kras_found = True
                        hgvsp = csq_dict.get("HGVSp", "")
                        consequence = csq_dict.get("Consequence", "")
                        
                        aa_change = ""
                        # Try to get from HGVSp first
                        if ":" in hgvsp:
                            aa_change = hgvsp.split(":")[1]
                        
                        # Fallback to Amino_acids and Protein_position if HGVSp is empty
                        if not aa_change and csq_dict.get("Amino_acids") and csq_dict.get("Protein_position"):
                            amino_acids = csq_dict.get("Amino_acids").split("/")
                            pos = csq_dict.get("Protein_position").split("/")[0]
                            if len(amino_acids) == 2:
                                ref_aa = amino_acids[0]
                                alt_aa = amino_acids[1]
                                aa_change = f"{ref_aa}{pos}{alt_aa}"
                        
                        # Prioritize missense variants with a protein change
                        if "missense_variant" in consequence and aa_change:
                            best_protein_change = aa_change
                            best_consequence = consequence
                            break
                                
            if kras_found and best_protein_change:
                is_activating, oncogenicity, effect = is_oncogenic_gof_local(best_protein_change)
                
                if is_activating:
                    # Get AF if available
                    af = ""
                    format_fields = parts[8].split(":")
                    if "AF" in format_fields:
                        af_idx = format_fields.index("AF")
                        tumor_sample_vals = parts[9].split(":") 
                        if len(tumor_sample_vals) > af_idx:
                            af = tumor_sample_vals[af_idx]
                            
                    results.append({
                        "Sample": pdo,
                        "Chrom": parts[0],
                        "Pos": parts[1],
                        "Ref": parts[3],
                        "Alt": parts[4],
                        "Protein_Change": best_protein_change,
                        "Consequence": best_consequence,
                        "Tumor_AF": af,
                        "Oncogenicity": oncogenicity,
                        "Mutation_Effect": effect
                    })

# Write outputs
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ["Sample", "Chrom", "Pos", "Ref", "Alt", "Protein_Change", "Consequence", "Tumor_AF", "Oncogenicity", "Mutation_Effect"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"Finished processing. Results written to {output_csv}")



import os
import csv
import glob

# List of PDOs
pdos = [
    "PDO_1070_vs_NT_1070",
    "PDO_1072_vs_NT_1072",
    "PDO_1090_vs_NT_1090",
    "PDO_1121_vs_NT_1121",
    "PDO_1141_vs_NT_1141",
    "PDO_1181_vs_NT_1181",
    "PDO_629",
    "PDO_727"
]

base_dir = "/rds/general/project/spatialtranscriptomics/live/sarek_mutect/variant_calling/cnvkit"
output_csv = "/rds/general/project/tumourheterogeneity1/live/PDOs_Pipeline/PDOs_outs/PDOs_KRAS_CNV_status.csv"

# KRAS locus (hg38)
kras_chr = "chr12"
kras_start = 25204789
kras_end = 25250929

results = []

for pdo in pdos:
    # CNVkit output files can be `.call.cns` or `.somatic.call.cns`. We'll look for `.call.cns`
    cns_file = os.path.join(base_dir, pdo, f"{pdo.split('_vs_NT')[0]}.call.cns")
    
    # Sometimes the prefix might just be the PDO name or the full string
    if not os.path.exists(cns_file):
        cns_file = os.path.join(base_dir, pdo, f"{pdo}.call.cns")
        
    if not os.path.exists(cns_file):
        print(f"Warning: Could not find CNS file for {pdo}")
        continue
        
    with open(cns_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            chrom = row['chromosome']
            start = int(row['start'])
            end = int(row['end'])
            
            # Check for overlap with KRAS
            if chrom == kras_chr and start <= kras_end and end >= kras_start:
                log2 = row['log2']
                cn = row['cn']
                
                # Determine amplification status based on absolute copy number
                # Normal diploid is 2. Amplified could be > 2.
                cn_val = int(cn)
                if cn_val > 2:
                    status = "Amplified"
                elif cn_val < 2:
                    status = "Deleted"
                else:
                    status = "Neutral"
                    
                results.append({
                    "Sample": pdo,
                    "KRAS_Log2": log2,
                    "KRAS_CN": cn,
                    "Status": status
                })
                break # Found the segment for this sample

# Write output
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ["Sample", "KRAS_Log2", "KRAS_CN", "Status"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in results:
        writer.writerow(row)

print(f"Finished checking KRAS CNV. Results written to {output_csv}")