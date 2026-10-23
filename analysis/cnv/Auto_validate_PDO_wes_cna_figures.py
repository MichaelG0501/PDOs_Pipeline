"""Validate and preview the six-PDO WES comparison PDFs on a PBS node.

Status: active terminal helper for Auto_PDO_wes_cna_mp_validation.R.
Methodology: analysis/methodology/cnv/Auto_PDO_wes_cna_mp_validation_methodology.md.
Input: output root argument and the exact pinned PyMuPDF library in scRef.
Outputs: live reports/Auto_pdf_validation.csv, Auto_*_preview.png.
Run: active gnmf Python analysis/cnv/Auto_validate_PDO_wes_cna_figures.py <output_root>.
No files outside the supplied output root are changed.
"""
import csv
import os
from pathlib import Path
import sys

if not os.environ.get("PBS_JOBID"):
    raise RuntimeError("Run on a PBS compute node")
sys.path.insert(0, "/rds/general/project/tumourheterogeneity1/live/scRef_Pipeline/ref_outs/TCGA/cna_gene_event_validation/intermediate/source_data/software/pymupdf_1.26.7")
import fitz

root = Path(sys.argv[1]).resolve()
expected = {"Auto_gene_cna_mp_rectangles.pdf": 4, "Auto_arm_cna_mp_rectangles.pdf": 2,
            "Auto_gene_cna_mp_boxplots.pdf": 10, "Auto_arm_cna_mp_boxplots.pdf": 10,
            "Auto_PDO_cna_mp_scatterplots.pdf": 20}
rows = []
for name, count in expected.items():
    path = root / "figures" / name
    if path.stat().st_size < 1000:
        raise RuntimeError(f"Empty PDF: {path}")
    with fitz.open(path) as doc:
        if len(doc) != count:
            raise RuntimeError(f"Wrong page count: {name}, {len(doc)} != {count}")
        for index, page in enumerate(doc):
            text = page.get_text()
            if len(text.strip()) < 80:
                raise RuntimeError(f"Missing text: {name} page {index + 1}")
            outside = []
            for block in page.get_text("dict")["blocks"]:
                for line in block.get("lines", []):
                    for span in line["spans"]:
                        box = fitz.Rect(span["bbox"])
                        if box.x0 < -1 or box.y0 < -1 or box.x1 > page.rect.width + 1 or box.y1 > page.rect.height + 1:
                            outside.append(span["text"])
            rows.append({"file": name, "page": index + 1, "text_characters": len(text),
                         "outside_text_spans": len(outside), "bytes": path.stat().st_size})
            if outside:
                raise RuntimeError(f"Clipped text: {name} page {index + 1}: {outside}")
            if index == 0:
                page.get_pixmap(matrix=fitz.Matrix(0.65, 0.65)).save(root / "reports" / f"{path.stem}_preview.png")
with (root / "reports" / "Auto_pdf_validation.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
print(f"PASS: {len(expected)} PDFs, {len(rows)} pages; text bounds and first-page previews verified")
