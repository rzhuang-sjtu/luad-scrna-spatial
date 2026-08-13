"""Combine 16 S12 lite PDFs into A4 multi-page composites (portrait + landscape)."""
from __future__ import annotations
from pathlib import Path
import math
from pypdf import PdfWriter, PdfReader, Transformation, PageObject
from pypdf.generic import RectangleObject

SRC_DIR = Path("${WORK_ROOT}/luad_figures/fig_s12/panels")
DST_DIR = Path("${WORK_ROOT}/luad_figures/fig_s12")
LETTERS = list("ABCDEFGHIJKLMNOP")
LITES = [SRC_DIR / f"S12{L}_lite.pdf" for L in LETTERS]

# A4 in points (72 dpi): 595.28 x 841.89
A4_W = 595.28
A4_H = 841.89
MARGIN = 14   # points (~5 mm)
GAP    = 6


def nup_pages(in_pdfs, page_w, page_h, ncol, nrow, margin=MARGIN, gap=GAP):
    """Return a list of PageObject, each an N-up composite of input PDF first-pages."""
    avail_w = page_w - 2 * margin - gap * (ncol - 1)
    avail_h = page_h - 2 * margin - gap * (nrow - 1)
    cell_w  = avail_w / ncol
    cell_h  = avail_h / nrow
    pages = []
    per_page = ncol * nrow
    for chunk_start in range(0, len(in_pdfs), per_page):
        chunk = in_pdfs[chunk_start:chunk_start + per_page]
        comp = PageObject.create_blank_page(width = page_w, height = page_h)
        for idx, src in enumerate(chunk):
            r = idx // ncol
            c = idx % ncol
            x0 = margin + c * (cell_w + gap)
            y0 = page_h - margin - (r + 1) * cell_h - r * gap
            src_page = PdfReader(str(src)).pages[0]
            sw = float(src_page.mediabox.width)
            sh = float(src_page.mediabox.height)
            scale = min(cell_w / sw, cell_h / sh)
            new_w = sw * scale
            new_h = sh * scale
            # centre within cell
            xoff = x0 + (cell_w - new_w) / 2
            yoff = y0 + (cell_h - new_h) / 2
            t = Transformation().scale(scale).translate(xoff, yoff)
            comp.merge_transformed_page(src_page, t)
        pages.append(comp)
    return pages


def build(out_path: Path, page_w, page_h, ncol, nrow):
    pages = nup_pages(LITES, page_w, page_h, ncol, nrow)
    w = PdfWriter()
    for p in pages:
        w.add_page(p)
    with out_path.open("wb") as f:
        w.write(f)
    print(f"  wrote {out_path.name}  ({len(pages)} pages, {ncol}x{nrow}/page)")


# Portrait A4: 4 cols x 4 rows -> 16/page = 1 page
build(DST_DIR / "S12_combined_A4_portrait.pdf",
      A4_W, A4_H, ncol=4, nrow=4)

# Landscape A4: 4 cols x 4 rows -> 16/page = 1 page
build(DST_DIR / "S12_combined_A4_landscape.pdf",
      A4_H, A4_W, ncol=4, nrow=4)

# Bonus: 2 cols x 8 rows portrait variant (taller cells, same A4)
build(DST_DIR / "S12_combined_A4_portrait_2x8.pdf",
      A4_W, A4_H, ncol=2, nrow=8)
# Bonus: 8 cols x 2 rows landscape variant (wider cells)
build(DST_DIR / "S12_combined_A4_landscape_8x2.pdf",
      A4_H, A4_W, ncol=8, nrow=2)

print("\nDONE.")
