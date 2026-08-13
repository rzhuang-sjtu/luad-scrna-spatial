"""Compose Figure S7 (panels a-f) into a single page PDF, no panel labels."""
from pathlib import Path
from pypdf import PdfReader, PdfWriter, PageObject, Transformation

MM = 72 / 25.4
W_MM, H_MM = 200, 460

src7 = Path("${WORK_ROOT}/luad_figures/fig_s7")
src8 = Path("${WORK_ROOT}/luad_figures/fig_s8")

# layout: list of (file, target_width_mm, target_height_mm, x_mm, y_mm)
# y measured from BOTTOM of page (PDF convention)
panels = []

def place(p, tw, th, x, y):
    panels.append((p, tw, th, x, y))

# Row 1: A | B  (UMAPs, 95 wide x 70 tall)
place(src7 / "figS7a_umap_tissue.pdf", 95, 70,    5, H_MM - 5  - 70)
place(src7 / "figS7b_umap_dataset.pdf", 95, 70, 105, H_MM - 5  - 70)
# Row 2: C | D  (ridges + confusion)
place(src7 / "figS7c_scanvi_ridges.pdf", 95, 75,   5, H_MM - 80 - 75)
place(src7 / "figS7d_confusion.pdf",     95, 65, 105, H_MM - 85 - 65)
# Row 3: E full width (heatmap, 190 x 111)
place(src8 / "figS7e_hallmark_heatmap.pdf", 190, 111, 5, H_MM - 165 - 111)
# Row 4: F full width (GO grid, 180 x 152)
place(src8 / "figS7f_go_enrichment.pdf",    190, 152, 5, H_MM - 281 - 152)

writer = PdfWriter()
canvas = PageObject.create_blank_page(width=W_MM * MM, height=H_MM * MM)
for path, tw, th, x_mm, y_mm in panels:
    if not path.exists():
        print(f"  MISS {path}")
        continue
    src = PdfReader(str(path))
    pg = src.pages[0]
    sw = float(pg.mediabox.width)
    sh = float(pg.mediabox.height)
    sx = (tw * MM) / sw
    sy = (th * MM) / sh
    t = Transformation().scale(sx=sx, sy=sy).translate(tx=x_mm * MM, ty=y_mm * MM)
    canvas.merge_transformed_page(pg, t)
    print(f"  placed {path.name} at ({x_mm:.0f}, {y_mm:.0f})  size {tw} x {th} mm")

writer.add_page(canvas)
out = Path("${WORK_ROOT}/luad_figures/fig_s7/FigureS7.pdf")
with open(out, "wb") as f:
    writer.write(f)
print(f"\nwrote {out}  canvas {W_MM} x {H_MM} mm")
