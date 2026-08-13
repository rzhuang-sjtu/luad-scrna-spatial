"""Combine S10 panel PNGs into a vertical preview."""
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
import os

PANEL_DIR = Path("${WORK_ROOT}/luad_figures/fig_s10/panels")
OUT       = Path("${WORK_ROOT}/luad_figures/fig_s10/fig_s10_preview.png")

panels = ["S10A","S10B","S10C","S10D","S10E","S10F","S10G","S10H","S10I","S10J"]

# Layout: 2 cols × 5 rows for compact preview
TILE_W = 1200
TILE_H = 700
GAP = 12
LABEL_H = 28
cols = 2
rows = 5
canvas_w = cols*(TILE_W+GAP) + GAP
canvas_h = rows*(TILE_H+LABEL_H+GAP) + GAP
canvas = Image.new("RGB", (canvas_w, canvas_h), "white")

font = None
for fp in ["/usr/share/fonts/truetype/dejavu/DejaVu-Sans-Bold.ttf",
          "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
    if os.path.exists(fp):
        try:
            font = ImageFont.truetype(fp, 22); break
        except Exception: pass
if font is None: font = ImageFont.load_default()

draw = ImageDraw.Draw(canvas)
for i, name in enumerate(panels):
    src = PANEL_DIR / f"{name}.png"
    if not src.exists(): continue
    img = Image.open(src).convert("RGB")
    ratio = min(TILE_W/img.width, TILE_H/img.height)
    nw, nh = int(img.width*ratio), int(img.height*ratio)
    img = img.resize((nw, nh), Image.LANCZOS)
    r, c = divmod(i, cols)
    x0 = GAP + c*(TILE_W+GAP)
    y0 = GAP + r*(TILE_H+LABEL_H+GAP)
    px = x0 + (TILE_W - nw)//2
    py = y0 + LABEL_H + (TILE_H - nh)//2
    canvas.paste(img, (px, py))
    draw.text((x0+6, y0+2), name, fill="black", font=font)
    print(f"  placed {name}")

canvas.save(OUT, "PNG", optimize=True)
print(f"[done] preview -> {OUT}  ({canvas.size[0]}x{canvas.size[1]})")
