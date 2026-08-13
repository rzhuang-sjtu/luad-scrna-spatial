"""Combine 16 panel PNGs into a 4x4 preview grid for Fig 7."""
from PIL import Image, ImageDraw, ImageFont
from pathlib import Path
import os

PANEL_DIR = Path("${WORK_ROOT}/luad_figures/fig7/panels")
OUT       = Path("${WORK_ROOT}/luad_figures/fig7/fig7_preview.png")

panels = ["7A","7B","7C","7D","7E","7F","7G","7H",
          "7I","7J","7K","7L","7M","7N","7O","7P"]

# uniform tile size (resize each panel to fit within 700x500)
TILE_W = 700
TILE_H = 500
GAP = 12
LABEL_H = 22

cols = 4
rows = 4
canvas_w = cols*(TILE_W+GAP) + GAP
canvas_h = rows*(TILE_H+LABEL_H+GAP) + GAP

canvas = Image.new("RGB", (canvas_w, canvas_h), "white")

# try to find an Arial font, fallback default
font = None
for fp in ["/usr/share/fonts/truetype/dejavu/DejaVu-Sans-Bold.ttf",
          "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
          "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf"]:
    if os.path.exists(fp):
        try:
            font = ImageFont.truetype(fp, 16)
            break
        except Exception:
            pass
if font is None:
    font = ImageFont.load_default()

draw = ImageDraw.Draw(canvas)
for i, name in enumerate(panels):
    src = PANEL_DIR / f"{name}.png"
    if not src.exists():
        print(f"  [skip] missing {src}")
        continue
    img = Image.open(src).convert("RGB")
    # fit into TILE_W x TILE_H preserving aspect
    ratio = min(TILE_W/img.width, TILE_H/img.height)
    nw, nh = int(img.width*ratio), int(img.height*ratio)
    img = img.resize((nw, nh), Image.LANCZOS)
    r, c = divmod(i, cols)
    x0 = GAP + c*(TILE_W+GAP)
    y0 = GAP + r*(TILE_H+LABEL_H+GAP)
    # paste centered
    px = x0 + (TILE_W - nw)//2
    py = y0 + LABEL_H + (TILE_H - nh)//2
    canvas.paste(img, (px, py))
    draw.text((x0+4, y0+2), name, fill="black", font=font)
    print(f"  placed {name} at row {r}, col {c}")

canvas.save(OUT, "PNG", optimize=True)
print(f"[done] preview saved -> {OUT}  ({canvas.size[0]}x{canvas.size[1]})")
