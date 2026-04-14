"""Generate icon.ico for the Shiny app launcher.

Draws a stylized volcano-plot motif (dark-blue square, dashed FC/FDR threshold
lines, red/blue/grey dots) at 256x256, then packs 16/32/48/64/128/256 into a
multi-resolution .ico.

Run from the app root:
    python3 scripts/make_icon.py
"""

from pathlib import Path
from PIL import Image, ImageDraw

BG = (8, 87, 120)          # #085778 — matches app palette
AXIS = (255, 255, 255, 160)
UP = (231, 76, 60)         # red
DOWN = (52, 152, 219)      # blue
NS = (204, 204, 204)       # grey

OUT = Path(__file__).resolve().parent.parent / "icon.ico"
SIZES = [(256, 256), (128, 128), (64, 64), (48, 48), (32, 32), (16, 16)]


def draw(size: int) -> Image.Image:
    img = Image.new("RGBA", (size, size), (0, 0, 0, 0))
    d = ImageDraw.Draw(img, "RGBA")

    radius = size // 6
    d.rounded_rectangle([0, 0, size - 1, size - 1], radius=radius, fill=BG)

    # Dashed threshold lines (FC verticals, FDR horizontal) — skip on tiny sizes
    if size >= 32:
        dash = max(2, size // 40)
        gap = dash
        # vertical lines at ~35% and ~65%
        for x_frac in (0.35, 0.65):
            x = int(size * x_frac)
            y = int(size * 0.15)
            while y < size * 0.9:
                d.line([(x, y), (x, min(y + dash, int(size * 0.9)))],
                       fill=AXIS, width=max(1, size // 128))
                y += dash + gap
        # horizontal line at ~40% from top
        y = int(size * 0.4)
        x = int(size * 0.1)
        while x < size * 0.9:
            d.line([(x, y), (min(x + dash, int(size * 0.9)), y)],
                   fill=AXIS, width=max(1, size // 128))
            x += dash + gap

    # Volcano dots — positions as fractions of size
    # (x_frac, y_frac, radius_frac, color)
    dots = [
        (0.22, 0.28, 0.06, DOWN),   # big down-regulated (top-left)
        (0.18, 0.45, 0.04, DOWN),
        (0.28, 0.38, 0.035, DOWN),
        (0.78, 0.25, 0.065, UP),    # big up-regulated (top-right)
        (0.72, 0.42, 0.045, UP),
        (0.82, 0.5, 0.035, UP),
        (0.5, 0.72, 0.05, NS),      # not-significant (bottom)
        (0.4, 0.78, 0.035, NS),
        (0.6, 0.78, 0.035, NS),
        (0.5, 0.88, 0.03, NS),
    ]

    # At small sizes, drop the tiny dots so the icon stays readable
    if size < 32:
        dots = [d_ for d_ in dots if d_[2] >= 0.05]
    elif size < 48:
        dots = [d_ for d_ in dots if d_[2] >= 0.04]

    for xf, yf, rf, color in dots:
        x = int(size * xf)
        y = int(size * yf)
        r = max(1, int(size * rf))
        d.ellipse([x - r, y - r, x + r, y + r], fill=color)

    return img


def main() -> None:
    base = draw(256)
    frames = [base.resize(sz, Image.LANCZOS) for sz in SIZES]
    frames[0].save(OUT, format="ICO", sizes=SIZES, append_images=frames[1:])
    print(f"wrote {OUT} ({OUT.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
