#!/usr/bin/env python3
"""Build the EukDetect 2.1 pipeline schematic as an editable SVG.

Palette and shape language are taken from the version 1 Figure S1 so the two
figures sit together: blue = {output_dir}/, crimson = {output_dir}/aln/,
gold = {output_dir}/filtering/.

    make_schematic.py out.svg
"""
import sys

BLUE = "#4C99D2"
PINK = "#CD1E65"
GOLD = "#FBCA07"
GREY = "#E9E9E9"
INK = "#1A1A1A"
FONT = "Helvetica, Arial, sans-serif"

W, H = 1180, 1120


def esc(s):
    return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def text(x, y, s, size=13, anchor="middle", weight="normal", fill=INK, style=""):
    return (f'<text x="{x}" y="{y}" font-family="{FONT}" font-size="{size}" '
            f'font-weight="{weight}" text-anchor="{anchor}" fill="{fill}" '
            f'{style}>{esc(s)}</text>')


def lines(x, y, rows, size=13, anchor="middle", weight="normal", fill=INK, lh=16):
    out = []
    for i, r in enumerate(rows):
        out.append(text(x, y + i * lh, r, size, anchor, weight, fill))
    return out


def doc(x, y, w, h, fill, rows, size=13, text_fill=INK):
    """A file: rectangle with the wavy bottom edge used in the v1 figure."""
    a = h * 0.12
    d = (f"M{x},{y} L{x+w},{y} L{x+w},{y+h-a} "
         f"C{x+w*0.75},{y+h-a*2.4} {x+w*0.60},{y+h+a*0.9} {x+w*0.5},{y+h-a*0.35} "
         f"C{x+w*0.38},{y+h-a*2.0} {x+w*0.22},{y+h+a*0.6} {x},{y+h-a*0.5} Z")
    g = [f'<path d="{d}" fill="{fill}" stroke="{INK}" stroke-width="1.6"/>']
    ty = y + h / 2 - (len(rows) - 1) * 8 + 4
    g += lines(x + w / 2, ty, rows, size, fill=text_fill)
    return g


def box(x, y, w, h, rows, size=13, fill="#FFFFFF", stroke=INK, rx=4, weight="normal",
        dash=None):
    sd = f' stroke-dasharray="{dash}"' if dash else ""
    g = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{rx}" fill="{fill}" '
         f'stroke="{stroke}" stroke-width="1.6"{sd}/>']
    ty = y + h / 2 - (len(rows) - 1) * 8 + 4
    g += lines(x + w / 2, ty, rows, size)
    return g


def step(x, y, w, rows, size=13):
    """An unboxed processing step, as in the v1 figure."""
    return lines(x + w / 2, y, rows, size)


def arrow(x1, y1, x2, y2, width=2.2):
    return (f'<line x1="{x1}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="{INK}" '
            f'stroke-width="{width}" marker-end="url(#ah)"/>')


def elbow(x1, y1, x2, y2, width=2.2):
    """Down, across, down."""
    ym = (y1 + y2) / 2
    return (f'<path d="M{x1},{y1} L{x1},{ym} L{x2},{ym} L{x2},{y2}" fill="none" '
            f'stroke="{INK}" stroke-width="{width}" marker-end="url(#ah)"/>')


def build():
    s = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'viewBox="0 0 {W} {H}">',
         '<defs><marker id="ah" viewBox="0 0 10 10" refX="9" refY="5" '
         'markerWidth="6" markerHeight="6" orient="auto-start-reverse">'
         f'<path d="M0,0 L10,5 L0,10 z" fill="{INK}"/></marker></defs>',
         f'<rect width="{W}" height="{H}" fill="#FFFFFF"/>']

    s.append(text(40, 42, "Figure S1. Schematic of the EukDetect 2 pipeline.",
                  16, anchor="start", weight="bold"))

    # ---- inputs -----------------------------------------------------------
    s += doc(52, 96, 150, 76, "#FFFFFF",
             ["Raw metagenomic", "fastq files"])
    s += [f'<path d="M232,96 h150 v76 a75,16 0 0 1 -150,0 z" fill="#FFFFFF" '
          f'stroke="{INK}" stroke-width="1.6"/>',
          f'<ellipse cx="307" cy="96" rx="75" ry="16" fill="#FFFFFF" '
          f'stroke="{INK}" stroke-width="1.6"/>']
    s += lines(307, 128, ["EukDetect marker", "database + ANI table"], 13)

    # ---- legend -----------------------------------------------------------
    s += box(700, 84, 420, 128, [], rx=3)
    s.append(text(910, 110, "File location key", 14, weight="bold"))
    for i, (lab, col) in enumerate([("{output_dir}/", BLUE),
                                    ("{output_dir}/aln/", PINK),
                                    ("{output_dir}/filtering/", GOLD)]):
        y = 134 + i * 26
        s.append(text(1010, y + 5, lab, 13, anchor="end"))
        s.append(f'<rect x="1026" y="{y - 9}" width="58" height="18" fill="{col}" '
                 f'stroke="{INK}" stroke-width="1.2"/>')

    s.append(arrow(127, 176, 127, 232))
    s.append(arrow(307, 176, 307, 232))

    # ---- outer container: mode all ----------------------------------------
    s += box(36, 236, 1108, 846, [], fill="#FFFFFF", rx=3)
    s.append(text(56, 272, "--mode all", 19, anchor="start"))

    # ---- mode aln ---------------------------------------------------------
    AX, AY, AW, AH = 56, 290, 290, 346
    s += box(AX, AY, AW, AH, [], fill=GREY, rx=3)
    s.append(text(AX + 14, AY + 32, "--mode aln", 17, anchor="start"))

    s += step(AX, AY + 62, AW, ["Align to marker database", "(bowtie2 --very-sensitive)"])
    s.append(arrow(AX + AW / 2, AY + 86, AX + AW / 2, AY + 108))
    s += step(AX, AY + 128, AW, ["Filter by alignment length", "(> 80% of read length)"])
    s.append(arrow(AX + AW / 2, AY + 152, AX + AW / 2, AY + 174))
    s += step(AX, AY + 194, AW, ["Quality filter", "(MAPQ ≥ 10)"])
    s.append(arrow(AX + AW / 2, AY + 218, AX + AW / 2, AY + 244))
    s += doc(AX + 70, AY + 248, 150, 72, PINK, ["Raw alignments", "(.bam)"],
             text_fill="#FFFFFF")

    # ---- mode analyze -----------------------------------------------------
    NX, NY, NW, NH = 386, 290, 738, 772
    s += box(NX, NY, NW, NH, [], fill=GREY, rx=3)
    s.append(text(NX + 14, NY + 32, "--mode analyze", 17, anchor="start"))

    s.append(arrow(AX + AW, AY + 284, NX, AY + 284, 3))

    cx = NX + 300

    s += step(cx - 150, NY + 62, 300, ["Low-complexity filtering",
                                       "(4-mer complexity < 0.5)"])
    s.append(arrow(cx, NY + 86, cx, NY + 104))
    s += step(cx - 150, NY + 124, 300, ["Duplicate read removal",
                                        "(alignment position)"])
    s.append(arrow(cx, NY + 148, cx, NY + 166))
    s += doc(cx - 75, NY + 170, 150, 64, GOLD, ["Filtered alignments", "(.bam)"])

    lx, rx_ = cx - 150, cx + 150
    s.append(elbow(cx, NY + 240, lx, NY + 272))
    s.append(elbow(cx, NY + 240, rx_, NY + 272))

    s += step(lx - 85, NY + 292, 170, ["Count reads, coverage,", "percent identity"], 12)
    s += step(rx_ - 85, NY + 292, 170, ["Per-read sequence", "identity"], 12)
    s.append(arrow(lx, NY + 316, lx, NY + 334))
    s.append(arrow(rx_, NY + 316, rx_, NY + 334))

    s += doc(lx - 75, NY + 338, 150, 60, GOLD, ["Reads per marker", "(.tsv)"], 12)
    s += doc(rx_ - 75, NY + 338, 150, 60, GOLD, ["Read identities", "(.tsv)"], 12)

    DX, DY, DW, DH = cx - 215, NY + 432, 430, 148
    s.append(elbow(lx, NY + 402, cx - 90, DY - 4))
    s.append(elbow(rx_, NY + 402, cx + 90, DY - 4))

    s += box(DX, DY, DW, DH, [], fill="#FFFFFF", rx=3, dash="6 4")
    s.append(text(cx, DY + 24, "Species disambiguation", 14, weight="bold"))
    s += lines(cx, DY + 48,
               ["Taxa that compete for reads are grouped by",
                "ANI \u2265 80 OR a shared NCBI genus"], 12)
    s += lines(cx, DY + 90,
               ["Within a group, a Mann-Whitney test on read identity",
                "(\u03b1 = 0.05) and an abundance-scaled median gap decide",
                "whether a taxon keeps its reads or hands them over"], 12)

    s.append(f'<path d="M{DX+DW},{DY+DH/2} L{DX+DW+54},{DY+DH/2}" fill="none" '
             f'stroke="{INK}" stroke-width="2.2" marker-end="url(#ah)"/>')
    s += doc(DX + DW + 58, DY + DH / 2 - 32, 150, 64, GOLD,
             ["Reassignment report", "(.tsv)"], 12)

    s.append(arrow(cx, DY + DH, cx, DY + DH + 28))
    s += step(cx - 150, DY + DH + 48, 300, ["Filter taxa with \u2265 4 reads",
                                            "on \u2265 2 marker genes"])

    bar = DY + DH + 86
    oy = DY + DH + 100
    s.append(f'<path d="M{cx},{DY+DH+72} L{cx},{bar}" fill="none" '
             f'stroke="{INK}" stroke-width="2.2"/>')
    s.append(f'<path d="M{cx-215},{bar} L{cx+215},{bar}" fill="none" '
             f'stroke="{INK}" stroke-width="2.2"/>')
    for dx in (-215, 0, 215):
        s.append(arrow(cx + dx, bar, cx + dx, oy - 4))

    s += doc(cx - 290, oy, 150, 68, GOLD,
             ["All hits,", "before filtering", "(.tsv)"], 12)
    s += doc(cx - 75, oy, 150, 68, BLUE,
             ["Filtered hits", "per species", "(.tsv)"], 12, text_fill="#FFFFFF")
    s += doc(cx + 140, oy, 150, 68, BLUE,
             ["RPKS and relative", "abundance", "(.tsv)"], 12, text_fill="#FFFFFF")

    s.append(text(W / 2, H - 16,
                  "Tabular output is tab-separated with a .tsv extension; "
                  "missing values are written as NA.", 12, fill="#555555"))

    s.append("</svg>")
    return "\n".join(s)


if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "eukdetect2_schematic.svg"
    with open(out, "w") as f:
        f.write(build())
    print("wrote", out)
