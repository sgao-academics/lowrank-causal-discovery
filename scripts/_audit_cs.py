# -*- coding: utf-8 -*-
"""Audit cs_main.tex against the Current Science Research Article brief."""
import os, re

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = open(os.path.join(HERE, "cs_main.tex"), encoding="utf-8").read()

# ---------------------------------------------------------------- word count
def strip_tex(s):
    s = re.sub(r"(?m)%.*$", "", s)                       # comments
    s = re.sub(r"\\begin\{(table|figure)\}.*?\\end\{\1\}", " ", s, flags=re.S)
    s = re.sub(r"\\begin\{thebibliography\}.*?\\end\{thebibliography\}", " ", s, flags=re.S)
    s = re.sub(r"\\caption\{", " ", s)
    s = re.sub(r"\\includegraphics[^}]*\}", " ", s)
    s = re.sub(r"\\cite\{[^}]*\}", " ", s)               # citations are not words
    s = re.sub(r"\\label\{[^}]*\}", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?(\[[^\]]*\])?", " ", s)    # commands
    s = s.replace("{", " ").replace("}", " ")
    s = re.sub(r"[~\\]", " ", s)
    return s

def wc(text):
    return len([w for w in re.split(r"\s+", text.strip()) if re.search(r"[A-Za-z0-9]", w)])

body = SRC.split(r"\section*{Introduction}")[1]
body = body.split(r"\section*{Acknowledgements}")[0]
abstract = SRC.split(r"\begin{abstract}")[1].split(r"\end{abstract}")[0]
intro_to_concl = SRC.split(r"\section*{Introduction}")[1].split(r"\section*{Acknowledgements}")[0]

print("=" * 66)
print("WORD COUNT")
print("  abstract               : %4d" % wc(strip_tex(abstract)))
print("  body (Intro..Conclusion): %4d   <- Current Science limit 4000" % wc(strip_tex(body)))
print("  abstract + body        : %4d" % (wc(strip_tex(abstract)) + wc(strip_tex(body))))

# ------------------------------------------------------- display-item order
seq = [(m.start(), m.group(1)) for m in re.finditer(r"\\ref\{([^}]*)\}", SRC)]
print("\nDISPLAY-ITEM FIRST-REFERENCE ORDER (as they appear in the text)")
seen = []
for _, k in seq:
    if k not in seen:
        seen.append(k)
for i, k in enumerate(seen, 1):
    print("   %d. %s" % (i, k))

# labels declared, and where
decl = [(m.group(1), m.start()) for m in re.finditer(r"\\label\{([^}]*)\}", SRC)]
order_ok = True
pos = {k: i for i, k in enumerate(seen)}
for k, at in decl:
    if k not in pos:
        print("   !! label never referenced:", k); order_ok = False
print("   every float referenced:", order_ok)

# ------------------------------------------------------------- figure widths
print("\nFLOAT WIDTHS")
for m in re.finditer(r"\\includegraphics\[[^\]]*\{([^}]*)\}", SRC):
    print("   ", m.group(1))

# --------------------------------------------------------------- page count
log = open(os.path.join(HERE, "cs_main.log"), encoding="utf-8", errors="ignore").read()
m = re.search(r"Output written on cs_main\.pdf \((\d+) pages?, (\d+) bytes\)", log)
print("\nPDF: %s pages, %s bytes" % (m.group(1), m.group(2)) if m else "\nPDF: page line not found")
undef = re.findall(r"Citation `([^']*)' .*undefined", log) + \
        re.findall(r"Reference `([^']*)' .*undefined", log)
print("undefined citations/refs:", undef if undef else "none")
print("=" * 66)
