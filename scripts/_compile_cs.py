# -*- coding: utf-8 -*-
"""Compile cs_main.tex / cs_supp.tex and verify the PDF is newer than the call.

pdflatex is invoked through subprocess from the manuscript directory so that
relative \\includegraphics paths resolve, and so that no PowerShell quoting or
`cd /d` behaviour is involved.
"""
import os, subprocess, sys, time, shutil

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def find_pdflatex():
    p = shutil.which("pdflatex")
    if p:
        return p
    for c in (r"pdflatex.exe",
              r"C:/Program Files/MiKTeX/miktex/bin/x64/pdflatex.exe",
              r"C:/texlive/2024/bin/windows/pdflatex.exe",
              r"C:/texlive/2025/bin/windows/pdflatex.exe"):
        if os.path.exists(c):
            return c
    raise SystemExit("pdflatex not found")


def compile_one(pdftex, stem, passes=2):
    t0 = time.time()
    for i in range(passes):
        r = subprocess.run([pdftex, "-interaction=nonstopmode",
                            "-halt-on-error", stem + ".tex"],
                           cwd=HERE, capture_output=True)
        logf = os.path.join(HERE, stem + ".log")
        log = open(logf, encoding="utf-8", errors="ignore").read() if os.path.exists(logf) else ""
        errs = [l for l in log.splitlines() if l.startswith("!")]
        if r.returncode != 0 or errs:
            print("  pass %d FAILED (rc=%d)" % (i + 1, r.returncode))
            for l in errs[:12]:
                print("   ", l)
            tail = "\n".join(log.splitlines()[-25:])
            print(tail)
            return False
    pdf = os.path.join(HERE, stem + ".pdf")
    if not os.path.exists(pdf):
        print("  no PDF produced"); return False
    if os.path.getmtime(pdf) < t0 - 1:
        print("  STALE PDF: mtime older than this run -> treat as failure"); return False
    log = open(os.path.join(HERE, stem + ".log"), encoding="utf-8", errors="ignore").read()
    warns = [l for l in log.splitlines() if "Warning" in l]
    print("  OK  %s.pdf  %d KB  warnings=%d  (%.1fs)" %
          (stem, os.path.getsize(pdf) // 1024, len(warns), time.time() - t0))
    for w in warns[:8]:
        print("     warn:", w.strip()[:120])
    return True


if __name__ == "__main__":
    pdftex = find_pdflatex()
    print("pdflatex:", pdftex)
    ok = True
    for stem in (sys.argv[1:] or ["cs_main"]):
        ok &= compile_one(pdftex, stem)
    sys.exit(0 if ok else 1)
