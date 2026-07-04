import subprocess, os, sys
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

tex_dir = r'C:\Users\高帅东\Desktop\ICLR2027_Submission_Clean'
tex_file = os.path.join(tex_dir, 'main.tex')
pdflatex = r'C:\Users\高帅东\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe'
bibtex = r'C:\Users\高帅东\AppData\Local\Programs\MiKTeX\miktex\bin\x64\bibtex.exe'

steps = [
    ('pdflatex #1', [pdflatex, '-interaction=nonstopmode', '-output-directory', tex_dir, tex_file]),
    ('bibtex', [bibtex, os.path.join(tex_dir, 'main')]),
    ('pdflatex #2', [pdflatex, '-interaction=nonstopmode', '-output-directory', tex_dir, tex_file]),
    ('pdflatex #3', [pdflatex, '-interaction=nonstopmode', '-output-directory', tex_dir, tex_file]),
]

for name, cmd in steps:
    result = subprocess.run(cmd, capture_output=True, cwd=tex_dir, timeout=120)
    ok = result.returncode == 0
    if not ok:
        with open(os.path.join(tex_dir, 'main.log'), 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if '!' in line and not line.strip().startswith('('):
                    print(f'  {line.rstrip()}')
    print(f'{name}: OK' if ok else f'{name}: FAIL')

pdf = os.path.join(tex_dir, 'main.pdf')
if os.path.exists(pdf):
    print(f'PDF: {os.path.getsize(pdf)/1024:.0f} KB')
else:
    print('ERROR: no PDF')
