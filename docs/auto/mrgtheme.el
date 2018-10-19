(TeX-add-style-hook
 "mrgtheme"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("amsbook" "11pt" "reqno")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "left=1in" "top=1in" "bottom=1.75in" "right=1in" "headsep=.4in") ("hyperref" "colorlinks=true" "citecolor=green" "urlcolor=green" "linkcolor=MRGGreen") ("textpos" "absolute") ("mdframed" "framemethod=TikZ" "skipabove=10pt" "skipbelow=10pt" "backgroundcolor=gray!10" "roundcorner=5pt" "linewidth=1pt") ("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("minted" "newfloat")))
   (add-to-list 'LaTeX-verbatim-environments-local "minted")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "mrgtitlepage"
    "amsbook"
    "amsbook11"
    "imakeidx"
    "geometry"
    "xcolor"
    "hyperref"
    "datetime2"
    "textpos"
    "graphicx"
    "float"
    "bm"
    "lastpage"
    "fancyhdr"
    "pdfpages"
    "amssymb"
    "epstopdf"
    "siunitx"
    "booktabs"
    "subcaption"
    "caption"
    "mdframed"
    "inputenc"
    "fontenc"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "amsmath"
    "textcomp"
    "capt-of"
    "minted")
   (TeX-add-symbols
    "todaymetrum")
   (LaTeX-add-amsthm-newtheorems
    "example"
    "remark"
    "codex")
   (LaTeX-add-xcolor-definecolors
    "MRGGreen"))
 :latex)

