(TeX-add-style-hook
 "refactor"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("amsbook" "11pt" "reqno")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "left=1in" "top=1in" "bottom=1.75in" "right=1in" "headsep=.4in") ("hyperref" "colorlinks=true" "citecolor=MRGGreen" "urlcolor=MRGGreen" "linkcolor=MRGGreen") ("textpos" "absolute") ("mdframed" "framemethod=TikZ" "skipabove=10pt" "skipbelow=10pt" "backgroundcolor=gray!10" "roundcorner=5pt" "linewidth=1pt") ("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("minted" "newfloat")))
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
    '("subtitle" 1)
    "mrgsubtitle"
    "mrgproject"
    "mrgtitle"
    "mrgcategory"
    "mrgPI"
    "todaymetrum")
   (LaTeX-add-labels
    "sec:org7cf529a"
    "sec:org114f78d"
    "sec:org25b8612"
    "sec:orgcb0bad7"
    "sec:orgc04ea48"
    "sec:org9ad9d0d"
    "sec:org9171e50"
    "sec:org7353f4e"
    "sec:orgd1f0253"
    "sec:org9ac72a2")
   (LaTeX-add-amsthm-newtheorems
    "example"
    "remark"
    "codex")
   (LaTeX-add-xcolor-definecolors
    "MRGGreen"))
 :latex)

