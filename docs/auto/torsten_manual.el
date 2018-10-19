(TeX-add-style-hook
 "torsten_manual"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("amsbook" "11pt" "reqno")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "letterpaper" "width=6.5in" "height=9in") ("hyperref" "colorlinks=true" "citecolor=MRGGreen" "urlcolor=MRGGreen" "linkcolor=MRGGreen") ("mdframed" "framemethod=TikZ" "skipabove=10pt" "skipbelow=10pt" "backgroundcolor=black!10" "roundcorner=10pt" "linewidth=1pt") ("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem")))
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-environments-local "minted")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (TeX-run-style-hooks
    "latex2e"
    "mrgtitlepage"
    "amsbook"
    "amsbook11"
    "imakeidx"
    "geometry"
    "graphicx"
    "pdfpages"
    "amssymb"
    "epstopdf"
    "xcolor"
    "hyperref"
    "courier"
    "listings"
    "siunitx"
    "booktabs"
    "mdframed"
    "subcaption"
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
    "mrgtitle")
   (LaTeX-add-labels
    "sec:org3550bd3"
    "sec:org9568a75"
    "sec:org85a9989"
    "sec:org3bcb751"
    "sec:org72e48d4"
    "sec:org6d3af50"
    "sec:orgce5f3c7"
    "sec:orge4fbaa2"
    "sec:orgee4f14a"
    "sec:org3ca4883"
    "sec:org959d934"
    "sec:orgec71834"
    "sec:orge3b95cb"
    "sec:org8dcf413"
    "sec:org64b05a8"
    "sec:orgd4036b3"
    "sec:org9139461"
    "sec:org8590bc5"
    "sec:orga4360cd"
    "sec:org0c02946"
    "sec:org734392b"
    "fig:orgf68d5c5"
    "sec:orga93f567"
    "tab:orgd533f27"
    "fig:orgc8b4df2"
    "fig:org9c3e632"
    "fig:orgf3926cc"
    "fig:org4e3d4b1"
    "sec:orgc879054"
    "sec:org49ef30f"
    "sec:orgc229a17"
    "sec:orgdfde26a"
    "eq:FK")
   (LaTeX-add-bibliographies
    "torsten")
   (LaTeX-add-index-entries
    "One Compartment Model"
    "Two Compartment Model"
    "General linear model"
    "General ODE Model"
    "Mixed ODE Model"
    "Friberg-Karlsson Model")
   (LaTeX-add-amsthm-newtheorems
    "example"
    "remark")
   (LaTeX-add-xcolor-definecolors
    "MRGGreen"))
 :latex)

