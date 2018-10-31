(TeX-add-style-hook
 "torsten_manual"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("amsbook" "11pt" "reqno" "oneside")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "letterpaper" "width=6.5in" "height=9in") ("hyperref" "colorlinks=true" "citecolor=MRGGreen" "urlcolor=MRGGreen" "linkcolor=MRGGreen") ("mdframed" "framemethod=TikZ" "skipabove=10pt" "skipbelow=10pt" "backgroundcolor=black!5" "roundcorner=4pt" "linewidth=1pt") ("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("minted" "newfloat")))
   (add-to-list 'LaTeX-verbatim-environments-local "minted")
   (add-to-list 'LaTeX-verbatim-environments-local "lstlisting")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "lstinline")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
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
    "minted"
    "caption")
   (TeX-add-symbols
    '("subtitle" 1)
    "mrgsubtitle"
    "mrgproject"
    "mrgtitle")
   (LaTeX-add-labels
    "sec:org18a62b9"
    "sec:org5a2ec37"
    "sec:orgac7ffc0"
    "sec:org1eb3894"
    "sec:orgc5c13c8"
    "sec:org64ca667"
    "sec:orge6241bf"
    "sec:org12054b0"
    "sec:org9761c32"
    "sec:org8eee07f"
    "sec:org99655c0"
    "sec:orge98da5a"
    "sec:orgeb453ce"
    "sec:org8b2665b"
    "sec:org7c67c29"
    "sec:org0c7345a"
    "sec:org938f9d0"
    "sec:org42b28e6"
    "sec:org5446901"
    "sec:org5784b2b"
    "sec:orga98b8e4"
    "sec:orgdeb2e94"
    "sec:org2ace732"
    "sec:orgfe80ec5"
    "sec:orgc94c7ab"
    "sec:org3021bee"
    "sec:org0179d0f"
    "sec:orga428dc1"
    "sec:org07fd4c5"
    "sec:org896be65"
    "sec:orgba866ac"
    "sec:orgdb8fe1c"
    "sec:org463e26a"
    "sec:orgfbe159b"
    "sec:org06b801a"
    "sec:org15d2b1f"
    "sec:orgde79636"
    "sec:org8dad922"
    "fig:orgf43a1d6"
    "sec:org2285c7e"
    "tab:orgf256a50"
    "fig:org8b179c7"
    "fig:orgbb8c03a"
    "fig:orgb02d4bd"
    "fig:org188eedc"
    "sec:orge488035"
    "sec:orgc24bc1a"
    "sec:orgd0052cb"
    "sec:orge8e3f4a"
    "sec:org231bc2a"
    "sec:orga0b89a3"
    "sec:orgd45884c"
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

