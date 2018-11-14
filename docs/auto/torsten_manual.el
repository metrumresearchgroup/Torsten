(TeX-add-style-hook
 "torsten_manual"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("amsbook" "11pt" "reqno" "oneside")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("geometry" "letterpaper" "width=6.5in" "height=9in") ("hyperref" "colorlinks=true" "citecolor=MRGGreen" "urlcolor=MRGGreen" "linkcolor=MRGGreen") ("mdframed" "framemethod=TikZ" "skipabove=10pt" "skipbelow=10pt" "backgroundcolor=black!5" "roundcorner=4pt" "linewidth=1pt") ("placeins" "section") ("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("minted" "newfloat")))
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
    "placeins"
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
    "sec:orgba15e9f"
    "sec:org365cdab"
    "sec:orga82e59c"
    "sec:org62552bb"
    "sec:org99d4b11"
    "sec:orge6a0e3f"
    "sec:org1ce2b16"
    "sec:org2d07e64"
    "sec:orgdd7196b"
    "sec:org8648219"
    "sec:orgc05530e"
    "sec:orgb30183b"
    "sec:org380dc5f"
    "sec:org966fbf1"
    "sec:org1ae8c06"
    "sec:org81fca5c"
    "sec:org9a969ec"
    "sec:orga13d796"
    "sec:org3064cbb"
    "sec:org2ff976f"
    "sec:orgd3bd3c1"
    "sec:org4a7cdf0"
    "sec:orgc26a125"
    "sec:org14893f4"
    "sec:org71a3795"
    "sec:org10ec6bb"
    "sec:orgaaf62c6"
    "sec:org666091e"
    "sec:org715546a"
    "sec:org90e24d2"
    "sec:org5b05b62"
    "sec:org3d2fe5d"
    "sec:org0c344fb"
    "sec:orgb735d69"
    "sec:org1298669"
    "sec:org337922d"
    "sec:org49e67cb"
    "org508cda2"
    "fig:orgbe6371c"
    "eq:onecpt"
    "sec:org0147929"
    "orge037615"
    "eq:twocpt"
    "sec:org657464c"
    "tab:orgb058d64"
    "fig:org0c822f0"
    "fig:org514c5c2"
    "fig:orgc22c74c"
    "fig:org85b2cc9"
    "sec:org15b373c"
    "sec:orge17e28a"
    "sec:orgb60e784"
    "sec:org6e723c6"
    "sec:orgd5d2c97"
    "orgd6e813b"
    "sec:orgba37f3d"
    "eq:FK"
    "fig:org9b72a2b"
    "sec:org8b99740"
    "sec:org6809fc1"
    "sec:org8d847b4"
    "sec:orgf03d923"
    "sec:org2d27c6c"
    "sec:org41e1a30"
    "sec:orgf446826"
    "fig:org8cfda29"
    "sec:orgd4e5fa2"
    "sec:org6596fdb"
    "tab:org9ac148e"
    "effCptModelMCMC"
    "effCptModelDens"
    "effCptModelPredictionsPK"
    "effCptModelPredictionsPD"
    "sec:org84102b2"
    "sec:org3648344"
    "sec:org87b0bf7"
    "tab:org997322d"
    "FKMCMC"
    "FKDens"
    "FKPredictions")
   (LaTeX-add-bibliographies
    "torsten")
   (LaTeX-add-index-entries
    "One Compartment Model"
    "Two Compartment Model"
    "General linear model"
    "General ODE Model"
    "Mixed ODE Model"
    "Friberg-Karlsson Model"
    "univariate integral"
    "linear interpolation")
   (LaTeX-add-amsthm-newtheorems
    "example"
    "remark")
   (LaTeX-add-xcolor-definecolors
    "MRGGreen"))
 :latex)

