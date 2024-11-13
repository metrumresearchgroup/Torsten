# Repository `docs`

Repository for the sources and published documentation set, versioned for each Stan minor release.

* The Stan User's Guide - example models and techniques for coding statistical models in Stan and using them to do inference and prediction.
* The Stan Reference Manual - specification for Stan language and core inference algorithms.
* The Stan Functions Reference - functions and distributions built into the Stan language.
* The CmdStan Guide - guide to the reference command-line interface.


## Repository directory structure

* `src` : directory of source files for Stan and CmdStan guides and reference manuals, each in its own named subdirectory:
    + `src/cmdstan-guide` - CmdStan Guide
    + `src/functions-reference` - Stan Functions Reference
    + `src/reference-manual` - Stan Reference Manual
    + `src/stan-users-guide` - Stan Users Guide

* `docs`: the directory `docs` on branch `master` is the [publishing source](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site) for the project pages site.  Whenever a verified member of the Stan organization pushes to `docs` on branch `master`,
GitHub (re)builds and (re)deploys the website.

## Documentation toolset

We use [Quarto](https://quarto.org/) to build the HTML website and standalone pdfs;
previously, we used [bookdown](https://github.com/rstudio/bookdown).
[Download quarto](https://quarto.org/docs/download/)
To build the pdf version of the docs, you will need to [install LaTeX](https://www.latex-project.org/get/) as well.

Quarto accepts [`.qmd`](https://quarto.org/docs/authoring/markdown-basics.html) source files
and uses the [Pandoc](https://pandoc.org) conversion engine.

## Scripts to build and maintain the docset

**`build.py`**

The program `build.py` convert the markdown files under `src` to html and pdf and populates the `docs` dir with the generated documentation.
Requires Python 3.7 or higher, due to call to `subprocess.run`, kwarg `capture_output`.
  + 2 required arguments:  <Major> <minor> Stan version, expecting 2 positive integer arguments, e.g. `2 28`
  + 2 optional arguments:  <format> <document>.  The output format is either `website` or `pdf`.  The document name corresponds to the name of the `src` subdirectory or `all`.


**Build script examples**

* `python build.py 2 35` - creates directory `docs/2_42` as needed; populates it will all generated documentation.
* `python build.py 2 35 website` - builds the docs website in `docs/2_42`.
* `python build.py 2 35 pdf functions-reference` - builds only the pdf version of the Stan functions reference,  resulting document is `docs/2_35/functions-reference-2_35.pdf`
* `python build.py 2 35 pdf all` - builds all pdfs from the Stan documentation set, resulting pdfs are in `docs/2_35`.


**Additional scripts**

The release process generates a new documentation set and adds links and redirects across the docset.

* `add_redirects.py` manages the redirects from unversioned links to the latest version.
* `link_to_latest.py` adds the "latest version" link into a docset.

The Stan Functions Reference contains HTML comments which describe the function signature for all functions.  The script `extract_function_sigs.py` is used to scrape these signatures into a plain text file.


## GitHub Pages

This repository uses
[GitHub Pages](https://docs.github.com/en/pages/getting-started-with-github-pages)
to serve the
[project pages](https://docs.github.com/en/pages/getting-started-with-github-pages/about-github-pages#project-pages-sites) site
with URL https://mc-stan.org/docs.
The publishing strategy is to serve the contents of the directory `docs` on branch `master`.
The `docs` directory contains an empty file named `.nojekyll` so that GitHub will treat the contents
as pre-generated HTML instead of trying to run [jekyll](https://jekyllrb.com).

