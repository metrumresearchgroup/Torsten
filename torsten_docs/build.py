#!/usr/bin/env python
"""
Build docs website and/or pdfs.

Command line arguments:

    arg 1: MAJOR version number
    arg 2: MINOR version number
    optional arg 3:  format, either "pdf" or "website"; defaults to both
    optional arg 4: pdf doc, defaults to "all", ignored if arg3 is "website"

Quarto builds in subdirectory of `src`, and then
resulting files are moved to directory named `docs/MAJOR_MINOR`
Document pdf has subtitle "Version MAJOR dot MINOR".
Directory and filenames have version string "MAJOR underscore MINOR".
"""

import glob
import os
import os.path
import shutil
import subprocess
import sys
import contextlib

all_docs = (
    "torsten-users-guide",
    "torsten-developers-guide",
)


@contextlib.contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    yield
    os.chdir(previous_dir)


def shexec(command):
    ret = subprocess.run(
        command,
        shell=True,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )
    if ret.returncode != 0:
        print("Command {} failed".format(command))
        print(ret.stderr)
        raise Exception(ret.returncode)


def safe_rm(fname):
    if os.path.exists(fname):
        os.remove(fname)


def make_pdfs(docspath, version, document):
    path = os.getcwd()
    srcpath = os.path.join(path, "src", document)
    with pushd(srcpath):
        pdfname = "".join([document, "-", version, ".pdf"])
        print("render {} as {}".format(document, pdfname))
        command = " ".join(["quarto render", "--output-dir _pdf", "--output", pdfname])
        shexec(command)
        outpath = os.path.join(srcpath, "_pdf", pdfname)
        safe_rm(os.path.join(docspath, pdfname))
        shutil.move(outpath, docspath)


def main():
    if sys.version_info < (3, 8):  # required by shutil.copytree
        print("requires Python 3.8 or higher, found {}".format(sys.version))
        sys.exit(1)
    global all_docs
    build_web = True
    build_pdfs = True
    docset = all_docs

    if len(sys.argv) > 2:
        stan_major = int(sys.argv[1])
        stan_minor = int(sys.argv[2])
    else:
        print("Expecting version number args MAJOR MINOR")
        sys.exit(1)

    stan_version = f"{stan_major}_{stan_minor}"
    os.environ["STAN_DOCS_VERSION"] = f"{stan_major}.{stan_minor}"
    os.environ["STAN_DOCS_VERSION_PATH"] = stan_version
    path = os.getcwd()
    docspath = os.path.join(path, "docs", stan_version)
    if not (os.path.exists(docspath)):
        try:
            os.makedirs(docspath)
        except OSError:
            print("Failed to create directory %s" % docspath)
            sys.exit(1)
        else:
            print("Created directory %s " % docspath)

    if len(sys.argv) > 3:
        if sys.argv[3] == "pdf":
            build_web = False
        elif sys.argv[3] == "website":
            build_pdfs = False
        else:
            print("Bad arg[3], should be 'pdf' or 'website'".format(sys.argv[3]))
            sys.exit(1)

    if len(sys.argv) > 4:
        if sys.argv[4] not in docset:
            print("Bad arg[4], should be one of %s" % " ".join(docset))
            sys.exit(1)
        docset = (sys.argv[4],)

    if len(sys.argv) > 5:
        print("Unused arguments:  %s" % " ".join(sys.argv[5:]))

    if build_web:
        print("render website")
        with pushd(os.path.join(path, "src")):
            command = "quarto render"
            shexec(command)
            shutil.copytree(
                "_website", docspath, copy_function=shutil.move, dirs_exist_ok=True
            )

    if build_pdfs:
        for doc in docset:
            make_pdfs(docspath, stan_version, doc)


if __name__ == "__main__":
    main()
