#!/usr/bin/python

# update all html pages in a document
# bookdown-generated HTML page content is in <section> element
# for each file, insert link to latest at beginning of page section

import glob
import os
import sys
import contextlib

redirect_template = '''
<div>
This is an old version, <a href="https://mc-stan.org/docs/{DIRNAME}/{FILENAME}">view current version</a>.
</div>
'''
canonical_template = '<link rel="canonical" href="https://mc-stan.org/docs/{DIRNAME}/{FILENAME}" />'

@contextlib.contextmanager
def pushd(new_dir):
    previous_dir = os.getcwd()
    os.chdir(new_dir)
    yield
    os.chdir(previous_dir)

def main():

    dir = sys.argv[1]
    dirname = os.path.split(dir)[-1]
    print('dir: {}'.format(dirname))
    with pushd(dir):
        files = [x for x in glob.glob("*.html")]
        for file in files:
            print(file)
            canonical = canonical_template.format(DIRNAME=dirname,FILENAME=file)
            redirect_msg = redirect_template.format(DIRNAME=dirname,FILENAME=file)

            with open(file) as infile:
                with open('tmpfile', 'w') as outfile:
                    for line in infile:
                        outfile.write(line)
                        if '<head>' in line.strip():
                            outfile.write(canonical)
                        if '<main class="content" id="quarto-document-content">' in line.strip():
                            outfile.write(redirect_msg)
                            outfile.write('')
            os.replace('tmpfile',file)

if __name__ == '__main__':
    main()
