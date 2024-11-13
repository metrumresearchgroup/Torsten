#!/usr/bin/env bash

### USAGE
# add_links.sh last_docs_version_dir
# add_links.sh 2_21

directory="$1"

mkdir docs/"$directory"_tmp

(cd docs/"$directory"; tar cf - .) | (cd docs/"$directory"_tmp; tar xvf - )

python3 link_to_latest.py docs/"$directory"_tmp/stan-users-guide "Stan User's Guide"
python3 link_to_latest.py docs/"$directory"_tmp/reference-manual "Stan Reference Manual"
python3 link_to_latest.py docs/"$directory"_tmp/functions-reference "Stan Functions Reference"
python3 link_to_latest.py docs/"$directory"_tmp/cmdstan-guide "CmdStan User's Guide"

mv docs/"$directory" docs/"$directory"_bak
mv docs/"$directory"_tmp/ docs/"$directory"

rm -rf docs/"$directory"_bak
