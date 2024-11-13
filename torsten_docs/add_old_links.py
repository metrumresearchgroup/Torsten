"""
Adds an entry to the table listing old versions on the documentation homepage.
"""

from operator import le
import pathlib
import sys

HERE = pathlib.Path(__file__).parent
INDEX = HERE / 'src' / 'index.qmd'


ROW = "| {major}.{minor}    | [html](https://mc-stan.org/docs/{major}_{minor}/reference-manual/) [pdf](https://mc-stan.org/docs/{major}_{minor}/reference-manual-{major}_{minor}.pdf) | [html](https://mc-stan.org/docs/{major}_{minor}/stan-users-guide/) [pdf](https://mc-stan.org/docs/{major}_{minor}/stan-users-guide-{major}_{minor}.pdf) | [html](https://mc-stan.org/docs/{major}_{minor}/cmdstan-guide/) [pdf](https://mc-stan.org/docs/{major}_{minor}/cmdstan-guide-{major}_{minor}.pdf) | [html](https://mc-stan.org/docs/{major}_{minor}/functions-reference/) [pdf](https://mc-stan.org/docs/{major}_{minor}/functions-reference-{major}_{minor}.pdf) |"


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python add_old_links.py <version>')
        sys.exit(1)
    major = sys.argv[1]
    minor = sys.argv[2]
    row = ROW.format(major=major, minor=minor)
    index = INDEX.read_text()
    if row in index:
        print('Row already exists')
        sys.exit(1)

    new_index = []
    added = False
    for line in index.splitlines():
        new_index.append(line)
        if not added and line.startswith('|---------|'):
            new_index.append(row)
            added = True

    if not added:
        print('Failed to add row')
        sys.exit(1)

    INDEX.write_text('\n'.join(new_index))

