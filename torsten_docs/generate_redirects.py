import pathlib

HERE =  pathlib.Path(__file__).parent

TEMPLATE = """<!DOCTYPE html>
<html lang="en-US">
  <meta charset="utf-8">
  <title>Redirecting&hellip;</title>
  <link rel="canonical" href="{REDIRECT_TO}">
  <script>location="{REDIRECT_TO}"</script>
  <meta http-equiv="refresh" content="0; url={REDIRECT_TO}">
  <meta name="robots" content="noindex">
  <h1>Redirecting&hellip;</h1>
  <a href="{REDIRECT_TO}">Click here if you are not redirected.</a>
</html>
"""

redirects = (HERE/'redirects.txt').read_text().strip()

for line in redirects.split('\n'):
    if line.startswith('#'):
        continue

    FROM, TO = line.split()
    content = TEMPLATE.format(REDIRECT_TO=TO)

    path = FROM.removeprefix("https://mc-stan.org/")

    file = HERE / path
    if file.exists():
        print(f"Skipping {file} as it already exists")
        continue

    file.write_text(content)
