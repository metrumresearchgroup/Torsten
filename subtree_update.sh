#!/bin/bash  
git rm -rf cmdstan/stan/lib/stan_math/stan/math/torsten && git commit -am "rm torsten_math"
git rm -rf cmdstan/stan/lib/stan_math && git commit -am "rm math"
git rm -rf cmdstan/stan && git commit -am "rm stan"
git rm -rf cmdstan && git commit -am "rm cmdstan"
git rm -rf stanc3 && git commit -am "rm stanc3"
rm -rf cmdstan stanc3

git subtree add --prefix=cmdstan --squash cmdstan torsten-develop
git rm -rf cmdstan/stan && git commit -am "rm stan"
git subtree add --prefix=cmdstan/stan --squash stan torsten-develop
git rm -rf cmdstan/stan/lib/stan_math && git commit -am "rm math"
git subtree add --prefix=cmdstan/stan/lib/stan_math --squash math torsten-develop
git rm -rf cmdstan/stan/lib/stan_math/stan/math/torsten && git commit -am "rm torsten_math"
git subtree add --prefix=cmdstan/stan/lib/stan_math/stan/math/torsten --squash torsten_math develop

git subtree add --prefix=stanc3 --squash stanc3 torsten-develop
