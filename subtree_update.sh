#!/bin/bash  
echo "Removing previous submodules"
git rm -rf cmdstan/stan/lib/stan_math/stan/math/torsten cmdstan/stan/lib/stan_math cmdstan/stan cmdstan stanc3 && git commit -am "rm previous container"
rm -rf cmdstan stanc3

if [ $# -eq 0 ]
then
    cmdstan_branch="torsten-develop"
    stan_branch="torsten-develop"
    math_branch="torsten-develop"
    torsten_math_branch="develop"
else
    cmdstan_branch=$1
    stan_branch=$2
    math_branch=$3
    torsten_math_branch=$4
fi

echo "Adding cmdstan ${cmdstan_branch} branch"
git subtree add --prefix=cmdstan --squash cmdstan ${cmdstan_branch}
git rm -rf cmdstan/stan && git commit -am "rm stan"
echo "Adding stan ${stan_branch} branch"
git subtree add --prefix=cmdstan/stan --squash stan ${stan_branch}
git rm -rf cmdstan/stan/lib/stan_math && git commit -am "rm math"
echo "Adding math ${math_branch} branch"
git subtree add --prefix=cmdstan/stan/lib/stan_math --squash math ${math_branch}
git rm -rf cmdstan/stan/lib/stan_math/stan/math/torsten && git commit -am "rm torsten_math"
echo "Adding torsten_math ${torsten_math_branch} branch"
git subtree add --prefix=cmdstan/stan/lib/stan_math/stan/math/torsten --squash torsten_math ${torsten_math_branch}


# include stanc3 for debugging
echo "Adding stanc3 torsten-develop branch"
git subtree add --prefix=stanc3 --squash stanc3 torsten-develop
