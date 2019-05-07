# setupTorsten-dev
# v 1.6
# Run this bash file in the directory where you want to install cmdstan
# with the Torsten functions. 
#
# update 1.6: use v2.17.1 for CmdStan, torsten-master for Stan and Stan Math
# update 1.5: changed clone of math so that deletion is not required
# update 1.4: download torsten from metrumresearchgroup git repositary.
# update 1.3: download cmdStan v2.16.0.
# update dev: download dev version of Torsten instead of last release.
# update 1.2: download cmdStan v2.14.0 and the master version of 
#             torsten-stan and torsten-math.
# update 1.1: download last version of cmdStan/dev with which Torsten
#             was tested.

#!/bin/bash
git clone https://github.com/stan-dev/cmdstan.git
cd cmdstan
git checkout release/v2.17.1
git clone https://github.com/metrumresearchgroup/stan.git
cd stan
git checkout torsten-master
cd lib
cd stan_math
git clone https://github.com/metrumresearchgroup/math.git .
git checkout torsten-master

