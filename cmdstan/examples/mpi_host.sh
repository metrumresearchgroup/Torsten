#!/bin/bash
rm ~/.ssh/config
rm hostfile
rm hostfile_*
qstat -f | grep -oh "ip[0-9\-]*" | tee >hostfile >(while read -r line ; do echo -e "Host $line\n   StrictHostKeyChecking no\n   UserKnownHostsFile=/dev/null" >> ~/.ssh/config; done)
declare -i y=$(wc -l < hostfile)
declare -i z=$1
declare -i x=$((y/z))
echo "split "$y" hosts for "$z" chains, with "$x" hosts per chain"
split -d -l $x -a 1 hostfile hostfile_

