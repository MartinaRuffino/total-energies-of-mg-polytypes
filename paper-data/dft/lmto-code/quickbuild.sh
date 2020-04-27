#!/bin/bash

set -ex

#[ -z "$LC_CTYPE" ] && export LC_CTYPE='en_US.UTF-8'
[ -z "$NINJA_STATUS" ] && export NINJA_STATUS='[%c %r %s/%t %p] '

date

declare -A t

[ -f quickbuild ] && rm -rf quickbuild
mkdir -p quickbuild
[ -f quickbuild/s ] && rm -rf quickbuild/s
[ -f quickbuild/m ] && rm -rf quickbuild/m
mkdir -p quickbuild/s quickbuild/m

cwd="$PWD"

lscpu

t[g]=`date +%s`

cd "$cwd/quickbuild/m"
if [ ! -f flags.mk ]; then
    ../../genflags.py gcc dbg openmpi pic > flags.mk
    sed -i -e s:'\s-g':'':g -e s:'\s-fbacktrace':'':g -e s:'\s-finit-\w*=\w*':'':g -e s:'\s-fstack-protector-all':'':g -e s:O1:O0:g flags.mk

    ../../configure.py
fi

cd "$cwd/quickbuild/s"
if [ ! -f flags.mk ]; then
    ../../genflags.py gcc dbg pic > flags.mk
    sed -i -e s:'\s-g':'':g -e s:'\s-fbacktrace':'':g -e s:'\s-finit-\w*=\w*':'':g -e s:'\s-fstack-protector-all':'':g -e s:O1:O0:g flags.mk

    ../../configure.py
fi

cd "$cwd"



t[n]=`date +%s`

ninja -C "$cwd/quickbuild/m"
ninja -C "$cwd/quickbuild/s"

set +x

t[f]=`date +%s`

printf "av times: configuring: %ds, compiling: %ds\n" $(((${t[n]} - ${t[g]})/2)) $(((${t[f]} - ${t[n]})/2))

date

if [ "$1" == "--img" ]; then
  ninja test-blm
  ls -l checks/
fi


