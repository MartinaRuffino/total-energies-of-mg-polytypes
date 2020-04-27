#!/bin/bash

set -ex

docker pull dmt4/archlinux-lm-min # necessary to update existing local image
docker run --rm -ti -u `id -u`:`id -g` -e LC_CTYPE="en_US.UTF-8" -e NINJA_STATUS='[%c %r %s/%t %p] ' -v `pwd`:/w:rw -w /w dmt4/archlinux-lm-min ./quickbuild.sh
