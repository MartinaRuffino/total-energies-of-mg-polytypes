#!/bin/bash -l

# Prepare the directory end environment for a test and runn it there.
#
# Arguments:
#   1  : path where the test should be run
#   2  : number of cores required by the job, not used here but useful for queuing system integration
#   3* : the command line to be executed

rund="$1"

mkdir -p "$rund"
exec >"$rund"/.testcmd.log 2>&1
set -x

np=$2

cmd0="$3"
shift 3
cmd="$(pwd)/$cmd0 $*"

name="${rund##*/}"; name="${name%.d}"
log="../${name}.log"

cwd=$(pwd)

export PATH="$cwd:$PATH"

grep=grep
sed=sed
# use sane utils on the insane system
[ "$(uname)" != Darwin ] || grep=ggrep
[ "$(uname)" != Darwin ] || sed=gsed


mpirun="$($grep -oPe '^\s*mpirun\s*=(.*)#?'  build.ninja | head -n 1 | $sed s:'^\s*mpirun\s*=\s*':'':)"
qcmd="$($grep -oPe '^\s*qcmd\s*=(.*)#?' build.ninja | head -n 1 | $sed -e s:'^\s*qcmd\s*=\s*':'': -e s:'\bNPROC\b':$np:g -e s:'\bNAME\b':ql-$name:g -e s:'\bOUTFILE\b':'.testcmd.log':g)"
qhdr="$($grep -oPe '^\s*qhdr\s*=(.*)#?' build.ninja | head -n 1 | $sed -e s:'^\s*qhdr\s*=\s*':'':)"

cd $rund


printf "#!/bin/bash -l
# set -x

cd $(pwd)
exec >$log 2>&1

t0=\$(date +%%s); echo \"t0: \$t0 s since époque\"

$qhdr

export PATH=\"$cwd:\$PATH\"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OMPI_MCA_hwloc_base_binding_policy=none

$cmd
rc=\$?

t1=\$(date +%%s); echo \"t1: \$t1 s since époque\"
echo \"walltime: \$((t1-t0)) s\"
exit \$rc
" > .run.sh

chmod u+x .run.sh

$qcmd ./.run.sh

exit 0

