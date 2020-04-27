#!/bin/bash -l

# from where to pick input files (abs path needed, using ploughman's abspath)
i="$(dirname $(readlink -m $0))"

# run all the stuff away from the root git repo so 'git stat' doe not show it all
r=dc-gm6-spw

# probably want to use "aprun" on cray
l='mpirun'

# number of lmf mpi processes for each task, mostly k-point parallel, for n > 10 see note below
n=16

# elements in order of decreasing difficulty
els="ga bi pd mn as sb b na pb f br p cu i ag cl tl ir pt sr ca au ni n ru cr li w ta sc zn c rh cd tc os re co zr te ti sn lu ge se y hf mo nb in mg cs v hg al fe h rb ba k be rn s kr xe si po o ar ne he"
# els="pd mn as sb b na pb f br p cu i ag cl tl ir pt sr ca au ni n ru cr li w ta sc zn c rh cd tc os re co zr te ti sn lu ge se y hf mo nb in mg cs v hg al fe h rb ba k be rn s kr xe si po o ar ne he"

# scalings
scs="0.94 0.96 0.98 1.00 1.02 1.04 1.06"

cwd="$PWD"

for e in $els; do
    for s in $scs; do

        w="$r/${e}-${s}.d"

        [ ! -e "$w" ] || continue

        mkdir -p "$w"
        cd "$w"

# customise as necessary or maybe better define templates in separate files, also look into using job arrays for places like archer which severely restrict pending job numbers.

        echo "#!/bin/bash -l
#SBATCH -J ${e}${s}-dc
#SBATCH -n $n
#SBATCH -N 1-1
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-socket=16
#SBATCH --mem-per-cpu 4026M
#SBATCH -p henrya
#SBATCH -o %x-%j.log

set -xe

module load hdf5/1.10.5/openmpi/4.0.1/intel/19.0.3 libxc/4.3.4/intel/19.0.3 python/2.7.16/intel/19.0.3

export PATH=\"$cwd:$i:\$PATH\"


cp \"$i\"/init.$e ./

init2ctrl.py -k \"$i\"/kpts.lmto.summary -l \"$l -n 1\" $e

# note: for n > 10 need to do the whole --quit=ham | grep BZMESH ... song and dance..
n=\$($l -n 1 lmf -vsc=$s --quit=ham ctrl.$e | grep BZMESH | awk -vs=$n '{print (\$2>s ? s : \$2)}')
$l -n \$n lmf -vsc=$s ctrl.$e &> scflog.$s.$e.log

rm -f mixm.* moms.* rst.* atm.*

" > q.sh

        sbatch q.sh
        cd "$cwd"
        echo $e $s
    done
done
