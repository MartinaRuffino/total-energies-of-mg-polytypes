#!/bin/bash
#SBATCH -C scarf18
#SBATCH -t 0-05:00:00
#SBATCH -n 96
#SBATCH -N 4

# assume the path and module setup is satisfactory
#module load intel/18.0.2 intel/mpi/18.0.2 intel/mkl/18.0.2 python/3.6.5-intel-17.0.2 scipy/1.0.0-gcc-4.8.5-python-3.6.3
#export PATH=~/lm/bld/:$PATH


#---archer config---
# #PBS -l select=3
# #PBS -l walltime=5:59:0
# #PBS -A c01-band
# #PBS -m ae
# #PBS -M jerome.jackson@stfc.ac.uk
# 
# module swap PrgEnv-cray PrgEnv-intel
# module load questaal/vv
# module load anaconda/python3 # this is run on the service nodes
# 
# cd $PBS_O_WORKDIR

elements='ag al ar as au b ba be bi br c ca cd cl co cr cs cu f fe ga ge h he hf hg i in ir k kr li lu mg mn mo n na nb ne ni o os p pb pd po pt rb re rh rn ru s sb sc se si sn sr ta tc te ti tl v w xe y zn zr'

if [ ! -e $1 ]
then
	elements=$1
fi

for i in $elements
do
	./init2ctrl.py $i

	if [ -e ctrl.$i ]
	then
		nmin=`lmf --quit=ham $i |grep BZMESH |awk -vs=$SLURM_NTASKS '{print ($2>s ? s : $2)}'`
		for j in 0.94 0.96 0.98 1.00 1.02 1.04 1.06
		do 
			fn=scflog_$j.$i

			# if either the file is absent, or it didn't converge
                        if ( ! grep -E '^c|^x' $fn &>/dev/null )
			then 
				# -IB (option to intel mpirun) selects the OFA infiniband fabric
				mpirun -IB -np $nmin lmf -vsc=$j $i >$fn
				rm -f {mixm,rst}.$i
			else
				echo "$fn ok"
			fi
		done

		./fit.py $i >> eos.txt
		rm -f {log,mixm,moms,rst,save,atm,basp,wkp}.$i

	fi
done
