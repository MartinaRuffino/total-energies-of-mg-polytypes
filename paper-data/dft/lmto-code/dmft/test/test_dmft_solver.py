#!/bin/env python
import os, sys, subprocess


def h5diff(f1,f2,d=0) :
    """
    this subroutine run h5diff because we do not want to use h5py
    """
    f=open('out.h5diff','w')
    if(d>0) :
        subprocess.call('h5diff -d {delta} {0} {1} '.format(f1,f2,delta=d),shell=True,stdout=f)
    else  :
        subprocess.call('h5diff  {0} {1} '.format(f1,f2),shell=True,stdout=f)
    f.close()
    f=open('out.h5diff','r')
    message=f.read()

    FAILED=False
    error_clue=['not','error','differences found','Some objects are not comparable','unable']
    for err in error_clue :
        if (err in message)  :
            FAILED=True
    return FAILED

def test_solver_input(mpicmd) :
    cmd=mpicmd+' lmfdmft nio  --rs=1,0 -vnkabc=3 --ldadc=0 --job=1 '
    subprocess.call(cmd,shell=True)
    FAILED=h5diff('solver_input.h5','solver_input_bm.h5',d=1e-8)
    if(FAILED) :
        print('FAILED : solver_input.h5 and solver_input_bm.h5 are not equl')
    else :
        print('PASSED : solver_input.h5 is equal to solver_input_bm.h5')



if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='test switch cix ')
    p.add_argument('--mpi',type=int,nargs='?',default=1)
    p.add_argument('-v','--verbosity',action='store_true',help='print output of lmfdmft')
    args = p.parse_args()

    mpicmd=''
    if(args.mpi>1) :
        mpicmd='mpirun -np '+str(args.mpi)+' '
    dirdata = os.path.dirname(sys.argv[0])

    subprocess.call('cp  '+dirdata+'/nio/* .',shell=True)
    subprocess.call('cp   '+dirdata+'/../../gwd/test/nio.code2/{basp,rst,site}.nio  .',shell=True)
    test_solver_input(mpicmd)
    sys.exit()
