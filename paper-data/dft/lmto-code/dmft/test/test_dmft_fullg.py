#!/usr/bin/env python
import os ,sys, subprocess
import numpy as np
def test_mode3fullg_nio(nk,mpicmd,verbosity=False,refpath='.'):
    cmd="""
cp {refpath}/niso/* . ;
cp sig.inp sig.niso
{mpicmd} lmfdmft  niso -vso=1 -vnkabc={nk} --ldadc=82.2 --job=1 --gprt:mode=3:nom=20 --fullg
""".format(mpicmd=mpicmd,nk=nk,refpath=refpath)
    print(cmd)
    if verbosity :
        subprocess.call(cmd,shell=True)
    else :
        log=open('out.lmfdmft','w')
        subprocess.call(cmd,shell=True,stdout=log)
    Gref=np.loadtxt('gk.ref')
    Gtest=np.loadtxt('gkloc.niso')
    maxdiff = abs(Gref - Gtest).max()
    tol=1e-6; res = 'FAILED' if maxdiff > tol else 'PASSED'
    print('===> {res}: maxdiff between gkloc.niso and ref gk.ref ({diff}) within tol of {tol}.'.format(res=res,diff=maxdiff,tol=tol))

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='test switch cix ')
    p.add_argument('--mpi',type=int,nargs='?',default=1)
    p.add_argument('--nk','-nk',type=int,nargs='?',default=8)
    p.add_argument('--refpath',type=str,nargs='?',default=os.path.dirname(sys.argv[0]))
    p.add_argument('-v','--verbosity',action='store_true',help='print output of lmfdmft')
    args = p.parse_args()

    mpicmd=''
    if(args.mpi>1) :
        mpicmd='mpirun -np '+str(args.mpi)+' '

    test_mode3fullg_nio(args.nk, mpicmd, verbosity=args.verbosity, refpath=args.refpath)
