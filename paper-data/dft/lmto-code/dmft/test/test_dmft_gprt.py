#!/usr/bin/env python
import os ,sys, subprocess
import numpy as np

def readgkloc(f,nk,nom,norb,diag,comments='#') :
    gk=np.loadtxt(f,comments=comments).T
    W=gk[0,0:nom]
    gk=gk[1::2]+1j*gk[2::2]
    try :
         if diag :
             gk=gk.reshape(norb,nk,nom)
         else :
            gk=gk.reshape(norb,norb,nk,nom)
    except :
        message='''
{0}  does not have the expected shape
got : {1}
expected : {2},{3}
=======>FAILED
'''.format(f,gk.shape,nk,nom)
        raise Exception(message)

    return W,gk


def test_mode3_nio(nk,mpicmd,verbosity=False,refpath='.'):
    cmd="""
set -x
cp {refpath}/lscoq/* ./ ;
cp sig.inp sig.lscoq
{mpicmd} lmfdmft lscoq --mu=2  --rs=1,0 -vnkabc={nk} --ldadc=82.2 --job=1 --gprt~rdsig=sig~mode=5
""".format(mpicmd=mpicmd,nk=nk,refpath=refpath)

    if verbosity :
        subprocess.call(cmd,shell=True)
    else :
        log=open('out.lmfdmft','w')
        subprocess.call(cmd,shell=True,stdout=log)

    nk=nk**3 # regular 3d kmesh
    nom=999
    norb=5
    W,gkfs=readgkloc('gk.ref',nk,nom,norb,True,comments='%')
    W,gkfm=readgkloc('gkloc.lscoq',nk,nom,norb,False)
    iom=0
    ik=4
    equal=True
    for i in range(5) :
        if(equal) :
            equal=np.allclose(gkfs[i,:,:],gkfm[i,i,:,:],atol=1e-4)
    if(equal) :
        print('===>PASSED : --fullg  : gk_fullg_mpi  is equal to the reference')
    else :
        print('===>FAILED : --fullg  : gk_fullg_mpi  is not equal to the reference')


if __name__ == '__main__':
    # need to test mode 3,5,19,20 with mpi, several atom, spin resolved
    import argparse
    p = argparse.ArgumentParser(description='test switch cix ')
    p.add_argument('--mpi',type=int,nargs='?',default=1)
    p.add_argument('--nk','-nk',type=int,nargs='?',default=1)
    p.add_argument('--refpath',type=str,nargs='?',default=os.path.dirname(sys.argv[0]))
    p.add_argument('-v','--verbosity',action='store_true',help='print output of lmfdmft')
    args = p.parse_args()

    mpicmd=''
    if(args.mpi>1) :
        mpicmd='mpirun -np '+str(args.mpi)+' '

    test_mode3_nio(nk=args.nk, mpicmd=mpicmd, verbosity=args.verbosity, refpath=args.refpath)

    sys.exit()
