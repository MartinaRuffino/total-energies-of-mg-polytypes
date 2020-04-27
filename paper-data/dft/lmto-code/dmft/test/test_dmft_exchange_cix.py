#!/usr/bin/env python
import os ,sys, subprocess
import numpy as np




def loadfile(f, mode='r'):
    fd = open(f, mode)
    s = fd.read()
    fd.close()
    return s

def readeimp(f) :
    s=loadfile(f)
    s=s.split('\n')
    set='-0123456789.,'
    lstE=[]
    for i in [0,1,3] :
        E=s[i].split('#')[0] # not take into account the comment
        E=''.join(c for c in E if c in set) # delete the first part Eimp=
        E=E.replace('[','')
        E=E.replace(']','')

        lstE.append([float(a) for a in E.split(',')])

    return np.array(lstE)

def test_tab_cix(tab1,tab2) :
    #file structure : mesh Re1 Im1 .....
    # 1:10 raws cix 1
    # 11:21 raws cix 2
    #In the future, could be deduce from indmfl.ext ...
    dcix11=tab1[1:11]
    dcix12=tab2[11:]
    dcix21=tab1[11:]
    dcix22=tab2[1:11]
    equal=(np.allclose(dcix11,dcix12,atol=1e-6) and  np.allclose(dcix21,dcix22,atol=1e-6))
    return equal


def test_normal(mpicmd,nk=2,verbosity=False) :
    #===========>  test normal mode<========================#
    f=open('cmd_test_normal_mode','w')
    # cix1 cix2
    cmd='set -xe ;cp sig.inp.temoin sig.inp    ; '
    cmd='cp indmfl1.ybco indmfl.ybco  ; '
    cmd+='cp dc1.ybco dc.ybco  ;'
    cmd+=mpicmd+' lmfdmft ybco  --rs=1,0 -vnkabc='+str(nk)+' --ldadc~fn=dc --job=1 ; '
    cmd+='cp delta.ybco delta1.ybco ;'
    cmd+='cp eimp1.ybco E1.ybco ;'
    f.write(cmd)
    f.write('\n')
    if(verbosity) :
        res=subprocess.call(cmd,shell=True)
    else :
        log=open('log_cix12','w')
        res=subprocess.call(cmd,shell=True,stdout=log)
        log.close()
    # 1 <-->2
    S=np.loadtxt('sig.inp.temoin')
    Smiroir=S.copy()
    Smiroir[:,1:11]=S[:,11:]
    Smiroir[:,11:]=S[:,1:11]
    np.savetxt('sig_miroir.inp',Smiroir,fmt="%1.8f")
    np.savetxt('sig.inp',Smiroir,fmt="%1.8f")
    cmd='cp indmfl2.ybco indmfl.ybco ; '
    cmd+='cp dc2.ybco dc.ybco ;'
    cmd+=mpicmd+' lmfdmft ybco  --rs=1,0 -vnkabc='+str(nk)+' --ldadc~fn=dc --job=1 ; '
    cmd+='cp delta.ybco delta2.ybco ;'
    cmd+='cp eimp1.ybco E2.ybco ;'
    f.write(cmd)

    if(verbosity) :
        subprocess.call(cmd,shell=True)
    else :
        log=open('log_cix21','w')
        subprocess.call(cmd,shell=True,stdout=log)
        log.close()
    f.close()
    # test if delta.inp are equal
    data1=np.loadtxt('delta1.ybco').T
    data2=np.loadtxt('delta2.ybco').T
    equal=test_tab_cix(data2,data1)
    if equal :
        print('==> PASSED  : normal mode test')
    else :
        print('==>  FAILED : delta.inp non equal')
    equal=True
    lstE1=readeimp('E1.ybco')
    lstE2=readeimp('E2.ybco')
    for iE in range(len(lstE1)) :
        #cix 1
        Ecix11=lstE1[iE][0:5]
        Ecix12=lstE2[iE][5:]
        if(not np.allclose(Ecix11,Ecix12) ) :
            equal=False
        #cix 2
        Ecix21=lstE1[iE][5:]
        Ecix22=lstE2[iE][0:5]
        if(not np.allclose(Ecix11,Ecix12) ) :
            equal=False

    if equal :
        print('==> PASSED : Eimp equal')
    else :
        print('==> FAILED : Eimp not equal')




def test_mode_5(mpicmd,nk=2,verbosity=True) :
    # # #===========> test mode 5 <========================#
    #  --gprt~rdsig=sig~mode=5  --fullg
        # cix1 cix2
    cmd='cp sig.inp.temoin sig.inp ; \n'
    cmd+='cp sig.inp sig.ybco ; \n'
    cmd+='cp indmfl1.ybco indmfl.ybco ; \n'
    cmd+='cp dc1.ybco dc.ybco ; \n'
    cmd+=mpicmd+' lmfdmft ybco -vnkabc='+str(nk)+' --rs=1,0  --ldadc~fn=dc --job=1  --gprt~rdsig=sig~mode=5  --fullg; \n'
    cmd+='cp gloc_1.ybco gloc_11 ;\n'
    cmd+='cp gkloc_1.ybco gkloc_11 ;\n'
    cmd+='cp gloc_2.ybco gloc_21 ;\n'
    cmd+='cp gkloc_2.ybco gkloc_21 ;\n'
    print(cmd)

    if(verbosity) :
        subprocess.call(cmd,shell=True)
    else :
        log=open('log_cix12_mode5','w')
        subprocess.call(cmd,shell=True,stdout=log)
        log.close()
    # 1<-->2
    S=np.loadtxt('sig.inp')
    Smiroir=S.copy()
    Smiroir[:,1:11]=S[:,11:]
    Smiroir[:,11:]=S[:,1:11]
    np.savetxt('sig_miroir.inp',Smiroir,fmt="%1.8f")
    np.savetxt('sig.inp',Smiroir,fmt="%1.8f")
    np.savetxt('sig.ybco',Smiroir,fmt="%1.8f")
    cmd='cp indmfl2.ybco indmfl.ybco ;\n'
    cmd+='cp dc2.ybco dc.ybco ;\n'
    cmd+=mpicmd+' lmfdmft ybco  --rs=1,0 -vnkabc='+str(nk)+' --ldadc~fn=dc --job=1    --gprt~rdsig=sig~mode=5  --fullg; \n'
    cmd+='cp gloc_1.ybco gloc_12 ;\n'
    cmd+='cp gkloc_1.ybco gkloc_12 ;\n'
    cmd+='cp gloc_3.ybco gloc_22 ;\n'
    cmd+='cp gkloc_3.ybco gkloc_22 ;\n' # In this config, the first 2 atoms are equivalent, the third is the second indenpent cix
    log=open('log_cix21_mode5','w')
    print(cmd)
    if(verbosity) :
        subprocess.call(cmd,shell=True)
    else :
        subprocess.call(cmd,shell=True,stdout=log)
        log.close()

    g11loc=np.loadtxt('gloc_11')
    g22loc=np.loadtxt('gloc_22')
    equal1=np.allclose(g11loc,g22loc,atol=1e-4)
    g12loc=np.loadtxt('gloc_12')
    g21loc=np.loadtxt('gloc_21')
    equal2=np.allclose(g12loc,g21loc,atol=1e-4)
    if (equal1 and equal2 ) :
         print('==> PASSED : gloc equal')
    else :
        print('==> FAILED : gloc equal')

    g11loc=np.loadtxt('gkloc_11')
    g22loc=np.loadtxt('gkloc_22')
    equal1=np.allclose(g11loc,g22loc,atol=1e-4)
    g12loc=np.loadtxt('gkloc_12')
    g21loc=np.loadtxt('gkloc_21')
    equal2=np.allclose(g12loc,g21loc,atol=1e-4)
    if (equal1 and equal2 ) :
         print('==> PASSED : gkloc equal')
    else :
        print('==> FAILED : gkloc equal')

if __name__ == '__main__':
    """
    In this test, there are 2 independent cix.
    we switch the two icix in indmftl and check if the results are equivalent.
     """
    # import argparse
    # p = argparse.ArgumentParser(description='test switch cix ')
    # p.add_argument('--mpi',type=int,nargs='?',default=1)
    # args = p.parse_args()
    import sys
    import argparse
    p = argparse.ArgumentParser(description='test switch cix ')
    p.add_argument('--mpi',type=int,nargs='?',default=1)
    p.add_argument('-nk','--nk',type=int,nargs='?',default=1)
    p.add_argument('-v','--verbosity',action='store_true',help='print output of lmfdmft')
    args = p.parse_args()
    verbose=args.verbosity

    nk=args.nk
    mpicmd=''
    if(args.mpi>1) :
        mpicmd='mpirun -np '+str(args.mpi)+' '

    nk=2

    dirdata = os.path.dirname(sys.argv[0]) + '/ybco'

    subprocess.call('cp  '+dirdata+'/* .',shell=True)
    test_normal(mpicmd,nk=nk,verbosity=verbose)
    test_mode_5(mpicmd,nk=nk,verbosity=verbose)
