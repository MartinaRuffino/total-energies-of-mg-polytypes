#!/bin/env python2.7
from pytriqs.operators.util.op_struct import set_operator_structure
import pytriqs.operators.util as op
from pytriqs.operators.util.U_matrix import U_matrix,U_matrix_kanamori,spherical_to_cubic,cubic_names,subarray
from pytriqs.operators.util.hamiltonians import h_int_slater, h_int_kanamori
from pytriqs.gf import *
from pytriqs.applications.impurity_solvers.cthyb import Solver
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.gf.tools import *
import numpy as np
import csv,os

from impurity import *


def readinfo() :
    info={}
    f = h5py.File('solver_input.h5', 'r')

    beta = f['beta'][0]
    nomg = int(f['nomg'][0])
    W = np.array([ (2*n+1)*1j*np.pi/beta for n in range(nomg)])
    info['W']=W
    info['nsp']=f['nsp'][0]
    return info


def loadfile(f, mode='r'):
    fd = open(f, mode)
    s = fd.read()
    fd.close()
    return s

def readparam(f,verbose=True) :
    s=loadfile(f)
    s=s.split('\n')
    inputp={}
    p={}
    p["measure_G2_n_fermionic"]=0
    p["measure_G2_n_bosonic"]=0
    p['measure_G2_n_l']=0
    inputp['mode']="Full"

    intp=["n_warmup_cycles","n_cycles","length_cycle","measure_G2_n_fermionic","measure_G2_n_l","measure_G2_n_bosonic"]
    inti=["fit_max_moment","fit_start","fit_stop","n_l","NW"]

    floati=["U","J"]
    stri=["mode"]

    for l in s :
        if(l==''):
            break
        key,val=l.split(',')
        if(key in intp):
            p[key]=int(val)
            if mpi.is_master_node() and verbose:
                print(key,val)
        elif (key in inti) :
            inputp[key]=int(val)
            if mpi.is_master_node() and verbose:
                print(key,val)
        elif ( key in floati) :
            inputp[key]=float(val)
            if mpi.is_master_node() and verbose:
                print(key,val)
        elif (key in stri) :
            inputp[key]=str(val)
            if mpi.is_master_node() and verbose:
                print(key,val)
        else :
            if mpi.is_master_node() and verbose:
                print(str(key)+' ignored')
                continue

    if((p["measure_G2_n_fermionic"]>0) and (p["measure_G2_n_bosonic"]>0)) :
        p['measure_G2_iw_ph'] = True
    else :
        p['measure_G2_iw_ph'] = False

    for key in p :
        if mpi.is_master_node() and verbose:
            print('{0} \t {1}'.format(key,p[key]))
    return p,inputp


def G0_triqs(imp, G0, spin_names, orbital_names, orbitals_index) :

    n_ind_orb = len(orbital_names)
    for ispin, spin in enumerate(spin_names) :
        if imp.offdiag :

            for orb1, orb2 in product(range(n_ind_orb), range(n_ind_orb)) :
                #Spin polarized ?
                if imp.nsp == 1 :
                    G0[spin][orb1, orb2].data[imp.nomg:] = imp.g0[0, orb1, orb2,:]
                    G0[spin][orb1, orb2].data[:imp.nomg] = np.conjugate(imp.g0[0, orb1, orb2, ::-1])

                else :
                    G0[spin][orb1, orb2].data[imp.nomg:] = imp.g0[ispin, orb1, orb2,:]
                    G0[spin][orb1, orb2].data[:imp.nomg] = np.conjugate(imp.g0[ispin, orb1, orb2, ::-1])

        else :
            for orb in range(n_ind_orb) :
                if imp.nsp == 1 :
                    G0[spin+'_'+orbital_names[orb]].data[imp.nomg:,0,0] = imp.g0[0, orb, orb, :]
                    G0[spin+'_'+orbital_names[orb]].data[:imp.nomg,0,0] = np.conjugate(imp.g0[0, orb, orb, ::-1])

                else :
                    G0[spin+'_'+orbital_names[orb]].data[imp.nomg:,0,0] = imp.g0[ispin, orb, orb, :]
                    G0[spin+'_'+orbital_names[orb]].data[:imp.nomg,0,0] = np.conjugate(imp.g0[ispin, orb, orb, ::-1])


def feedsig(imp,sig_tab, Sig, spin_names, orbital_names, orbitals_index,legendre=False) :
    """
    feed from triqs to sig from imp.
    shift : takes just positive matsubara freq
    remark : it can be used to any quantity which has the same structure as sig
    """

    n_ind_orb = len(orbital_names)

    if legendre :
        len_mesh=len(Sig.mesh)
        shift=0
    else :
        len_mesh=len(Sig.mesh)/2
        shift=len_mesh



    for ispin, spin in enumerate(spin_names) :
        if imp.offdiag :

            for orb1, orb2 in product(range(n_ind_orb), range(n_ind_orb)) :
                #Spin polarized ?
                if imp.nsp == 1 :
                    sig_tab[0, orb1, orb2,:] += 0.5 * Sig[spin][orb1, orb2].data[shift:]
                else :
                    sig_tab[ispin, orb1, orb2,:] = Sig[spin][orb1, orb2].data[shift:]


        else :
            for orb in range(n_ind_orb) :
                if imp.nsp == 1 :
                    sig_tab[0, orb, orb, :] +=  0.5 * Sig[spin+'_'+orbital_names[orb]].data[shift:,0,0]
                else :
                    sig_tab[ispin, orb, orb, :] = Sig[spin+'_'+orbital_names[orb]].data[shift:,0,0]


def H(imp,orbital_name,spin_names,H_int) :
    from pytriqs.operators import *
    E = imp.matrixform(imp.Eimp)
    n_ind_orb = E.shape[1]
    fops = []
    H0 = H_int*0
    for ispin, spin in enumerate(spin_names) :
        if imp.offdiag :
            for orb1, orb2 in product(range(n_ind_orb), range(n_ind_orb)) :
                fops.append( (spin,orbital_name[orb1]))
                fops.append( (spin,orbital_name[orb1]))

                if imp.nsp == 1 :
                    H0 += E[0,orb1,orb2] * c_dag(spin,orbital_name[orb1]) * c(spin,orbital_name[orb2])
                else :
                    H0 += E[ispin,orb1,orb2] * c_dag(spin,orbital_name[orb1]) * c(spin,orbital_name[orb2])
        else :
            for orb1 in range(n_ind_orb) :
                fops.append((spin+'_'+orbital_name[orb1],0))
                if imp.nsp == 1 :
                    H0 += E[0,orb1,orb1] * n(spin+'_'+orbital_name[orb1],0)
                else :
                    H0 += E[ispin,orb1,orb1] * n(spin+'_'+orbital_name[orb1],0)
    return H0, fops


def solverH1(imp,verbose=False) :
    p,para=readparam('para_ctqmc'+str(imp.imp_ind)+'.dat',verbose=verbose)
    # p is the parameters used by the solver
    #para is the parameters used elsewhere
    beta=imp.beta


    U=para['U']
    J=para['J']
    typeH=para['mode']
    n_l=para['n_l']

    nomg=imp.nomg

    spin_names = ['up','down']

    cubic_name = cubic_names(imp.L)

    T = spherical_to_cubic(imp.L)
    Um = U_matrix(imp.L, U_int=U, J_hund=J, basis='other', T=T)


    #Is all the bands included in the dmft ?
    # A band is concedered as included is in the cix matrix, the diag term is non zero
    orbital_names=[]
    for ib in range(len(cubic_name)):
        if imp.sigind[0,ib,ib] > 0 :
            orbital_names.append(cubic_name[ib])

    orbital_index=[ cubic_name.index(orb) for orb in orbital_names]
    if mpi.is_master_node():
        print("bands select are : ")
        for ib,b in enumerate(orbital_names) :
            print('orbital {0} : {1}'.format(orbital_index[ib],b))
    from pytriqs.atom_diag import *
    Um=subarray(Um,len(Um.shape)*[ orbital_index ])
    H_int=h_int_slater(spin_names, orbital_names, Um, off_diag=imp.offdiag)
    H0,fops = H(imp,orbital_names,spin_names,H_int)
    print(H0)
    print(fops)
    # === construction of Gf

    gf_struct = set_operator_structure(spin_names, orbital_names, imp.offdiag)
    ad_r=AtomDiag(H_int, fops)
    ad0=AtomDiag(H0, fops)


    G_iw = atomic_g_iw(ad_r, beta, gf_struct, nomg)
    G0_iw = atomic_g_iw(ad0, beta, gf_struct, nomg)
    G0_triqs(imp, G0_iw, spin_names, orbital_names, orbital_index)

    Sig=dyson(G0_iw=G0_iw,G_iw=G_iw)
    imp.sig=np.zeros((imp.nsp, len(orbital_names), len(orbital_names), imp.nomg),dtype=complex)
    feedsig(imp,imp.sig, Sig, spin_names, orbital_names, orbital_index)


def H1() :
    f = h5py.File('solver_input.h5', 'r')
    nicix=f['nicix'][0]
    info=readinfo()
    sig_final=[info['W'].imag]
    gl_final=[]

    for icix in range(nicix) :
        imp=impurity(icix)
        imp.readh5()
        solverH1(imp)
        sig=imp.compress_form(imp.sig)

        for sp in range(info['nsp']) :
            for sig_orb in sig[sp] :
                sig_final.append(sig_orb.real)
                sig_final.append(sig_orb.imag)


    np.savetxt('sig.inp',np.array(sig_final).T,fmt='%1.8f')






if __name__ == '__main__':
    H1()
