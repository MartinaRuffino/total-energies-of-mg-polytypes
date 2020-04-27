import h5py
import numpy as np
from itertools import product

class impurity :
    """
    An impurity is a inequibalent cix block in lmfdmft.
    it's defined by :

    l      :    the Angular-Momentum Quantum Number.
                there are (2*l+1)^2 possible channels

    offdiag :   if false, offdiagonal element in term
                of orbital are not taken into account

    nsp    :    the number of independent spins
                if 1, spins are averaged
                if 2, spins are considered as independant.

    delta_LL' : hybridisation function

    E_LL''     : impurity energy level

    sigind :   tab from m1,m2 to chan.
               each pair m1,m2 is represented in the compressed from by a channel(chan).
               each m1 represent a band for the impurity solver.
               if the chan of m1,m1 is zero, m1 is not included in the impurity solver.
               if chan of m1,m2 is non zeros, off diag element are included in the impurity solver
               example for d orbitals xy, yz, z^2, xz, x^2-y^2 :

                  1 0 0 0 0
                  0 1 0 0 0
                  0 0 2 0 5
                  0 0 0 0 0
                  0 0 4 0 3
              xy,yz,z^2 and x^2-y^2 active and of diagonal element taken into account between z^2 and x^2-y^2.
              xy and yz are considered as equivalent i.e sig[chan=1] will be the mean of xy and yz.
     sigind_red : same as sigind for active channel i.e :
                  1 0 0 0
                  0 1 0 0
                  0 0 2 5
                  0 0 4 3
    lst_active : lst of diagonal index of sigind which are non zeros. in the example above, it's 0,1,2,3
    """

    def __init__(self,imp_ind) :
        self.imp_ind = imp_ind #0 < =  inmpurity_ind < =  nicix
        #read for h5 file
        self.L = 0
        self.sizeH = 0 # 2*l+1
        self.n_inde = 0 # nb of independent channel
        self.sigind = 0 # tab of index from compressed to matrix form [nsp,m1,m2]
        self.offdiag = False # true if off diag m1,m2 present
        self.nsp = 0 # nb of indepent spin
        self.delta = 0 # [nsp,chan,iw]

        self.g0 = 0 # [nsp,m1_active,m2_active,iw]
        self.G0 = 0    # [nsp,chan,iw]

        self.sig = 0   # [nsp,chan,iw]
        self.gl = []    # [nsp,chan,iw]

        self.Eimp = 0  # [nsp,chan]

        self.lst_active = [] # list of m witg sigind_mm non zeros
        self.n_active = 0 # nb of active m
        self.W=[]
        self.beta = 0
        self.nomg = 0

        self.lst_chan=0

    def readh5(self,DC=True) :
        """
        read necessary info from solver_input.h5
        see doc on this h5 file
        """
        f = h5py.File('solver_input.h5', 'r')

        self.beta = f['beta'][0]
        self.nomg = int(f['nomg'][0])
        self.nsp = int(f['nsp'][0])
        self.W = np.array([ (2*n+1)*1j*np.pi/self.beta for n in range(self.nomg)])

        self.L = f['l'][self.imp_ind]
        self.sizeH = 2*self.L + 1

        self.sigind = np.array(f['sigind'][:,self.imp_ind,:self.sizeH,:self.sizeH])

        self.sigind_red=[] # same as sigind but just for active chan

        # get active m
        self.active_m()

        self.lst_chan=self.getchan() # lst[chan] is a list of (m1,m2) where sigind_red[m1,m2]=chan. (usefull to iterate over chan.

        self.offdiag = not(isdiag(self.sigind[0]))

        self.delta = np.array(f['delta'])
        if DC :
            self.Eimp = np.array(f['eimp'])
        else :
            self.Eimp = np.array(f['eimp_wo_dc'])
        self.g0=self.getG0() # in mat form
        self.sig=0*self.g0 # in mat form


    def active_m(self) :
        """
        return non zeros diag index
        """

        for m in range(self.sizeH) :
            if self.sigind[0,m,m] > 0 :
                self.lst_active.append(m)
        self.n_active = len(self.lst_active)

        self.sigind_red=np.zeros((self.nsp,self.n_active,self.n_active), dtype=int)
        for (im1,m1), (im2,m2) in product(enumerate(self.lst_active), enumerate(self.lst_active)) :
            self.sigind_red[:,im1,im2]=self.sigind[:,m1,m2]


    def matrixform(self,tab) :
        """
        if tab is rank 1 : no omega
                       2 : depends of omega
        """
        if len(tab.shape) == 1  :
            omega=False
        elif len(tab.shape) == 2  :
            omega=True
        else :
            raise Exception(' matrixform : size of tab not managed')
        if omega :
            M=np.zeros((self.nsp, self.n_active, self.n_active,self.nomg), dtype=tab.dtype)
        else :
            M=np.zeros((self.nsp, self.n_active, self.n_active), dtype=tab.dtype)



        for (im1,m1), (im2,m2) in product(enumerate(self.lst_active), enumerate(self.lst_active)) :
            for spin in range(self.nsp) :
                chan = self.sigind[spin, m1, m2] - 1 # fortran convention !!!!!
                if chan < 0 : continue # channel inactive
                if omega :
                    M[spin, im1, im2,:] = tab[chan, :]
                else :
                    M[spin, im1, im2] = tab[chan]

        return M

    def getG0(self) :
        """
        G0= (iw_n-delta-Eimp)^-1
        """
        Emat=self.matrixform(self.Eimp)
        deltamat=self.matrixform(self.delta)

        g=np.zeros((self.nsp, self.n_active, self.n_active, self.nomg), dtype=complex)
        de=np.zeros((self.nsp, self.n_active, self.n_active, self.nomg), dtype=complex) # delta+ E

        Id=np.identity(self.n_active, dtype=complex)
        for spin in range(self.nsp) :
            for iw, w in enumerate(self.W) :
                g[spin,:,:,iw]= np.linalg.inv( self.W[iw] * Id - deltamat[spin,:,:,iw] - Emat[spin])
        return g


    def getchan(self) :
        """
        read sigind_red to get  :         lst[chan-1]  a list of index whith the channel chan ( the -1 is from fortran to python convention....)
        example if sigred =                   1 0 0 0
                                              0 1 0 0
                                              0 0 2 5
                                              0 0 4 3
        it returns [ [(0,0),(1,1)], [(2,2)], [(3,3)], [(3,2)], [(2,3)]]
        """
        mx=0

        nb_ind_chan=max([ a for a in self.sigind_red[0].reshape(self.n_active**2) if a>0]) - min([ a for a in self.sigind_red[0].reshape(self.n_active**2) if a>0]) +1

        lst=[ [] for i in range(nb_ind_chan)]
        shift=min([ a for a in self.sigind_red[0].reshape(self.n_active**2) if a>0])
        for i, j in product(range(self.n_active),range(self.n_active)) :
            if  self.sigind_red[0,i,j] > 0 :
                lst[self.sigind_red[0,i,j]-shift].append((i,j))
        return lst

    def compress_form(self,M) :
        """
        first index of M are M[nsp,m1,m2]
        return lst[nsp,chan]= average of M[nsp,sigind[m1,m2]] for sigind[m1,m2]=chan
        if nsp = 1 , [lst]
        if nsp = 2 , return [lstu,lstd]
        """
        lst=[]
        if self.nsp > 1 :
            lstd=[]
        for chan, lst_indices in enumerate(self.lst_chan) :
            mean=np.mean([ M[0, ind[0], ind[1]] for ind in lst_indices],axis=0)
            lst.append(mean)

            if self.nsp > 1 :
                mean=np.mean([ M[1, ind[0], ind[1]] for ind in lst_indices],axis=0)
                lstd.append(mean)

        if self.nsp == 1 :
            return [lst]
        else :
            return [lst,lstd]



def isdiag(M) :
    if  np.count_nonzero(M - np.diag(np.diagonal(M))) > 0 :
        return False
    else :
        return True
