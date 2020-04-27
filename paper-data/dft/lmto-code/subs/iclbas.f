      integer function iclbas(class,ipc,nbas)
C- Returns an index to iclbas atom in basis given class
C ----------------------------------------------------------------
Ci Inputs
Ci   class :class for which to find basis atom
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nbas  :size of basis
Co Outputs
Co   iclbas:index to first occurrence of site ib in ipc table
Co         :returns 0 if class has no corresponding site
Cu Updates
Cu   08 Nov 12  Size of ipc is now passed
C ----------------------------------------------------------------
      implicit none
C Passed parameters
      integer nbas,class,ipc(nbas)
C Local parameters
      integer ibas

      do  ibas = 1, nbas
        iclbas = ibas
        if (ipc(ibas) == class) return
      enddo
      iclbas = 0
      end
      integer function iclbsj(ic,ipc,nbas,nrbas)
C- Returns an index to nrbas atom in basis given the class
C ----------------------------------------------------------------------
Ci Inputs:
Ci   ic    :class index
Ci   ipc   :class index: site ib belongs to class ipc(ib) (mksym.f)
Ci   nbas  : abs  = number of atoms in the basis
Ci         : sign <0 to return with -n if there are fewer than nrbas
Ci         : members of class ic, where n=number members of class ic
Ci   nrbas :(nrbas>0) Find the nrbas-th occurrence of class ic
Ci         :(nrbas<0) Count the number of sites belonging to class ic
Co Outputs:
Co   iclbsj:(nrbas>0) index to the nrbas-th atom belonging to class ic
Co         :(nrbas<0) the number of sites belonging to class ic
Cb Bugs
Cb   iclbsj always returns 1 if nrbas=0 ---
Cb          routine does not confirm that nrbas is nonzero.
Cu Updates
Cu   08 Apr 11 Added nrbas<0 option
Cu   10 May 01 Returns -n rather than -1 when nrbas-th atom not found
C ----------------------------------------------------------------------
      implicit none
C Passed parameters
      integer ic,nbas,ipc(abs(nbas)),nrbas
C Local parameters
      integer ibas,n,nbasa

      iclbsj = 1
      n = 0
      nbasa = iabs(nbas)
      do  ibas = 1, nbasa
        if (ipc(ibas) == ic) n = n+1
        if (n == nrbas) then
          iclbsj = ibas
          return
        endif
      enddo

      if (nrbas < 0) then
        iclbsj = n
        return
      endif

      if (nbas < 0) then
        iclbsj = -n
        return
      endif

      call fexit3(-1,111,' Exit -1 ICLBSJ: sought atom no. %i'//
     .  ' in class %i but only %i atoms exist',nrbas,ic,n)

      end
