      subroutine newalp(ldim,nl,nbas,ipc,nsp,isp,pp,alpha)
C- Interactive input of new screening constants
      implicit none
      integer nl,nbas,ldim,ldmax,nsp,ipc(*)
      parameter (ldmax = 360)
      double precision alpha(0:nl**2-1,0:nbas),pp(6,0:nl-1,nsp,1)
C Static
      double precision copy(ldmax)
      integer i,iprint,lm,ibas,ll,isp,stdo,nglob
      external iprint
      double precision const
      logical done
      save copy, done
      data copy /ldmax*0d0/ done /.false./


      stdo = nglob('stdo')
      if (ldim > ldmax) stop 'NEWALP : Increase static dimension'
      call dcopy(ldim,copy,1,alpha,1)
      if (done) return
      print *, 'NEWALP : set all alpha''s to a constant ? 1 for yes'
      read *, i
      if (i == 1) then
        print *, 'Enter the constant'
        read *, const
        call dcopy(ldim,const,0,alpha,1)
        call dcopy(ldim,alpha,1,copy,1)
        done = .true.
        if (iprint() > 30) write (stdo,10) const
        return
      endif
      print *, 'NEWALP : change to ''gamma'' representation ? 1 for yes'
      read *, i
      if (i == 1) then
        do  ibas = 1, nbas
          do  lm = 0, nl**2-1
          alpha(lm,ibas) = pp(5,ll(lm+1),isp,ipc(ibas))
          enddo
        enddo
        call dcopy(ldim,alpha,1,copy,1)
        done = .true.
        return
      endif
      print *, 'Input the vector alpha'
      read *, alpha
      call dcopy(ldim,alpha,1,copy,1)
      done = .true.
      if (iprint() > 30) write (stdo,20)
      if (iprint() > 30) write (stdo,30)
     .                        ((alpha(lm,ibas),lm=0,nl-1),ibas=1,nbas)
   10 format(8x,' All screening constants set to ',f10.6)
   20 format(8x,' New screening constants :')
   30 format(9f8.3)
      end
