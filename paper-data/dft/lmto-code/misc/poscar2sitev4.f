      subroutine fmain
      implicit none
      character strn*100,direct*10
      integer ifi,jfi,fopng,awrite,i,j,nbas,ibas(20),it(20),a2vec,nspec
      double precision vol,plat(3,3),vol0,avwsr,a,alat,pos(3)

      jfi = fopng('POSCAR',-1,1)
      read(jfi,*)
      read(jfi,*) vol
      vol = -vol
      read(jfi,*) plat(:,1)
      read(jfi,*) plat(:,2)
      read(jfi,*) plat(:,3)

      a = avwsr(plat,1d0,vol0,1)
C     Convert alat to a.u.
      alat = (vol/vol0)**(1d0/3d0) / 0.529177d0

      read(jfi,'(100a)') strn
      i = 0
      ibas = -1
      nspec = a2vec(strn,len(strn),i,2,', ',2,-3,9,it,ibas)
      if (nspec <= 0) call rx('poscar2site: no species')
      read(jfi,'(100a)') direct

      nbas = 0
      do  i = 1, nspec
        do  j = 1, ibas(i)
          nbas = nbas + 1
          read(jfi,*) pos,strn
          print 333, nbas,pos,trim(strn)
  333     format(i4,3f14.8,2x,a)
        enddo
      enddo

      rewind jfi
      read(jfi,*)
      read(jfi,*)
      read(jfi,*)
      read(jfi,*)
      read(jfi,*)
      read(jfi,'(100a)') strn
      read(jfi,'(100a)') direct


      ifi = fopng('site',-1,0)
      rewind ifi
      i = awrite('%% site-data vn=3.0 xpos fast io=15 '//
     .  'nbas=%i alat=%d plat=%9:1;8,1d',' ',180,ifi,
     .  nbas,alat,plat,a,a,a,a,a)

      do  i = 1, nspec
        do  j = 1, ibas(i)
          read(jfi,*) pos,strn
          write(ifi,321) strn,pos
  321     format(1x,a8,3f17.10)
        enddo
      enddo




      end
