C --- MvS's patches to MSM mfdef,mopen,mclose ---
      subroutine mfdef(ipr)
C- General initialization, pick up file extension
C (MSM version: file extension from environment variable MOLID)
      call finits(2,0,0,i)
      end
      subroutine mopen(ifi,ext,type)
C- File opening, using ext for fnam
C  For now, map the following files: 'lg' => 'log'
      implicit none
      integer ifi
      character ext*(*), type*1
      character fnam*4
      integer fopn,fadd,fhndl,jfi,mode,iprint,i1mach

      fnam = ext
      if (fnam == 'lg') fnam = 'log'
      mode = 0
      if (type == 'u' .or. type == 'U') mode = 4

      jfi = fhndl(fnam)
C      if (jfi == -1) jfi = fadd(fnam,ifi,mode)-1
      if (jfi == -1) jfi = fadd(fnam,ifi,mode)
      if (jfi == ifi) then
        jfi = fopn(fnam)
      else
        call fexit2(-1,111,' Exit -1 mopen: '//
     .  'sought to assign unit %i to '''//fnam//'%a'' when assigned %i',
     .    ifi,jfi)
      endif
      if (iprint() > 100) call awrit3(' mopen: fnam='//fnam//
     .                       ' ifi=%i jfi=%i mode=%i',' ',256,i1mach(2),
     .                       ifi,jfi,mode)
      end
      subroutine mclose(ifi)
      call fclose(ifi)
      end
C      subroutine fmain
C      implicit none
C      call mfdef(1)
C      call mopen(82,'pd','f')
C      call mopen(82,'pd','f')
C      call mopen(71,'lg','f')
C      call mclose(82)
C      call mopen(83,'p3','f')
C      call mopen(84,'pd','u')
C      end
