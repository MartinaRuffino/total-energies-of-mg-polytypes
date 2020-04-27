      subroutine rdQGcou(mode,nqc,npwmbx,QpGcut_cou,qp,npwmb,ngvecc)
C- Read from file QGcou
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :1s digit
Ci         :0 skip reading header.
Ci         :  (header must have been already read)
Ci         :1 rewind file and read header nqc,npwmbx,QpGcut_cou
Ci         :10s digit
Ci         :0 => return without reading body
Ci         :1 for each q in 1:nq read qp,npwmb,ngvecc
Ci         :  If nq is 1, only data for the next point is read
Ci         :  In this case the file pointer should not be rewound
Ci         :  between calls
Ci  QpGcut_cou
Cio Inputs/Output
Cio  nqc   :number of qp at which coulomb interaction is calculated
Cio        :(called nqnumc in old GW code)
Cio        :Output if 1s digit mode = 0, otherwise input
Cio  npwmbx:max number of G vectors for Coulomb interaction and
Cio        :max number of interstitial product functions for any q
Cio        :(called ngcmx in old GW code)
Cio        :Output if 1s digit mode = 0, otherwise input
Co Outputs
Co  qp     :vector of qp, for each of nqc points
Co  npwmb  :number of G vectors for Coulomb interaction and
Co         :number of interstitial product functions (q dependent)
Co         :(called ngc in old GW code)
Co  ngvecc :G-vectors for coulomb interaction, as multiples of qlat
Cr Remarks
Cr
Cu Updates
Cu   22 Dec 17 Adapted from readqg.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nqc,npwmbx
      double precision QpGcut_cou
      integer :: ngvecc(3,npwmbx,nqc),npwmb(nqc)
      real(8) :: qp(3,nqc)
C ... Local parameters
      integer ifi,iq
      procedure(integer) :: fopnx

      ifi = fopnx('QGcou',2,4,-1)
      if (mod(mode,10) > 0) then
        rewind(ifi)
        read(ifi) nqc, npwmbx, QpGcut_cou
      endif
      if (mode < 10) return

      do  iq = 1, nqc
        read(ifi) qp(1:3,iq), npwmb(iq)
        read (ifi) ngvecc(1:3,1:npwmb(iq),iq)
      enddo
      end

      subroutine rdQIBZ(mode,nqibz,nq0i,qibz,wibz)
C- Read qp from files QIBZ and Q0P
C ----------------------------------------------------------------------
Ci Inputs
Ci  mode   :1s digit
Ci         :0 skip reading header.
Ci         :  (header must have been already read)
Ci         :1 rewind file QIBZ and read nqibz from header
Ci         :2 rewind file Q0P  and read nq0i from header
Ci         :3 combination of 1+2
Ci         :10s digit
Ci         :0 => return without reading body
Ci         :1 for each q in 1:nqibz read qp,wibz from file QIBZ
Ci         :1 If nq0i > 0 from file Q0P
Co Outputs
Co  qibz   :vector of qp, for each of nqibz+nq0i points
Co  wibz   :weights corresponding to q
Cr Remarks
Cu Updates
Cu   22 Dec 17 Adapted from readqg.F
C ----------------------------------------------------------------------
      implicit none
C ... Passed parameters
      integer mode,nqibz,nq0i
      real(8) :: qibz(3,nqibz+nq0i),wibz(nqibz+nq0i)
C ... Local parameters
      double precision xv(4)
      integer ifi,iq
      procedure(integer) :: fopnx

      if (mod(mod(mode,10),2) > 0) then
        ifi = fopnx('QIBZ',2,1,-1)
        rewind ifi
        read(ifi,*) nqibz
      endif

      if (mod(mode,10)/2 > 0) then
        ifi = fopnx('Q0P',2,1,-1)
        rewind ifi
        read(ifi,*) nq0i
      endif

      if (mode < 10) return

      ifi = fopnx('QIBZ',2,1,-1)
      do  iq = 1, nqibz
        read(ifi,*) xv
        qibz(1:3,iq) = xv(1:3)
        wibz(iq) = xv(4)
      enddo

      if (nq0i == 0) return

      ifi = fopnx('Q0P',2,1,-1)
      do  iq = 1, nq0i
        read(ifi,*) xv  ! weight comes before qp in this file!
        qibz(1:3,nqibz+iq) = xv(2:4)
        wibz(nqibz+iq) = xv(1)
      enddo

      end
