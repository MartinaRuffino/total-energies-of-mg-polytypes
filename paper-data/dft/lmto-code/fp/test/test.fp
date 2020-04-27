#!/bin/tcsh -f

# A shell script testing operation of fp suite
# set verbose

alias call 'set retcall = \!\!:2 ; set callarg = \!\!:3 ; goto \!\!:1'
alias runjob 'set retcall = \!\!:1; set outfile = \!\!:2 ; set callarg = \!\!:3 ; goto runjob'
alias runrdcmd 'set retcall = \!\!:1; set rdcmdfmt = \!\!:2 ; set outfile = \!\!:3 ; set callarg = \!\!:4 ; goto runrdcmd'
alias findcmd  'set retcall = \!\!:1 ; set prog_cmd = \!\!:2 ; set path_name = \!\!:3 ; set make_path = \!\!:4 ; goto findcmd'
alias extract_res_n 'set retcall = \!\!:1; set testvar = \!\!:2 ; set refvar = \!\!:3 ; set keyword = \!\!:4  ; set arg_number = \!\!:5 ; set occur_number = \!\!:6 ; set sed_strn = \!\!:7 ; goto extract_res_n'
alias compare_res 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set refvar = \!\!:4 ; set tol = \!\!:5 ; set passvar = \!\!:6 ; goto compare_res'
alias compare_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set tol = \!\!:4 ; set passvar = \!\!:5 ; goto compare_res_0'
alias compare_resf 'set retcall = \!\!:1; set testvar = \!\!:2 ; set refvar = \!\!:3 ; set keyword = \!\!:4  ; set arg_number = \!\!:5 ; set occur_number = \!\!:6 ; set sed_strn = \!\!:7 ; goto compare_resf'
#alias zcmpmfiles_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set tol = \!\!:4 ; set passvar = \!\!:5 ; set ndig = \!\!:6 ; set srcfile = \!\!:7 ; set reffile = \!\!:8 ; goto zcmpmfiles_res_0 '
alias cnvt_d_fmt  'set retcall = \!\!:1; set testvar = \!\!:2 ; set testval = \!\!:3 ; goto cnvt_d_fmt'
alias query 'set retcall = \!\!:1 ; set retcall2 = \!\!:2 ; set callarg = \!\!:3 ; goto query'
alias zcmpmfiles_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set tol = \!\!:3 ; set passvar = \!\!:4 ; set ndig = \!\!:5 ; set srcfile = \!\!:6 ; set reffile = \!\!:7 ; goto zcmpmfiles_res_0 '
alias zcmpmfiles_res_tol 'set retcall = \!\!:1; set keyword = \!\!:2 ; set tol = \!\!:3 ; set passvar = \!\!:4 ; set ndig = \!\!:5 ; set srcfile = \!\!:6 ; set reffile = \!\!:7 set nlines = \!\!:8; goto zcmpmfiles_res_tol '
alias zcmpmcx 'set retcall = \!\!:1; set keyword = \!\!:2 ; set toldif = \!\!:3 ; set ndiff = \!\!:4 ; set prec = \!\!:5 ; set passvar = \!\!:6 ;set srcfile = \!\!:7 ; set reffile = \!\!:8 ; goto zcmpmcx '

set allargs = ($argv)

set a
set slow
set testfile = $0
set testdir = $testfile:h
set topdir  = `cd $testdir/../..; pwd`
set tmpdir = $cwd
set space = '        '
set jobspassed
set failed = 0
alias zcat 'gunzip -c'
alias zcat 'cat'

#alias mpix mpirun

# Prepend current working-directory, top-level dir to path
set path = ($cwd $topdir $topdir/utils $topdir/testing $path)

set plot = `which fplot`
if (-x "$plot") then
  if `$plot -h | sed -n 1,1p | awk '{print ($1 == "fplot")}'` set have_fplot
endif
set mcx = `which mcx`
if (-x "$mcx") then
  if `$mcx --h |& sed -n 1,1p | awk '{print ($7 == "(vsn" && ($8 * 1 >= 1.04))}'` set have_mc
endif
set pldos = `which pldos`
if (-x "$pldos") then
  if `$pldos -h | sed -n 1,1p | awk '{print ($2 == "pldos")}'` set have_pldos
endif
# see if ghostscript is available
set gs = `which gs`
if (-x "$gs") then
  if `$gs --help | sed -n 1,1p | awk '{print ($2 == "Ghostscript")}'` set have_ghostscript
endif

# see if gnu grep is available
echo X | grep -A 1 X > & /dev/null
set retval = $status
if ($retval == 0) set gnu_grep

# --- Pick off switches ---
while (`echo $1 | sed -e 's/\(.\).*/\1/' `  ==  "-")

  set arg1 = $1; shift
  if ($?verb) echo test.fp: parsing switch $arg1
  switch ($arg1)
    case "--quiet":
      set quiet
      unset slow
      breaksw
    case "--add0":
      set add0 = `which add0`
      if (! -x "$add0") then
        echo "test.dmft (abort): missing add0"
        exit
      endif
      breaksw
    case "--poszer":
      set poszer = `which poszer`
      if (! -x "$poszer") then
        echo "test.dmft (abort): missing poszer"
        exit
      endif
      breaksw
    case "--MPIK":
      set MPI
      set MPIK
      unset MPI
      breaksw
    case "--MPI":
      set MPI
      set MPIK
      unset MPIK
      breaksw
    case "--clean":
      set clean
      breaksw
    case "--veryclean":
      set clean
      set veryclean
      breaksw
    case "--no-iact*":
      unset slow
      breaksw
    case "--haveout":
      set haveout
      breaksw
    case "--list":
      goto showtests
      breaksw
    case "--whichexec"
      set quiet; unset quiet
      findcmd chk00 lmf "$path" "$topdir"
      chk00:
#       Doesn't work
#       if ($?MPIK) then
#       findcmd chk01 lmf "$path" "$topdir"
#       chk01:
#       endif
      exit 0
      breaksw
    case "--libxc":
      set libxc
      breaksw
    case "--noplot*":
      set noplot
      set have_pldos
      unset have_pldos
      set have_fplot
      unset have_fplot
      breaksw
    case "--verb*":
      set verb = 1
      breaksw
    case "--so4":
      set so4
      breaksw

    case "--all":
      set mater_lst = (copt te zrt co cr3si6 fe cu au gas srtio3 al cdte felz gasls zbgan gdn eras c al fept mgo crn cs)
      set joblist
      while (`echo $1 | sed -e 's/\([0-9][0-9]*\)/-/'`  ==  "-")
        set joblist = ($joblist $1)
        shift
      end
      set pass
      set failed
      foreach i ($mater_lst)
        $testfile `echo $allargs | sed s/--all//g | sed -e 's/\([0-9][0-9]*\)//g' | sed -e 's/-add/-add0/g'` $i $joblist
        set retval = $status
        if ($retval != 0) then
          unset pass
          set failed = ($failed $i)
#         echo " $testfile : failed test $i ... aborting"
#            exit -1
        endif
      end
      if ($?clean) then
        exit
      else if ($?pass) then
        echo "$space $testfile : all tests PASSED ($mater_lst)"
        exit
      else
        echo "$space $testfile : These tests FAILED:" $failed
        exit -1
      endif

    default:
      echo unrecognized switch $arg1
      goto usage
  endsw

end

echo ' '
echo "         ---- test.fp: test FP program lmf ---"

# --- use copt as default in the absence of specific choice ---
if ($#argv == 0) then
  set ext = copt
  echo "$space .... no file extension specified; use input file ctrl.copt"
else
  set ext = $argv[1]
  shift
endif

if (! -e $testdir/ctrl.$ext) then
   echo ' '
   echo " test.fp aborting ... missing file $testdir/ctrl.$ext"
   goto usage
endif

if ($ext == "copt") then
  echo '         Case CoPt: a distorted L12 environment with four atoms.'
#    echo '         Other checks:'
#    echo '         spin-pol, tetrahedron+metal=3, forces(mode 12),'
#    echo '         2-kappa basis, inequiv kmxa,lmax,rmt'
  set cplst = ($testdir/{ctrl.copt,spec.prop})
# set fitbas  = ' Co RSMH= 2.425 2.717 1.017 EH= -0.360 -0.200 -0.222'
# set fitbas2 = ' Pt RSMH= 2.461 3.042 1.085 EH= -0.441 -0.200 -0.200'
  set dfmax1tol1 = 0.1
  set dfmaxntol1 = 0.1
else if ($ext == "bzt") then
  echo '         Case Ba3ZnTa2O9: an oxide.  Tests charge density plotting mode'
  set cplst = ($testdir/{ctrl.bzt,site.bzt})
else if ($ext == "coptso") then
  echo '         Case coptso: a distorted L12 environment with four atoms.'
#    echo '         Other checks:'
#    echo '         Spin-orbit coupling, and scaling of SO,'
  set cplst = ($testdir/{ctrl.coptso,socscl.coptso,spec.prop})
  set dfmax1tol1 = 0.1
  set dfmaxntol1 = 0.1
else if ($ext == "fe") then
  echo '         Case Fe: spin-polarized Fe in bcc structure'
  set cplst = ($testdir/{ctrl.fe,fs.fe,shfac.fe})
  set dosmulltol = 7e-3
  set dosclstol = .0015
  set jdostol = .0008

else if ($ext == "fept") then
  set cplst = ($testdir/atparms $testdir/fept/{ctrl,site,site2,site3}.fept)
# lmfa ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2
# lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-6 -vbeta=.3 -vmet=5 --iactiv
# make FM supercell and show it is approximately self-consistent
# lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 '--rsedit~rs~scell site2'
# lmf ctrl.fept -vnsp=2 -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 '--rsedit~rs~scell site2'; cp rst2.fept rst.fept
# lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-4 -vbeta=.3 -vmet=5 --rs=1,2 -vnk2=nk/2 -vafm=0 -vfile=2 --iactiv
# spin flip reproduces same result
# lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-4 -vbeta=.3 -vmet=5 --rs=1,2 -vnk2=nk/2 -vafm=0 -vfile=2 '--rsedit~rs~set all 1 flip~savea~a'
# lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-4 -vbeta=.3 -vmet=5 --rs=2,0 -vnk2=nk/2 -vafm=0 -vfile=2 --iactiv
# next 2 lines demonstrate that FA data is completely overwritten with 'add all 0 1'
# lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-4 -vbeta=.3 -vmet=5 --rs=1,2 -vnk2=nk/2 -vafm=0 -vfile=2 '--rsedit~rsfa~read'
# lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-4 -vbeta=.3 -vmet=5 --rs=2,0 -vnk2=nk/2 -vafm=0 -vfile=2 --iactiv
# next 3 lines do the same with afm on (suppress AFM:: in ctrl file, remove specialspec1)
# lmfa ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-4 -vbeta=.3 -vmet=5 -vnk2=nk/2 -vafm=1 -vfile=3
# rsa
# reada
# set all 1 zers
# set asite 2 flip 2
# add asite 1 1 1,2 1,2
# save
#  lmf ctrl.fept -vtpp1=3 -vtpp2=3 -vtpd1=4 -vtpd2=4 -vbigbas=2 -vehmax=-.45 -vnk=10 -vnsp=2 -vconvc=1e-4 -vbeta=.3 -vmet=5 --rs=1,0 -vnk2=nk/2 -vafm=1 -vfile=3 --iactiv
else if ($ext == "kfese") then
  echo '         Case KFeSe: striped AFM phase'
  set cplst = ($testdir/{ctrl.kfese,shfac.kfese,semi.mater})
else if ($ext == "felz") then
  echo '         Case felz: spin-polarized Fe spin-orbit coupling'
  set cplst = ($testdir/lzsz/{ctrl,rsta}.felz)
  set dorbmtol = 0.00001
else if ($ext == "ni") then
  echo '         Case Ni: fixed-spin-moment method'
  set cplst = ($testdir/{ctrl.ni,site.ni,atparms})
else if ($ext == "gasls") then
  echo '         Case GaAs: GaAs with spin-orbit coupling'
  set cplst = ($testdir/ls/{ctrl,rsta,syml}.gasls)
  set gmtol = 0.0001
# lineeval specifies which line evals of interest are in AFTER bndfp: tag
# eval1 eval2 specifies ranges of columns evals are in
  set lineeval=3  eval1=4  eval2=9  evalso=5
  set lineevald=2  eval1d=7  eval2d=9  evalsod=1
else if ($ext == "gaslc") then
  echo '         Case GaAs: GaAs with spin-orbit coupling, local orbitals and QSGW self-energy'
  set cplst = ($testdir/{ls/ctrl.gaslc,ls/rst.gaslc,ls/syml.gaslc,ls/sigm.gaslc,semi.mater})
  set gmtol = 0.0001
# lineeval specifies which line evals of interest are in AFTER bndfp: tag
# eval1 eval2 specifies ranges of columns evals are in
  set lineeval=4  eval1=5  eval2=9  evalso=6
else if ($ext == "cs") then
  set cplst = ($testdir/{ctrl,basp,site}.$ext)
  set gmtol = 0.0001
  set dehf1tol1 = 1e-5
else if ($ext == "gdn") then
#    echo '         Case GdN: Test of LDA+U, and also LDA with spin polarized 4f core'
  set cplst = ($testdir/{ctrl.gdn,occnum.gdn,bfield.gdn,syml.gdn})
  set gmtol = 0.0001
  set dosmulltol = 2e-4
else if ($ext == "cdte") then
  echo '         Case CdTe: Test of LDA+U in different modes'
  set cplst = ($testdir/{ctrl.cdte,occnum.cdte,semi.mater})
  set gmtol = 0.0001
else if ($ext == "eras") then
  echo '         Case ErAs: Test of LDA+U'
  set cplst = ($testdir/{ctrl.eras,occnum.eras})
  set gmtol = 0.0001
else if ($ext == "er") then
  echo '         Case Er: Test of LDA+U'
  set cplst = ($testdir/{ctrl.er,site.er,occnum.er,syml.er,specialspec1,atparms})
  set gmtol = 0.0001
  set dehf1toln = 6e-6
else if ($ext == "zrt") then
  echo '         Case zrt:'
  set cplst = ($testdir/{ctrl.zrt})
  set dfmax1tol1 = 0.01
  set dehf1toln = 6e-6
else if ($ext == "co") then
  echo '         Case Co: a hexagonal environment with two equivalent atoms.'
#    echo '         Other checks:'
#    echo '         spin-polarized, tetrahedron+metal=4'
  set cplst = ($testdir/{ctrl,syml,qp}.co)
  set lmfargs1 = "-vmet=4 -vlmf=1 -vnk=8 -vnit=10 --pr31"
# set fitbas = ' A RSMH= 2.375 2.722 1.047 EH= -0.342 -0.200 -0.236'
  set drmsqtol1  = 1e-6
  set pdostol    = 0.01
  set dostol = .007
else if ($ext == "cr3si6") then
  echo '         Case Cr3Si6: a hexagonal environment with several atoms, two classes'
  set cplst = ($testdir/{ctrl.cr3si6})
  set lmfargs1 = "--pr51 -vnit=2 --no-iactiv --time=6"
# set fitbas  = ' Cr RSMH= 2.959 3.303 1.275 EH= -0.307 -0.200 -0.200'
# set fitbas2 = ' Si RSMH= 1.729 1.727 -1.000 EH= -0.609 -0.200 0.000'
  set pdostol = 4e-4
else if ($ext == "te") then
  echo '         Case Te: molecular statics in an open structure'
#    echo '         There is only one symmetry-allowed degree of freedom.'
  set cplst = ($testdir/{ctrl.te})
  set lmfargs2 = '-vminx=t --rs=1,1 -vnk=3 -vnit=20 -vlf1=4 -vlmxl=4 -vnk=3 -vngd=20 -vkmx=3 -vconv=1d-4 -verefc=0 -verefa=0'
  set dfmaxntol1 = 0.002
  set drmsqtol1  = 1e-6
else if ($ext == "srtio3") then
  echo '         Case SrTiO3: an oxide with local orbitals.'
  set cplst = ($testdir/{ctrl.srtio3})
  set lmfargs1 = " "
  set dfmax1tol1 = 1e-5
  set drmsqtol1 = 2e-6
else if ($ext == "al") then
  echo '         Case Al: simple fcc metal'
  set cplst = ($testdir/{ctrl.al})
  set lmfargs1 = " "
  set dfmax1tol1 = 1e-5
  set drmsqtol1 = 2e-6
else if ($ext == "tio2") then
  echo '         Case tio2: example of relaxation'
  set cplst = ($testdir/{ctrl,site}.tio2)
  set lmfargs1 = " "
  set dfmax1tol1 = 1e-5
  set drmsqtol1 = 2e-6
else if ($ext == "cu") then
  echo '         Case Cu: illustration of high-lying local orbitals'
  echo '                  and bands of Cu up to ~50 eV.'
  set cplst = ($testdir/{ctrl.cu,syml.cu})
  set lmfargs1 = " "
  set drmsqtol1 = 1e-6
else if ($ext == "na") then
  echo '         Case na: illustration of low- and high-lying local orbitals'
  set cplst = ($testdir/{ctrl.na})
  set lmfargs1 = " "
  set drmsqtol1 = 1e-6
else if ($ext == "c") then
  echo '         Case C: test of homogeneous background'
  set cplst = ($testdir/{ctrl.c})
  set lmfargs1 = " "
  set drmsqtol3 = 5e-6
else if ($ext == "crn") then
  echo '         Case CrN: test of CLS with core hole'
  set cplst = ($testdir/{ctrl.crn})
  set lmfargs1 = " "
  set drmsqtol1 = 1e-6
else if ($ext == "gas") then
  echo '         Case GaAs: GaAs with local 3d orbital.'
  set cplst = ($testdir/{ctrl.gas})
else if ($ext == "zbgan") then
  echo '         Case zincblende GaN: calculate joint DOS and Im eps'
  set cplst = ($testdir/{ctrl.zbgan})
else if ($ext == "au") then
  echo '         Case Au: Set up for high accuracy calculation of the equilibrium volume'
  set cplst = ($testdir/{ctrl,site}.au)
else if ($ext == "mgo") then
  echo '         Case MgO: test of external potential'
  set cplst = ($testdir/{ctrl,basp,site}.mgo)
else
  echo test.fp: No test case for $ext
  exit -1
endif
endif

if ( $?joblist == 0 ) then
set joblist = ($argv)
if ( $#joblist == 0 ) set joblist = (1 2 3 4)
endif

if ($ext == "fept") goto chk2e

echo $joblist | egrep '\b1\b' >/dev/null
if ($status) goto chk1e
set jobid = 1

cat <<EOF

         --- Test 1.  Basic check of programs lmfa,lmf ---
         Checks that program lmfa produces a sensible atm file
         and that program lmf iterates to the proper energy.

EOF
endif

if (0) then
  echo oops ..
  exit
else if ($ext == "fe") then
if ($?quiet) then
cat <<EOF
         The Fe test illustrates:

         1. how to generate a 2D Fermi surface

         2. how to generate k-resolved spin 1 DOS

         3. the optics editor --popted to extract DOS(k) near the Fermi level

         4. Im eps calculated with linear optics (no local fields)

EOF
else
cat <<EOF
       * The Fe test illustrates how to invoke lmf in order to
         generate a 2D Fermi surface (minority spin).

         After lmf completes, the following steps will draw a Fermi surface:
         (you need the FPLOT package, and the mcx calculator for these steps)
           mcx -r:open bnds.fe -shft=0 -w b2 -r:open bnds.fe -shft=0 -w b3 -r:open bnds.fe -shft=0 -w b4 -r:open bnds.fe -shft=0 -w b5
           fplot -f fp/test/plot.fs0
           gv fplot.ps

       * This test also generates Im eps(omega).
         The collinear case is stored in opt.fe
         Case SO coupling (SO=3) is stored in opt-so3.fe
         Case SO coupling (SO=1) is stored in opt-so1.fe
         To compare the three, try:
           fplot -lt 1 -ord x2+x5 opt-so3.fe -lt 2,col=1,0,0 opt-so1.fe -lt 3,col=0,1,1 -ord x2+x5 opt.fe

       * You can test the PBE functional using --libxc
         You can compare the result against Kieron Burke's easypbe (internally compiled),
         uncomment the line .TESTXC in fp/test/ctrl.fe :
         # PBE using easypbe functional : XCFUN={lxcf} GGA={lxcg}
         .TESTXC  lmfa -vlxcf=4 -vlxcg=3 fe
         You should find that the results are nearly identical.

       * This test also generates k-resolved spin 1 DOS, using the joint DOS code optinq.f.
         The total DOS can be compared against the standard DOS generated with SAVDOS=T.
         This example has 140 irreducible k-points.

         It also projects the DOS onto the Fe t2g and eg orbitals (through switches --jdosw=5,6,8,21,22,24 --jdosw2=7,9,23,25)
         Thus:
           jdos.fe  contains both total DOS (col 2) and the t2g and eg Mulliken projections (cols 3 and 4)
           poptb.fe contains k-resolved DOS (cols 2:141), the t2g and eg projections of it (cols 142:281 and 282:421)

       * To compare the optinq-generated DOS to the standard SAVDOS=T do:
           echo 100 40 -9 9 |  pldos -escl=13.6 -fplot -lst='1' -lst2 dos.fe
           fplot -disp -frme 0,1,0,.7 -y 0,30 -p0 -ab 'x1*13.6' jdos.fe -ab 'x1*13.6' \
             -colsy 3 -lt 2,col=0,1,0 jdos.fe -ab 'x1*13.6' -colsy 4 -lt 2,col=0,0,1 jdos.fe \
             -lt 2,col=1,0,0 dosp.dat  -lt 2 -tp 2~-0.0908,0,-0.0908,30
         Black and red lines are DOS computed through optinq and via SAVDOS=T, respectively.
         Blue and green lines are projections of the total DOS onto the Fe t2g and eg d states.

       * To confirm that the resolved DOS sum to the total DOS, do:
           mcx -br poptb.fe -split a 1,nr+1 1,2,'(nc-1)/3'+2,'(nc-1)/3*2'+2,'(nc-1)/3*3'+2 a11 a12 -csum -ccat jdos.fe -coll 1,2 -- -px:5
           mcx -br poptb.fe -split a 1,nr+1 1,2,'(nc-1)/3'+2,'(nc-1)/3*2'+2,'(nc-1)/3*3'+2 a11 a13 -csum -ccat jdos.fe -coll 1,3 -- -px:5
           mcx -br poptb.fe -split a 1,nr+1 1,2,'(nc-1)/3'+2,'(nc-1)/3*2'+2,'(nc-1)/3*3'+2 a11 a14 -csum -ccat jdos.fe -coll 1,4 -- -px:5

       * To plot the first k-point contribution to the total dos and the corresponding t2g- and eg- projections do:
           fplot -disp -frme 0,1,0,.7 -x -.3,.3 -y 0,.1 -lt 1 -colsy 2 -br poptb.fe \
           -lt 2,col=0,1,0 -colsy 2+140 -br poptb.fe -lt 2,col=1,0,0 -colsy 2+140+140 -br poptb.fe

       * After the DOS calculation the optics editor is run (--popted) to extract DOS(k) near the Fermi level (file pka.fe)

EOF
endif
else if ($ext == "gas") then
cat <<EOF
         1. Self-consistency with overlapping, space-filling spheres

         2. Energy levels tabulated for an autogenerated list of k-points ("qp" mode)
            that combines a dense mesh near k=0 and a uniform mesh over the BZ.

EOF
else if ($ext == "zbgan") then
if ($?quiet) then
cat <<EOF
         The ZB GaN test illustrates:

         1. Several joint density-of-states are generated

EOF
else
cat <<EOF
         1. The density is made self-consistent

         2. The following are generated and moved to these files:
                   file             contents
              jdos-samp.zbgan       Total joint DOS, by sampling integration
              jdos-tet3.zbgan       Total joint DOS, using tetwtq (accurate with ability to decompose DOS, but slow)
              pjdos-tet3.zbgan      Joint DOS resolved by (occ,unocc) pair (4 occ states * 7 unocc states)
              pjdos-tetk.zbgan      Joint DOS resolved by k (29 irreducible points)
              opt-tet3.zbgan        Im eps, tetrahedron integration
              opt-samp.zbgan        Im eps, sampling integration
              opt-abs.zbgan         Im eps with separate conduction and valence quasi-Fermi levels, similar to opt-samp.zbgan
              opt-lum.zbgan         Im eps simulating luminescence: scale transition rate by [1-f(E-imrefp)] * f(E-imrefn)

         You can verify that the resolved DOS sum to the total dos, e.g.
           mcx pjdos-tet3.zbgan -split a 1,nr+1 1,2,nc+1 a11 a12 -csum -ccat jdos.zbgan --

         After the test passes, try plotting the 29 k-resolved contributions to the JDOS with:
           fplot -frme 0,1,0,.7 -y 0,2 -lt 1 -colsy 2:nc pjdos-tetk.zbgan
         Plot the 28 (ib,jb)-resolved contributions to the JDOS with:
           fplot -frme 0,1,0,.7 -y 0,9 -lt 1 -colsy 2:nc pjdos-tet3.zbgan

         To confirm that the band-resolved contributions to Im eps (popt.zbgan)
         sum to the total for each polarization (opt-tet3.zbgan), do:
           mcx -vnpr=28 -qr popt.zbgan -split a 1,nr+1 1,2,2+npr,2+npr+npr,2+npr+npr+npr a11 a12 -csum a13 -csum a14 -csum -ccat -ccat -ccat opt-tet3.zbgan --

         The nonequilibrium tests employ distinct hole and electron quasi-Fermi levels (imref)
         (just below the valence band and above the conduction band, respectively).
         For absorption   (OPTICS_MODE=8) opt-abs.zbgan looks similar to opt-samp.zbgan (a finer k-mesh)
         For luminescence (OPTICS_MODE=9) opt-lum.zbgan has a sharp peak near the fundamental gap.

EOF
endif
if ($?MPI) then
  echo "$space not implemented in lmf-MPI ... skipping"
  goto chk1e
endif
else if ($ext == "gdn") then
cat <<EOF
         The GdN test also illustrates:

         1. the LDA+U implementation

         2. Artifically render f blocks of hamiltonian diagonal,
            with energy -1.5Ry in the majority channel and +1.4Ry in the minority

         3. Application of an additional external magnetic field with SO coupling
            in the LDA+U implementation.
            Bands are created with SO coupling, with no external B field,
            and again with a B field of 0.037Ry = 0.5eV added along the z axis.

            If you have plbnds and fplot installed, you can draw a picture
            comparing the bands in the two cases with the following after this test executes:

              echo -8,8,5,10 | plbnds -fplot -ef=0 -scl=13.6 -lbl=L,G,X,W,G -dat=blue bnds.nobfield.gdn
              echo -8,8,5,10 | plbnds -fplot -ef=0 -scl=13.6 -lbl=L,G,X,W,G -dat=dat  bnds.bfield.gdn
              mv plot.plbnds plot.plbnds~
              echo "% char0 colr=3,bold=4,clip=1,col=1,.2,.3" >>plot.plbnds
              echo "% char0 colb=2,bold=2,clip=1,col=.2,.3,1" >>plot.plbnds
              awk '{if (\$1 == "-colsy") {sub("-qr","-lt {colr} -qr");print;sub("dat","blue");sub("colr","colb");print} else {print}}' plot.plbnds~ >>plot.plbnds
              fplot -disp -f plot.plbnds

              See doc/External-Bfield.html for the figure generated by these bands,
              and some discussion of the shifts generated by the B field.

         4. An LDA calculation with the majority 4f levels in the core, and
         the minority 4f levels in the core with no charge (tests partial core occupancy)

         5. Charge mixing with no screening

         6. Eigenvectors saved in ASCII disk file eveca

EOF
else if ($ext == "tio2") then
cat <<EOF
         The TiO2 test illustrates:

         1. reading positions data from a site file

         2. a lattice relaxation, using Fletcher-Powell

EOF
else if ($ext == "cdte") then
if ($?quiet) then
cat <<EOF
         The CdTe test illustrates the LDA+U implementation in various modes.

EOF
else
cat <<EOF
         The cdte test illustrates the LDA+U implementation in various modes:

         It begins with a self-consistent LDA+U calculation with U=4 eV on Cd d, FLL limit.
         The main effects are:
           * to push down the Cd d levels by 2 eV relative to the LDA
           * reduce total energy is about 20mRy less binding relative to LDA (-0.389970 Ry)
             (there is presumably a corresponding shift in the free atom energy, but
             no attempt was made to calculate it)
           * increase the bandgap by 0.1 eV relative to LDA (0.52 eV) to 0.62 eV

         Starting from this density and density-matrix, three one-shot calculations are performed.

         1. A potential shift -2 eV on Cd d (IDU=4) is used in place of U=4 eV.
            The test verifies that both the total energy and the bandgap are unchanged,
            and that the density is almost self-consistent.

         2. An additional potential shift +1 eV on Cd s (IDU=4) is included.
            It has the effect of increasing the gap to 0.99 eV, and reducing the energy by 70 mRy.
            Ueff is determined to be 0.337 Ry to generate this potential shift.
            (Note that this pass somewhat reduces the Cd s density-matrix)

         3. A normal LDA+U is included on the Cd s orbital, U=0.337 in the FLL.
            The test verifies that the gap, small change in output density, and EHK
            are essentially identical to the results of test 2.

        As an additional check,
        a pass is made reading rst file with shear added and with --rs=101.


EOF
endif
else if ($ext == "cs") then
cat <<EOF
         The Cs illustrates the use of large spheres and small KMXV.
         Results of this test can be compared against the delta codes project.
         The volume-dependence of the total energy is in excellent agreement with Wien2k.

EOF
else if ($ext == "eras") then
if ($?quiet) then
cat <<EOF
         The ErAs test illustrates the LDA+U implementation in various modes.

EOF
else
cat <<EOF
         The ErAs test illustrates the LDA+U implementation for:

         1. a case when the LDA puts f orbitals at the Fermi energy.

         2. U is applied to phi,phidot on both d and f orbitals,
            (Compare to U applied to phi only, file fp/test/out.lmf.eras.idu2)

         3. convergence to a metastable solution with a reasonable spin moment
            but wrong orbital moment.

            For a much more stable solution, continue with
               cp $testdir/occnum2.eras occnum.eras
               rm mixm.eras wkp.eras dmats.eras
               lmf eras -vnit=30
            This aligns the orbital moment antiparallel to the spin moment.

            For a still more stable solution, continue with
               cp $testdir/occnum3.eras occnum.eras
               rm mixm.eras wkp.eras dmats.eras
               lmf eras -vnit=30
            This aligns the orbital moment parallel to the spin moment.

EOF
endif
else if ($ext == "er") then
if ($?quiet) then
cat <<EOF
         The Er test illustrates the LDA+U implementation in various modes.

EOF
cat <<EOF
         The Er test also illustrates the LDA+U implementation with:

           1. spin-orbit coupling

           2. the initial occupation number matrix given in spherical harmonics

           3. use of local orbitals in LDA+U

EOF
endif
else if ($ext == "cu") then
cat <<EOF
         The cu test illustrates:

         1.  high-lying local orbitals (Cu 5s,5p,4d are included as local orbitals)

         2.  METAL=5 for BZ integration

         3.  bands mode with APWs, pwmode 11 (see command-line argument --band in last lmf invocation)

         4.  Total DOS is generated.

EOF
if (! $?quiet) then
cat <<EOF
         To see the effect of APWs on the Cu energy band structure, do this following:
           zcat fp/test/bnds.cu-pwmode11 >bnds.lapw
           echo -10 40 5 10 | plbnds -lbl=L,G,X,W,G,K -ef=0 -scl=13.6 -fplot -dat=apw bnds.lapw
           zcat fp/test/bnds.cu >bnds.lmto
           echo -10 40 5 10 | plbnds -lbl=L,G,X,W,G,K -ef=0 -scl=13.6 -fplot -dat=lmto bnds.lmto
           awk 'BEGIN {print "%char0 lta=2,bold=3,col=1,0,0"} {if ($1 == "-colsy") {print; sub("lmto","apw"); sub("ltb","lta");print} else {print}}' plot.plbnds > plot.plbnds~
           fplot -disp -f plot.plbnds~

EOF
endif
else if ($ext == "na") then
cat <<EOF
         The na test illustrates:

         1.  comparison of the total energy for conventional and extended Na 2p orbitals

         2.  Resizing the interstitial mesh density

EOF
if (! $?quiet) then
cat <<EOF
             After the test finishes, compare the first three energies with:
               grep ^c out.lmf.na
             The third and fourth energy employ the same LMTO basis set,
             except that the fourth has a finer interstitial mesh density.

EOF
endif
else if ($ext == "srtio3") then
cat <<EOF
         The srtio3 test illustrates:

         1.  how to freeze augmented wave functions for a particular
             species (cf token FRZWF=t for species Sr).

         2.  local orbitals (Sr 4p and Ti 4p are included as local orbitals)

         3.  Compare conventional to extended local orbitals
             (last step recalculated with extended Sr 4p and Ti 4p local orbitals)

         4.  APW addition to basis

         5.  Renormalized free O atom for better starting density

         6.  use of restart file with shifted atom positions

         7.  the --rhopos switch

         8.  the --newvesmt switch

EOF
else if ($ext == "copt") then
cat <<EOF
         The copt test illustrates:

         1.  a spin-polarized case

         2.  Sampling with Fermi function at fairly high temperature (5mRy)
             Also demonstrates the LDA estimate of specific heat with the --cvK switch.

         3.  forces with correction mode 12 (FORCES=12)

         4.  two-kappa basis, with second kappa consisting of Co p only

         5.  inequivalent kmxa,lmxa,rmt

         6.  Perturbation treatment of core (LFOCA=2), and use spin-averaged potential for core (NMCORE=t)

         7.  Two-mode charge density mixing treatment

         For additional checks, modify fp/test/ctrl.copt and try the followng calculations:
             a. Set -vfancym=t .  Compare output to fp/test/out.lmf.fancymix.copt
                Invokes a complex charge mixing pattern (note: doesn't converge properly)

             b. Set -vnmcore=f .  Compare output to fp/test/out.lmf.nonmcore.copt
                Calculation is set to calculate the core in a spin averaged way, so that
                the up and down core densities are identical.
                This test shows the effect of treating the core spin polarized:
                The energy should drop by ~1 mRy.

             c. Set -vpwmode=1 .  Compare output to fp/test/out.lmf.pwmode1.copt
                Adds plane waves at Gamma=0 to basis.

             d. Set -vpwmode=11.  Compare output to fp/test/out.lmf.pwmode11.copt
                Adds q-dependent plane waves to basis.

                You should find the following total energies:
                  pwmode   etot
                     0    -1.9696
                     1    -2.0064
                    11    -2.0097

EOF
else if ($ext == "coptso") then
cat <<EOF
         The coptso test extends the CoPt test to:

         1.  test SO coupling, and evaluate site-resolved contribution to L.S

         2.  test SO coupling with L.S matrix elements scaled by 1/2

         The various estimates for the SO energy
         (change in H-F total energy, change in band sum, perturbation estimates) are all similar.
         A converged calculation requires a much finer k-mesh.

EOF
else if ($ext == "te") then
cat <<EOF
         The Te test illustrates:

         1.  a simple static relaxation of elemental Se with one symmetry-allowed degree of freedom

         2.  an spd-sp basis augmented by (1) floating orbitals, and (2) q-dependent APWs

         3.  Comparison of forces and total energy with and without floating orbitals,
             and with and without APWs.

         4.  PBE functional (XCFUN=4 GGA=4) through Kieron Burke's easypbe (internally compiled).
             Alternatively run this test with --libxc, which invokes the PBE functional via the libxc library.
             There are slight differences owing to differences in how gradients are constructed.
             The easypbe-libxc total energies differ by ~2mRy, probably because of differing core treatments,
             since the difference disappears if the same core density is used (see test 3).

         5.  Demonstrate local k-point shortening (SHORBZ=F)

EOF
if (! $?quiet) then
cat <<EOF
         lmf will first relax the atoms with the basis including floating orbitals.

         After relaxation, a new calculation is performed that removes floating orbitals
         while adding plane waves to the basis, so the total energy and forces may be compared.
         The basis size is variable, but averages ~80 orbitals, a little more than the floating
         orbitals case (~70 orbitals). About 3 mRy is gained relative to the floating orbitals case.

         Note that KMXA=5 is used with the PW calculation.  It isn't necessary in this case; still the
         user is cautioned to monitor this parameter when using energy cutoffs higher than 3 Ry or so.

         As a last step the calculation is repeated at the relaxed position with only atom-centered
         MTO's (neither floating orbitals nor plane waves).  Elimination of floating orbitals reduces
         the basis from 75 to 39 orbitals, and reduces the LDA total energy by about 4 mRy.

         The forces are not strongly affected by the local orbitals (or APWs), as can be seen by
         looking at the maximum force after the last step.

EOF
endif
else if ($ext == "zrt") then
cat <<EOF
         The zrt test also checks and illustrates the following:

         1.  the use of restart files in both binary and ascii forms

         2.  two-kappa basis

EOF
else if ($ext == "co") then
set strn = "Install the fplot plotting package to have this script create"
if ($?have_pldos && $?have_fplot) set strn = "This script will prompt you to see whether it should create "
cat <<EOF
         The co test also checks and illustrates the following:

         1.  a spin-polarized case

         2.  mixed tetrahedron/sampling method (METAL=4)

         3.  Constrained mixing (first spin is kept frozen, then charge, then both are allowed to change)
             Also SYMGRP_RHOPOS is used to ensure density remains positive definite

         4.  Broyden mixing

         5.  bands mode and spin texture weights (see command-line argument --band~col=5:9...)
EOF
if (! $?quiet) then
cat <<EOF
             SO coupling and color weights for projection onto Co d channels are included.
             Two separate weights are made (one each for d majority and minority bands.)
             $strn
             a figure from the bands file.

         6.  bands mode and spin texture weights (see command-line argument --band~colst=5:9...)
             The script will prompt you for a figure, if you want to see it.
             In this case the system is almost entire collinear; only a small region
             where majority and minority band come together will there be any deviation the z axis.
             (shown as red in the figure).
             To distinguish spin-up and spin-down bands, invoke the spin texture feature this way:
               echo -10,5,5,10 | plbnds -fplot~lt=1,col=0,0,0,colw=1,0,0,colw2=0,1,0~spintexture=z -ef=0 -scl=13.6 -lbl=M,G,A,L,G,K bnds.co-st

         7.  a few states are generated in q-list mode, with q rotated by an arbitrary rotation.
             The test is performed without SO coupling to demonstrate that eigenvalues are
             essentially identical whether rotated or not.
             If the test is repeated with SO coupling turned on, slight differences will appear.

EOF
else
cat <<EOF

EOF
endif
else if ($ext == "cr3si6") then
cat <<EOF
         The cr3si6 test also checks and illustrates the following:

         1.  insulator

         2.  forces with correction mode 1 (FORCES=1)

         3.  verbose output

         4.  Linearization energy floated using energy-weighted density matrix (Kotani style)

         5.  A case with low kmxa (kmxa=2)

         6.  Uses a 7-point quadrature for radial mesh (RQUAD=2)

         7.  switch --efrnge

         7.  OPTIONS_WRONSK=t

EOF

endif

set refout=$testdir/out.lmf.$ext testout=out.lmf.$ext
if ($?libxc) set refout=$testdir/out.lmf.pbe.$ext
if ($?MPIK && -e $testdir/out.lmf.$ext && ! ($?libxc)) set refout=$testdir/out.lmf.$ext
if ($?MPI && -e $testdir/out.lmf-MPI.$ext && ! ($?libxc)) set refout=$testdir/out.lmf-MPI.$ext

if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk1e
endif

set pass
#  echo DEBUG ; goto chk12
query chk11 chk1e 'run this test'
chk11:
# ... Look for executables
findcmd chk11a rdcmd "$path" "optional"
chk11a:
findcmd chk11b lmf "$path" "$topdir"
chk11b:
findcmd chk11c lmfa "$path" "optional"
chk11c:

# If output already available, just run through checks
if ($?haveout) goto chk12

# ... Setup: remove existing files and copy new ones
echo "$space rm -f {atm,qpp,eula,moms,mixm,rst,rsta,save,log,hssn,wkp,cv,bsmv,syml,bnds,dmats,dmats-save,occnum,sigm,qpp,jdos,popt,poptb,shfac,site,socscl,dos}.$ext specialspeca specialspecc $testout"
             rm -f {atm,qpp,eula,moms,mixm,rst,rsta,save,log,hssn,wkp,cv,bsmv,syml,bnds,dmats,dmats-save,occnum,sigm,qpp,jdos,popt,poptb,shfac,site,socscl,dos}.$ext specialspeca specialspecc $testout
if (! $?clean) then
echo "$space cp $cplst ."
             cp $cplst .
endif

# ... Run lmf program
if ($?MPIK && ! $?clean) then
  if ($?libxc) then
  egrep  '^TSTXCMK' ctrl.$ext >/dev/null
  if ($status) then
    echo "$space ... no category TSTXCMK ... skipping MPIK calculation"
    goto chk1e
  endif
  runrdcmd chk12 %11f $testout "-cat:TSTXCMK --noerr ctrl.$ext"
  else
  egrep  '^TMPKLMF' ctrl.$ext >/dev/null
  if ($status) then
    echo "$space ... no category TMPKLMF ... skipping MPI calculation"
    goto chk1e
  endif
  runrdcmd chk12 %11f $testout "-cat:TMPKLMF --noerr ctrl.$ext"
  endif
else if ($?MPI && ! $?clean) then
  if ($?libxc) then
  egrep  '^TSTXCMI' ctrl.$ext >/dev/null
  if ($status) then
    echo "$space ... no category TSTXCMI ... skipping MPI calculation"
    goto chk1e
  endif
  runrdcmd chk12 %11f $testout "-cat:TSTXCMI --noerr ctrl.$ext"
  else
  egrep  '^TMPILMF' ctrl.$ext >/dev/null
  if ($status) then
    echo "$space ... no category TMPILMF ... skipping MPI calculation"
    goto chk1e
  endif
  runrdcmd chk12 %11f $testout "-cat:TMPILMF --noerr ctrl.$ext"
  endif
else if (! $?clean) then
  if ($?libxc) then
  runrdcmd chk12 %11f $testout "-cat:TESTXC --noerr ctrl.$ext"
  endif
  runrdcmd chk12 %11f $testout "-cat:TESTLMF --noerr ctrl.$ext"
else
  if (-e ctrl.$ext) then
    runrdcmd chk11e %11f . "-cat:CLEAN --noerr ctrl.$ext"
chk11e:
  endif
  echo "$space rm -f ctrl.$ext semi.mater"
               rm -f ctrl.$ext semi.mater
  goto chk1e
endif
chk12:

# ... Extract total energies, forces, magnetic moments 1st and last iter
extract_res_n chk12a efa erfa "etot=" 2 0 etot=
chk12a:
set ehf1  =  `cat $testout | grep ehf= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1  =  `cat $testout | grep ehk= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehf1r =  `zcat $refout | grep ehf= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1r =  `zcat $refout | grep ehk= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set dq1   =  `cat $testout | grep 'RMS DQ=' | sed -n 1,1p | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dq1r   = `zcat $refout | grep 'RMS DQ=' | sed -n 1,1p | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqn   =  `cat $testout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`

if ($ext == "cr3si6") then
  set ik1  =  `cat $testout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("ik=","",$5); print $5}} '`
  set ik1r =  `zcat $refout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("ik=","",$5); print $5}} '`
  set nb1  =  `cat $testout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("band=","",$4); print $4}} '`
  set nb1r =  `zcat $refout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("band=","",$4); print $4}} '`
  set ek1  =  `cat $testout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("e1=","",$2); print $2}} '`
  set ek1r =  `zcat $refout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("e1=","",$2); print $2}} '`

  set ik2  =  `cat $testout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("ik=","",$9); print $9}} '`
  set ik2r =  `zcat $refout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("ik=","",$9); print $9}} '`
  set nb2  =  `cat $testout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("band=","",$8); print $8}} '`
  set nb2r =  `zcat $refout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("band=","",$8); print $8}} '`
  set ek2  =  `cat $testout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("e2=","",$6); print $6}} '`
  set ek2r =  `zcat $refout | grep 'BZWTS:' | awk '{ if (NR == 1) {gsub("e2=","",$6); print $6}} '`
endif

if ($ext == "coptso") then
   set ehf0 = `cat $testout | grep '^c' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehf1 = `cat $testout | grep '^c' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set ehf2 = `cat $testout | grep '^c' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,3p | tail -1`
   set ehk0 = `cat $testout | grep '^c' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehk1 = `cat $testout | grep '^c' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set ehk2 = `cat $testout | grep '^c' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,3p | tail -1`
   set mmom0 = `cat $testout | grep '^c' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set mmom1 = `cat $testout | grep '^c' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set mmom2 = `cat $testout | grep '^c' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,3p | tail -1`

   set sev0 = `extract-lines "START LMF " "Exit 0 LMF" 1 out.lmf.coptso | grep 'occ. bands' | tail -1 | awk '{print $4}'`
   set sev1 = `extract-lines "START LMF " "Exit 0 LMF" 3 out.lmf.coptso | grep 'occ. bands' | sed -n 1,1p | awk '{print $4}'`

   set ehfr0 = `zcat $refout | grep '^c' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehfr1 = `zcat $refout | grep '^c' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set ehfr2 = `zcat $refout | grep '^c' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,3p | tail -1`
   set ehkr0 = `zcat $refout | grep '^c' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehkr1 = `zcat $refout | grep '^c' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set ehkr2 = `zcat $refout | grep '^c' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,3p | tail -1`
   set mmomr0 = `zcat $refout | grep '^c' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set mmomr1 = `zcat $refout | grep '^c' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set mmomr2 = `zcat $refout | grep '^c' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,3p | tail -1`

  if ($?add0) then
    echo -n "         ..." ; $add0 $testout
  else if ($?poszer) then
    echo -n "         ..." ; $poszer $testout
  endif

  if (! $?quiet) then

  call zdiffiles chk1d0 "CPU -1 $testout $refout"
  chk1d0:

  echo "$space Harris energy,  no SO coupling    = $ehf0"
  echo "$space corresponding reference energy    = $ehfr0"
  set ediff = `echo $ehf0 $ehfr0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space H-K energy,     no SO coupling    = $ehk0"
  echo "$space corresponding reference energy    = $ehkr0"
  set ediff = `echo $ehk0 $ehkr0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space band sum,       no SO coupling    = $sev0"
  echo "$space band sum, 1st it with SO coupling = $sev1"
  set ediff = `echo $sev1 $sev0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space mag. moment,    no SO coupling   = $mmom0"
  echo "$space corresponding reference moment   = $mmomr0"
  set ediff = `echo $mmom0 $mmomr0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space Harris energy, with SO coupling   = $ehf1"
  echo "$space corresponding reference energy    = $ehfr1"
  set ediff = `echo $ehf1 $ehfr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space H-K energy,    with SO coupling   = $ehk1"
  echo "$space corresponding reference energy    = $ehkr1"
  set ediff = `echo $ehk1 $ehkr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space mag. moment,  with SO coupling   = $mmom1"
  echo "$space corresponding reference moment   = $mmomr1"
  set ediff = `echo $mmom1 $mmomr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space Harris energy, scaled SO coupling = $ehf2"
  echo "$space corresponding reference energy    = $ehfr2"
  set ediff = `echo $ehf2 $ehfr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space H-K energy,    scaled SO coupling = $ehk2"
  echo "$space corresponding reference energy    = $ehkr2"
  set ediff = `echo $ehk2 $ehkr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space mag. moment, scaled SO coupling  = $mmom2"
  echo "$space corresponding reference moment   = $mmomr2"
  set ediff = `echo $mmom2 $mmomr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  set ediff = `echo $ehf1 $ehf0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo -n "$space H-F diff,   full SO - no coupling = $ediff"
  set ediff = `echo $ehf2 $ehf0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo  "$space scaled SO - no coupling = $ediff"

  set ediff = `echo $ehk1 $ehk0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo -n "$space H-K diff,   full SO - no coupling = $ediff"
  set ediff = `echo $ehk2 $ehk0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo  "$space scaled SO - no coupling = $ediff"
  set ediff = `echo $mmom1 $mmom0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  set ediff = `echo $sev1 $sev0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space band sum change (force theorem)   = $ediff"

  set tso = `extract-lines "START LMF " "Exit 0 LMF" 3 out.lmf.coptso | grep 'Total SO coupling:' | tail -1 | awk '{print $9}'`
  echo -n "$space 2nd order pert calculation tso   = $tso"
  set tso = `extract-lines "START LMF " "Exit 0 LMF" 4 out.lmf.coptso | grep 'Total SO coupling:' | tail -1 | awk '{print $9}'`
  echo "       scaled pert calculation = $tso"

  set tso = `extract-lines "START LMF " "Exit 0 LMF" 3 out.lmf.coptso | grep 'Total SO coupling:' | tail -1 | awk '{print $5}'`
  echo -n "$space <L.S> 1st order estimate tso     = $tso"
  set tso = `extract-lines "START LMF " "Exit 0 LMF" 4 out.lmf.coptso | grep 'Total SO coupling:' | tail -1 | awk '{print $5}'`
  echo "       scaled 1st order est.   = $tso"

  echo
  echo -n "$space mom diff,   full SO - no coupling = $ediff"
  set ediff = `echo $mmom2 $mmom0  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo  "$space  scaled SO - no coupling = $ediff"

  endif # quiet

# pass checks

if ($?defatol1 == 0) set defatol1 = 2e-6
if ($?dehf1tol1 == 0) set dehf1tol1 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol1 == 0) set dmom1tol1 = 1e-5
if ($?dmomntol1 == 0) set dmomntol1 = 1e-5
if ($?dfmax1tol1 == 0) set dfmax1tol1 = 0.1
if ($?dfmaxntol1 == 0) set dfmaxntol1 = 0.1
if ($?drmsqtol1 == 0) set drmsqtol1 = 1e-5

echo

compare_res chk1da "ehf with   no SO coupling" $ehf0 $ehfr0 $dehf1tol1 pass
chk1da:

compare_res chk1db "ehk with   no SO coupling" $ehk0 $ehkr0 $dehf1tol1 pass
chk1db:

compare_res chk1dc "ehf with full SO coupling" $ehf1 $ehfr1 $dehf1tol1 pass
chk1dc:

compare_res chk1dd "ehk with full SO coupling" $ehk1 $ehkr1 $dehf1tol1 pass
chk1dd:

compare_res chk1de "ehf with full SO coupling" $ehf1 $ehfr1 $dehf1tol1 pass
chk1de:

compare_res chk1df "ehk with full SO coupling" $ehk1 $ehkr1 $dehf1tol1 pass
chk1df:

compare_res chk1dg "ehf with 1/2  SO coupling" $ehf2 $ehfr2 $dehf1tol1 pass
chk1dg:

compare_res chk1dh "ehk with 1/2  SO coupling" $ehk2 $ehkr2 $dehf1tol1 pass
chk1dh:

compare_res chk1di "mom with   no SO coupling" $mmom0 $mmomr0 $dmomntol1 pass
chk1di:

compare_res chk1dj "mom with full SO coupling" $mmom1 $mmomr1 $dmomntol1 pass
chk1dj:

compare_res chk1dk "mom with 1/2  SO coupling" $mmom2 $mmomr2 $dmomntol1 pass
chk1dk:

goto chk1e0

endif # coptso

egrep ' pwmode=[^0]' $testout >/dev/null
if (! $status) then
  set epw  =  `cat  $testout | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set epwr =  `zcat $refout  | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
endif

grep 'JDOS using' $testout >/dev/null
if (! $status) then
  set jdos
endif

grep 'Etot(LDA+U)' $testout >/dev/null
if (! $status) then
  set eldau  =  `cat $testout | grep 'Etot(LDA+U)' | tail -1 | awk '{print $NF}'`
  set eldaur =  `zcat $refout | grep 'Etot(LDA+U)' | tail -1 | awk '{print $NF}'`
endif

grep 'Maximum Harris force' $testout >/dev/null
if (! $status) then
  set fmax1  = `cat $testout | grep 'Maximum Harris force' | sed -n 1,1p | awk '{print $5}'`
  set fmax1r = `zcat $refout | grep 'Maximum Harris force' | sed -n 1,1p | awk '{print $5}'`
  set fmaxn  = `cat $testout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
  set fmaxnr = `zcat $refout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
endif

grep mmom= $testout >/dev/null
if (! $status) then
set mmom1  =  `cat $testout      | grep mmom= | sed -n 1,1p | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmom1r =  `zcat $refout | grep mmom= | sed -n 1,1p | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomn  =  `cat $testout      | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomnr =  `zcat $refout | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
endif

if ($ext == "te") then
set ebig  = `cat $testout | grep ^C | sed -n 1,1p | awk '{sub(".*ehf=","");sub("ehk=.*",""); print $0}'`
set ebigr = `zcat $refout | grep ^C | sed -n 1,1p | awk '{sub(".*ehf=","");sub("ehk=.*",""); print $0}'`
endif

set ediff = `echo $efa $erfa  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} printf "%10.2E", k}'`
if (! $?quiet) then
  echo " "
  echo "$space Total energy last free atom      = $efa"
  echo "$space Total energy of reference        = $erfa"
  echo "$space                    difference    =  $ediff"
  echo ' '

  echo "$space first iteration Harris energy    = $ehf1"
  echo "$space first iteration reference energy = $ehf1r"
  set ediff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo ' '

  echo "$space first iteration K-Sham energy    = $eks1"
  echo "$space first iteration reference energy = $eks1r"
  set ediff = `echo $eks1 $eks1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehf1 $eks1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`


  if ($?ik1) then
    echo
    echo "$space highest eval in band (iter 1)    = $ek1"  for band $nb1  found at ik=$ik1
    echo "$space corresponding reference          = $ek1r" for band $nb1r found at ik=$ik1r
    echo "$space lowest eval in band (iter 1)     = $ek2"  for band $nb2  found at ik=$ik2
    echo "$space corresponding reference          = $ek2r" for band $nb2r found at ik=$ik2r
  endif

  echo ' '
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmax1) then
  echo "$space first iteration maximum force    = $fmax1"
  echo "$space first iteration reference force  = $fmax1r"
  endif
  if ($?mmom1) then
    echo "$space first iteration magnetic moment  = $mmom1"
    echo "$space first iteration reference moment = $mmom1r"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif
  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo " "

  if ($?ebig) then
  echo "$space K-Sham energy, big basis         = $ebig"
  echo "$space corresponding reference energy   = $ebigr"
  set ediff = `echo $ebig $ebigr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo " "
  endif

  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  echo " "
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?ebig) then
  set ediff = `echo $ebig $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference large->small basis    = $ediff"
  echo " "
  endif

  if ($?epw) then
    echo "$space last iteration E(MTO + APW)      = $epw"
    echo "$space last iteration ref E(MTO + APW)  = $epwr"
    set ediff = `echo $epw $epwr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $ediff"
    echo " "
  endif

  if ($?eldau) then
    echo "$space last iteration Etot(LDA+U)       = $eldau"
    echo "$space last iteration ref E(LDA+U)      = $eldaur"
    set ediff = `echo $eldau $eldaur  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $ediff"
  endif

  if ($?fmaxn) then
  echo "$space last iteration maximum force     = $fmaxn"
  echo "$space last iteration reference force   = $fmaxnr"
  endif
  if ($?mmom1) then
  echo "$space last iteration magnetic moment   = $mmomn"
  echo "$space last iteration reference moment  = $mmomnr"
  set mdiff = `echo $mmomn $mmomnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"
  endif
  echo "$space last iter RMS input-output drho  = $dqn"
  echo "$space last iter reference RMS drho     = $dqnr"
  echo " "

  if ($ext == "te") then
    set esmall  = `cat $testout | grep ^c | tail -1 | awk '{sub(".*ehf=","");sub("ehk=.*",""); print $0}'`
    echo "$space Energy, MTO's only             = $esmall"
    echo "$space Energy, MTO's + floating orbs  = $ebig"
    echo "$space Energy, MTO's + APWs           = $epw"
  endif

  zcat $refout | grep RELAX >/dev/null
  if ($status == 0) then
    call showout chk13 RELAX
chk13:
    echo ' '
  endif

  if ($?add0) then
    echo -n "         ..." ; $add0 $testout
  else if ($?poszer) then
    echo -n "         ..." ; $poszer $testout
  endif
  call zdiffiles chk14 "CPU -1 $testout $refout"
chk14:
endif

if ($?have_pldos && $?have_fplot && $?slow) then
  egrep '^PLOTBND' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk14c
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk14a chk14c "generate a picture for energy bands"
  chk14a:
  runrdcmd chk14b %11f . "-cat:PLOTBND --noerr ctrl.$ext"
  chk14b:
   echo "$space $plot -disp -pr10 -f plot.plbnds"
                $plot -disp -pr10 -f plot.plbnds
  chk14c:

  egrep '^PLOTBST' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk14g
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk14d chk14g "generate a picture for energy bands with spin texture"
  chk14d:
  runrdcmd chk14e %11f . "-cat:PLOTBST --noerr ctrl.$ext"
  chk14e:
   echo "$space $plot -disp -pr10 -f plot.plbnds"
                $plot -disp -pr10 -f plot.plbnds
  chk14g:

endif


# ... Check that FA fit basis set is within tol of reference
#  echo "$fitbas" > tmp.$ext.ref
#  if ($?fitbas2) then
#    echo "$fitbas2" >> tmp.$ext.ref
#  endif
#  grep RSMH atm.$ext  > tmp.$ext.lmfa
#  cmp tmp.$ext.lmfa tmp.$ext.ref >/dev/null
#  set retval = $status
#  echo -n "$space lmfa fit basis identical to reference ? ..."
#  if ($retval == 0) then
#   echo yes
#  else
#    echo no
#    unset pass
#  endif

if ($?defatol1 == 0) set defatol1 = 2e-6
if ($?dehf1tol1 == 0) set dehf1tol1 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol1 == 0) set dmom1tol1 = 1e-4
if ($?dmomntol1 == 0) set dmomntol1 = 1e-4
if ($?dfmax1tol1 == 0) set dfmax1tol1 = 0.1
if ($?dfmaxntol1 == 0) set dfmaxntol1 = 0.1
if ($?drmsqtol1 == 0) set drmsqtol1 = 1e-4

# pass checks
chk1c:

echo "$space ... automatic pass checks :"

# ... Check that FA total energy is within tol of reference
compare_res chk1ca "FA etot (last species)" $efa $erfa $defatol1  pass
chk1ca:

compare_res chk1cb "1st  iter ehf" $ehf1 $ehf1r $dehf1tol1 pass
chk1cb:

if (! $?fmax1) goto chk1cc
compare_res chk1cc "1st  iter max force" $fmax1 $fmax1r $dfmax1tol1 pass
chk1cc:

if (! $?mmom1) goto chk1cd
compare_res chk1cd "1st  iter mmom" $mmom1 $mmom1r $dmom1tol1 pass
chk1cd:

compare_res chk1ce "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chk1ce:

if ($?eldau) then
compare_res chk1ci "last iter E(LDA+U)" $eldau $eldaur $dehf1toln pass
chk1ci:
endif

if ($?epw) then
compare_res chk1cj "last iter E(MTO+PW)" $epw $epwr $dehf1toln pass
chk1cj:
endif

if ($?fmaxn) then
compare_res chk1cf "last iter max force" $fmaxn $fmaxnr $dfmaxntol1 pass
chk1cf:
endif

if ($?mmomn) then
compare_res chk1cg "last iter mmom" $mmomn $mmomnr $dmomntol1 pass
chk1cg:
endif

compare_res chk1ch "last iter RMS dq" $dqn $dqnr $drmsqtol1 pass
chk1ch:

if ($?ik1) then
compare_res chk1cl "k index to highest occupied band" $ik1 $ik1r 0 pass
chk1cl:

compare_res chk1cm "k index to lowest unoccupied band" $ik2 $ik2r 0 pass
chk1cm:
endif

if ($ext == "copt") then
call zcmpnfiles chk1cn "6 cv.$ext $testdir/cv.$ext"
chk1cn:
echo -n "$space files cv.$ext and $testdir/cv.$ext equivalent to 6 digits? ... "
if ($retval == 0) then
  echo  yes
else
  echo no
  unset pass
endif
endif

# compare jdos to reference
if ! ($?jdostol) set jdostol = .0001
if ! ($?epstol) set epstol = .0001
if (-e $testdir/jdos-samp.$ext) then
zcmpmfiles_res_0 chk1cka "Max deviation in jdos-samp.$ext from reference $testdir/jdos-samp.$ext" $jdostol pass 4 jdos-samp.$ext $testdir/jdos-samp.$ext
chk1cka:
endif
if (-e $testdir/jdos-tet.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 jdos-tet.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer jdos-tet.$ext
endif
zcmpmfiles_res_0 chk1ckb "Max deviation in jdos-tet.$ext from reference $testdir/jdos-tet.$ext" $jdostol pass 4 jdos-tet.$ext $testdir/jdos-tet.$ext
chk1ckb:
endif
if (-e $testdir/jdos-tet3.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 jdos-tet3.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer jdos-tet3.$ext
endif
zcmpmfiles_res_0 chk1ckc "Max deviation in jdos-tet3.$ext from reference $testdir/jdos-tet3.$ext" $jdostol pass 4 jdos-tet3.$ext $testdir/jdos-tet3.$ext
chk1ckc:
endif
if (-e $testdir/pjdos-tet3.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 pjdos-tet3.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer pjdos-tet3.$ext
endif
zcmpmfiles_res_0 chk1ckd "Max deviation in pjdos-tet3.$ext from reference $testdir/pjdos-tet3.$ext" $jdostol pass 4 pjdos-tet3.$ext $testdir/pjdos-tet3.$ext
chk1ckd:
endif
if (-e $testdir/pjdos-tetk.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 pjdos-tetk.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer pjdos-tetk.$ext
endif
zcmpmfiles_res_0 chk1cke "Max deviation in pjdos-tetk.$ext from reference $testdir/pjdos-tetk.$ext" $jdostol pass 4 pjdos-tetk.$ext $testdir/pjdos-tetk.$ext
chk1cke:
endif
if (-e $testdir/dos-tet3.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 dos-tet3.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer dos-tet3.$ext
endif
set refb = $testdir/dos-tet3.$ext
if ($?libxc) set refb = $testdir/dos-tet3.pbe.$ext
zcmpmfiles_res_0 chk1ckf "Max deviation in dos-tet3.$ext from reference $refb" $jdostol pass 4 dos-tet3.$ext $refb
chk1ckf:
endif
if (-e $testdir/pdos-tetk.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 pdos-tetk.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer pdos-tetk.$ext
endif
zcmpmfiles_res_0 chk1ckg "Max deviation in pdos-tetk.$ext from reference $testdir/pdos-tetk.$ext" $jdostol pass 4 pdos-tetk.$ext $testdir/pdos-tetk.$ext
chk1ckg:
endif

if (-e $testdir/opt-tet3.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 opt-tet3.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer opt-tet3.$ext
endif
zcmpmfiles_res_0 chk1ckh "Max deviation in opt-tet3.$ext from reference $testdir/opt-tet3.$ext" $epstol pass 4 opt-tet3.$ext $testdir/opt-tet3.$ext
chk1ckh:
endif

if (-e $testdir/opt-samp.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 opt-samp.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer opt-samp.$ext
endif
zcmpmfiles_res_0 chk1cki "Max deviation in opt-samp.$ext from reference $testdir/opt-samp.$ext" $epstol pass 4 opt-samp.$ext $testdir/opt-samp.$ext
chk1cki:
endif

if (-e $testdir/opt-abs.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 opt-abs.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer opt-abs.$ext
endif
zcmpmfiles_res_0 chk1cki2 "Max deviation in opt-abs.$ext from reference $testdir/opt-abs.$ext" $epstol pass 4 opt-abs.$ext $testdir/opt-abs.$ext
chk1cki2:
endif

if (-e $testdir/opt-lum.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 opt-lum.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer opt-lum.$ext
endif
zcmpmfiles_res_0 chk1cki3 "Max deviation in opt-lum.$ext from reference $testdir/opt-lum.$ext" $epstol pass 4 opt-lum.$ext $testdir/opt-lum.$ext
chk1cki3:
endif

if ! ($?dostol) set dostol = .0001
if (-e dos.$ext && -e $testdir/dos.$ext) then
if ($?add0) then
  echo -n "         ..." ; $add0 dos.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer dos.$ext
endif
zcmpmfiles_res_0 chk1ckj "Max deviation in dos.$ext from reference $testdir/dos.$ext" $dostol pass 4 dos.$ext $testdir/dos.$ext
chk1ckj:
endif

# if ! ($?bndtol) set bndtol = .0001
# if (-e bnds.$ext && -e $testdir/bnds.$ext) then
# if ($?add0) then
#   echo -n "         ..." ; $add0 bnds.$ext
# else if ($?poszer) then
#   echo -n "         ..." ; $poszer bnds.$ext
# endif
# zcmpmfiles_res_0 chk1ckk "Max deviation in bnds.$ext from reference $testdir/bnds.$ext" $bndstol pass 4 bnds.$ext $testdir/bnds.$ext
# chk1ckk:
# endif

if ($ext == "fe") then
  if ($?have_mc) then
    echo ' '
    echo "$space ... compare DOS at E=.05 Ry (near Ef) to sum of k-resolved DOS, to sum of k-resolved DOS in FBZ:"
    echo  \
   'mcx -f6f12.6 jdos.fe -rowl 501 -coll 1,2 -br poptb.fe -rowl 501 -coll 2:141 -csum pka.fe -csum -rsum -s1/16^3 -ccat -ccat -ccat -e6 x1 x2 x3 x4 x3-x2 x4-x3| grep -v rows'
    echo '       E          DOS        sum(k)     sum(FBZ)  sum(k)-DOS  sum(FBZ)-sum(k)'
    mcx -f6f12.6 jdos.fe -rowl 501 -coll 1,2 -br poptb.fe -rowl 501 -coll 2:141 -csum pka.fe -csum -rsum -s1/16^3 -ccat -ccat -ccat -e6 x1 x2 x3 x4 x3-x2 x4-x3| grep -v rows
    set sumk_diff   = `mcx -f6f12.6 jdos.fe -rowl 501 -coll 1,2 -br poptb.fe -rowl 501 -coll 2:141 -csum pka.fe -csum -rsum -s1/16^3 -ccat -ccat -ccat -e6 x1 x2 x3 x4 x3-x2 x4-x3 | tail -1 | awk '{print $5}'`
    set sumfbz_diff = `mcx -f6f12.6 jdos.fe -rowl 501 -coll 1,2 -br poptb.fe -rowl 501 -coll 2:141 -csum pka.fe -csum -rsum -s1/16^3 -ccat -ccat -ccat -e6 x1 x2 x3 x4 x3-x2 x4-x3 | tail -1 | awk '{print $6}'`
  endif
endif

# compare qp list and bnds file to reference
if ($ext == "gas") then
set ndig = 8
call zcmpnfiles chk1ckl "$ndig qp.$ext $testdir/qp.$ext"
chk1ckl:
echo -n "$space ... files qp.$ext and $testdir/qp.$ext equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
if ($?add0) then
  echo -n "         ..." ; $add0 bnds.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer bnds.$ext
endif
set ndig = 4
call zcmpnfiles chk1ckm "$ndig bnds.$ext $testdir/bnds.$ext"
chk1ckm:
echo -n "$space ... files bnds.$ext and $testdir/bnds.$ext equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif

if ($ext == "gdn") then
echo " "
set ndig = 4
if ($?add0) then
  echo -n "         ..." ; $add0 bnds.nobfield.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer bnds.nobfield.$ext
endif
call zcmpnfiles chk1ckn "$ndig bnds.nobfield.$ext $testdir/bnds.nobfield.$ext"
chk1ckn:
echo -n "$space ... files bnds.nobfield.$ext and $testdir/bnds.nobfield.$ext equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
if ($?add0) then
  echo -n "         ..." ; $add0 bnds.bfield.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer bnds.bfield.$ext
endif
call zcmpnfiles chk1cko "$ndig bnds.bfield.$ext $testdir/bnds.bfield.$ext"
chk1cko:
echo -n "$space ... files bnds.bfield.$ext and $testdir/bnds.bfield.$ext equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
set ndig = 4
if ($?add0) then
  echo -n "         ..." ; $add0 bnds.fdiagonal.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer bnds.fdiagonal.$ext
endif
call zcmpnfiles chk1ckp "$ndig bnds.fdiagonal.$ext $testdir/bnds.fdiagonal.$ext"
chk1ckp:
echo -n "$space ... files bnds.fdiagonal.$ext and $testdir/bnds.fdiagonal.$ext equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
echo " "
endif

if ($ext == "fe" && ! $?libxc) then
if ($?add0) then
  echo -n "         ..." ; $add0 pka.$ext
  echo -n "         ..." ; $add0 opt.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer pka.$ext
  echo -n "         ..." ; $poszer opt.$ext
endif
zcmpmfiles_res_0 chk1ckq "Max deviation in pka.$ext from reference $testdir/pka.$ext" 0.03 pass 4 pka.$ext $testdir/pka.$ext
chk1ckq:
if ($?sumk_diff) then
compare_res_0 chk1ckr "Difference sum(k)-total DOS" $sumk_diff .0001 pass
chk1ckr:
endif
if ($?sumfbz_diff) then
compare_res_0 chk1cks "Difference FBZ IBZ DOS" $sumfbz_diff .00001 pass
chk1cks:
endif
zcmpmfiles_res_0 chk1ckt "Max deviation in opt.$ext from reference $testdir/opt.$ext" 0.03 pass 4 opt.$ext $testdir/opt.$ext
chk1ckt:
zcmpmfiles_res_0 chk1cku "Max deviation in opt-so1.$ext from reference $testdir/opt-so1.$ext" 0.03 pass 4 opt-so1.$ext $testdir/opt-so1.$ext
chk1cku:
zcmpmfiles_res_0 chk1ckv "Max deviation in opt-so3.$ext from reference $testdir/opt-so3.$ext" 0.03 pass 4 opt-so3.$ext $testdir/opt-so3.$ext
chk1ckv:
endif

# compare bnds to reference
if (-e bnds.$ext) then
if ! ($?bndstol) set bndstol = 1e-4
# echo 'DEBUG'; # set verbose
#  zcmpmfiles_res_0 chk1ck "... Max deviation in bnds.$ext from reference" $bndstol pass 4 bnds.$ext $testdir/bnds.$ext
#  chk1ck:
set ndig = 4
set refb = $testdir/bnds.$ext
if ($?libxc) set refb = $testdir/bnds.pbe.$ext
call zcmpnfiles chk1ck "$ndig bnds.$ext $refb"
chk1ck:
echo -n "$space ... files bnds.$ext and $refb equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot < 1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif

if (-e bnds.$ext-pwmode11) then
if ! ($?bndstol) set bndstol = 1e-4
set ndig = 4
call zcmpnfiles chk1ck2 "$ndig bnds.$ext-pwmode11 $testdir/bnds.$ext-pwmode11"
chk1ck2:
echo -n "$space ... files bnds.$ext-pwmode11 and $testdir/bnds.$ext-pwmode11 equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot < 1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif
endif

# Comparison must be tested against absolute value ...
# if (-e bnds.$ext-st) then
# if ! ($?bndstol) set bndstol = 1e-4
# set ndig = 4
# call zcmpnfiles chk1ck3 "$ndig bnds.$ext-st $testdir/bnds.$ext-st"
# chk1ck3:
# echo -n "$space ... files bnds.$ext-st and $testdir/bnds.$ext-st equivalent to $ndig digits? ... "
# if ($retval == 0) then
#   echo  yes
# else
#   if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot < 1.)}'` == 1) then
#     echo ok "($retval difference(s) of $ncharfile)"
#   else
#     echo no "($retval difference(s) remaining of $ncharfile)"
#     unset pass
#   endif
# endif
# endif
# endif

if (-e bnds.$ext-rot) then
if ! ($?bndstol) set bndstol = 1e-4
set ndig = 4
call zcmpnfiles chk1ck4 "$ndig bnds.$ext-rot $testdir/bnds.$ext-rot"
chk1ck4:
echo -n "$space ... files bnds.$ext-rot and $testdir/bnds.$ext-rot equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot < 1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif
endif

if (-e bnds.$ext-st) then
if ! ($?bndstol) set bndstol = 1e-4
set ndig = 4
call zcmpnfiles chk1ck5 "$ndig bnds.$ext-st $testdir/bnds.$ext-st"
chk1ck5:
echo -n "$space ... files bnds.$ext-st and $testdir/bnds.$ext-st equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot < 1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif
endif

chk1e0:

if ($?pass) then
    echo "$space test 1 PASSED ($ext)"
    set jobspassed = ($jobspassed 1)
else
    echo "$space test 1 FAILED ($ext)"
    set failed = ($failed 1)
endif

chk1e:


echo $joblist | egrep '\b2\b' >/dev/null
if ($status) goto chk2e
set jobid = 2
if ($ext == "srtio3") then
cat <<EOF

         --- Test 2.  Check rotation of lattice, and density ---

         SrTiO3 case tests lattice rotation, including rotation of local densities.
         The lattice is rotated 90 degrees around y, mapping x into z and z into -x.
         The density should remain self-consistent and forces be unchanged, mutatis mutandis.
         Run this test only after completing test 1.

EOF
else
cat <<EOF

         --- Test 2.  Core-level spectroscopy (EELS), Mulliken analysis, partial DOS ---

EOF
endif
set strn = "Install the fplot plotting package to have this script create"
if ($?have_pldos && $?have_fplot) set strn = "After the calculation completes, this script will ask whether you want to generate"
if ($?quiet) then
else if ($ext == "crn") then
cat <<EOF
         The CrN test case generates core-level spectroscopy for the
         1s state in N.  The self-consistent calculation proceeds with
         an electron missing from the N 1s core, which corresponds to
         the 'sudden approximation' (system relaxes instantanously
         from electron ejected out of a core hole).

         $strn
         a figure showing the core level spectroscopy

EOF
else if ($ext == "gas") then
cat <<EOF

         Demonstration and check of total and partial DOS. This switch is passed to lmf:
           --pdos~lcut=2,2,1,1~mode=1
         It generates l-resolved DOS in the Ga spd, As spd and empty sphere sp channels,
         from partial waves in augmentation spheres.
         Partial DOS are stored in file pdos.gas and total DOS in tdos.gas.

         Overlapping space-filling spheres enable comparison of partial DOS to ASA.

EOF
if ($?have_pldos && $?have_fplot && $?slow) then
cat <<EOF
         To see what fraction of the total DOS is carried by the Ga and As spd channels,
         invoke the following after this script executes:
           rdcmd '-f:#rdcmd:%2f' -cat:CMPDOS --noerr ctrl.gas
           ghostview fplot.ps

         To see how well the sum of spd partial waves on Ga+As+ES sites with overlapping
         spheres account for the total integrated DOS, you can redo generating the DOS
         including d orbitals on the empty spheres.  The following script does this, and
         compares the integrated DOS to the total integrated DOS
           rdcmd '-f:#rdcmd:%2f' -cat:TSTIDOS --noerr ctrl.gas
           ghostview fplot.ps

EOF
endif

else if ($ext == "cr3si6") then
cat <<EOF
         The cr3si6 test demonstrates lmf's ability to generate m-resolved site DOS using a
         symmetry-reduced k-mesh.  To accumulate the DOS, lmf rotates the eigenvectors through
         points in the star of each irreducible point

         1.  Partial DOS on the Cr and Si sites are initially generated with symmetry operations suppressed.
             (Note: the calculation is set up to average over the 3 Cr sites and 6 Si sites ... not true DOS!)

         2.  The calculation is repeated with symmetry operations included.  For numerical reasons the
             two DOS are slightly different.  However they can be seen to be nearly identical
             by comparing pictures.

        To compare the different DOS pictorially, after the test executes do
          echo 15 4 -8 8 | pldos -ef=0 -escl=13.6 -fplot~ext=1 -lst='1;2,3,4;5;6;7;8;9;10;11;12;13' pdos1.cr3si6
          echo 15 4 -8 8 | pldos -ef=0 -escl=13.6 -fplot~ext=2 -lst='1;2,3,4;5;6;7;8;9;10;11;12;13' pdos.cr3si6
          fplot -pr10 -f fp/test/plot.cr3si6dos
          gs fplot.ps
        The panels correspond to :
              1      Cr s
              2      Cr px+py+pz
              3..7   the 5 different Cr d states
              8      Si s
              9..11  the three Si p states.
        Note that (-m,m) DOS appear the same : compare panels (3,7) or (4,6) or (9,11).

EOF

else if ($ext == "co") then
cat <<EOF
         The Co test case illustrates partial dos resolved by both l and m.

         $strn
         a figure of l- and spin- resolved DOS

EOF
else if ($ext == "gdn") then
cat <<EOF
         The GdN test case illustrates a special use of Mulliken analysis with spin-orbital coupling

         The following files are created:
           dos-mull.gdn  : total DOS resolved into spin components
           tdos.gdn      : total DOS, generated with SO coupling, without resolution into partial dos
           dos.nosym.gdn : site-, l-, and m- resolved Mulliken DOS, without assuming symmetry operations
           dos.sym       : site-, l-, and m- resolved Mulliken DOS, using symmetry operations
           tdos.sym      : sum of all partial DOS (should be equivalent to tdos.gdn)
           tdos1.gdn     : extracted from tdos.gdn: data in 2-D array format to compare with tdos.sym

         Several checks/demonstrations are performed:

         1. The total DOS can be drawn with color weights to distinguish spins. It makes use of
            spin information in dos-mull.gdn.  A prompt will ask you whether to draw the DOS.
            (Figure shows the majority f portion in red, minority f portion in green)

         2. The m-resolved Mulliken DOS is calculated with and without symmetry operations,
            to confirm the symmetrizer works in the SO case.
            A prompt will ask you whether to draw the DOS in the xy channel

         3. The channels in the partial DOS is are summed to compare against the total DOS.

EOF
else if ($ext == "fe") then
cat <<EOF
         The fe test case illustrates both Mulliken analysis resolved by l and m and core-level spectroscopy.

         This script will create a picture of the DOS resolved into spin and three groups of orbitals:
         (s+p), (d states of t_2 symmetry), (d states of e_g symmetry).
         (Caution: potential is not self-consistent!)

         Mulliken analysis is also useful for the spin-coupled case, enabling
         the resolution of total DOS into spin components.  The following commands will generate
         a picture of the total DOS, where the e2 part of the d channel is colored in red.

           lmf --rs=0 -vso=t --mull:mode=2 -vnk=6 -vnit=1 fe
           mv dos.fe tdos.fe
           lmdos --nosym -vso=t --mull:mode=2 --dos:npts=1001:window=-.7,.8 -vnk=6 fe
           mv dos.fe dos-mull.fe
           echo 40 7 -9 10 | pldos -ef=0 -escl=13.6 -fplot '-lst=13,17' -ref:fn=tdos.fe:chan=1:scale dos-mull.fe

EOF
endif

if ($ext == "srtio3") then
set refout=$testdir/out.lmf-rot.$ext testout=out.lmf-rot.$ext
else
set refout=$testdir/out.lmf-dos.$ext testout=out.lmf-dos.$ext
if ($?MPIK && -e $testdir/out.lmf-dos.$ext) set refout=$testdir/out.lmf-dos.$ext
if ($?MPI && -e $testdir/out.lmf-MPI-dos.$ext) set refout=$testdir/out.lmf-MPI-dos.$ext
endif
#set pass; echo 'TEST'; goto chk256e
set pass
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk2e
endif

if ($ext == "gdn" && ! $?clean) then
if (! -e dmats-save.gdn || ! -e rst-save.gdn) then
  echo "$space ... skipping test : missing one or both of files dmats-save.gdn rst-save.gdn "
  echo "$space     Run test 1 for GdN to generate them automatically."
  goto chk2e
endif
endif
if ($ext == "crn" || $ext == "co" || $ext == "fe" || $ext == "cr3si6") then
else
if ($?MPIK || $?MPI) then
  echo "$space ... no MPI-specific checks ... skipping this calculation"
  goto chk2e
endif
endif

query chk21 chk2e 'run this test'
chk21:
# ... Look for executables
findcmd chk21a rdcmd "$path" "optional"
chk21a:
findcmd chk21b lmf "$path" "$topdir"
chk21b:
findcmd chk21c lmfa "$path" "optional"
chk21c:
findcmd chk21d lmdos "$path" "optional"
chk21d:

# If output already available, just run through checks
if ($?haveout) goto chk23

if ($ext == "srtio3") then
  if (! -e rst.$ext || ! -e ctrl.$ext) then
    echo "$space ... missing ctrl.$ext or rst.$ext ... complete test 1 before running this one"
    echo "$space ... skipping test"
    goto chk2e
  endif
else
echo "$space cp $cplst ."
             cp $cplst .
echo "$space rm -f {atm,dos,mixm,rst,save,log,hssn,wkp,dos,tdos,pdos,dos-mull,qpp}.$ext $testout"
             rm -f {atm,dos,mixm,rst,save,log,hssn,wkp,dos,tdos,pdos,dos-mull,qpp}.$ext $testout
endif

# For CLS
if ($ext == "srtio3") then
  runrdcmd chk22 %11f $testout "-cat:TESTROT --noerr ctrl.$ext"
else if ($?MPIK && ! $?clean) then
  egrep '^MPIKCLS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk25
  runrdcmd chk22 %11f $testout "-cat:MPIKCLS --noerr ctrl.$ext"
else if ($?MPI && ! $?clean) then
  egrep '^MPICLS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk25
  runrdcmd chk22 %11f $testout "-cat:MPICLS --noerr ctrl.$ext"
else if (! $?clean) then
  egrep '^TESTCLS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk25
  runrdcmd chk22 %11f $testout "-cat:TESTCLS --noerr ctrl.$ext"
else
  if (-e ctrl.$ext) then
  runrdcmd chk2e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk2e
endif

chk22:
echo ' '
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif
call zdiffiles chk23 "CPU -1 $testout $refout"
chk23:

if ($ext == "srtio3") goto chk28

if ($?have_pldos && $?have_fplot && $?slow && ! $?clean) then
  egrep '^PLOTCLS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk24c
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk24a chk24c "generate a picture for CLS"
  chk24a:
  runrdcmd chk24b %11f . "-cat:PLOTCLS --noerr ctrl.$ext"
  chk24b:
   echo "$space $plot -disp -pr10 -f plot.dos"
                $plot -disp -pr10 -f plot.dos
  chk24c:

  egrep '^PLOTMUL' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk24f
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk24d chk24f "generate a picture for Mulliken DOS"
  chk24d:
  runrdcmd chk24e %11f . "-cat:PLOTMUL --noerr ctrl.$ext"
  chk24e:
  echo "$space $plot -disp -pr10 -f plot.dos"
               $plot -disp -pr10 -f plot.dos

  if ($ext != "gdn") goto chk24f
  query chk24i chk24f "generate a picture comparing xy channel of Mulliken DOS, with & without symops"
  chk24i:
  echo "$space $plot -disp -colsy 10,11 dos.gdn -lt 2,col=1,0,0 -colsy 10,11 dos.nosym.gdn"
               $plot -disp -colsy 10,11 dos.gdn -lt 2,col=1,0,0 -colsy 10,11 dos.nosym.gdn
  chk24f:

if ($ext == "srtio3") then
else if ($ext == "fe") then
cat <<EOF

   Consider contracting the lm-resolved FP Mulliken dos, and comparing to the ASA Mulliken dos, e.g:
   (Caution: fp case is not self-consistent)
      echo 20 10 -0.7 .8 | pldos -fplot -lst="1;3:7:2;9:17:2" -lst2 dos-mull.fe
      fplot -disp -pr10 -f plot.dos
      cp fplot.ps fp.ps
      zcat $testdir/../../testing/dos-mull.fe >dos.fe
      echo 20 10 -0.7 .8 | pldos -fplot -lst="1;3;5" -lst2 dos.fe
      fplot -disp -pr10 -f plot.dos
      cp fplot.ps asa.ps

EOF
endif # materials

endif # plotting branch

echo "$space ... automatic pass checks :"

# compare dos-cls to reference
if (! -e "$testdir/dos-cls.$ext") goto chk256e
if ! ($?dosclstol) set dosclstol = 1e-4
call zcmpmfiles chk256 "6 dos-cls.$ext $testdir/dos-cls.$ext"
# call zcmpnfiles chk256 "6 dos-cls.$ext $testdir/dos-cls.$ext"
chk256:
compare_res_0 chk25a "... Max deviation in dos-cls from reference" $retval $dosclstol pass
chk25a:
# echo -n "$space ... files dos-cls.$ext and $testdir/dos-cls.$ext equivalent to 6 digits? ... "
# if ($retval == 0) then
#   echo  yes
# else
#   call zcmpnfiles chk253 "4 dos-cls.$ext $testdir/dos-cls.$ext"
#   chk253:
#   echo -n "no ... to 4 digits? ... "
#   if ($retval == 0) then
#     echo yes
#   else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
#     echo ok "($retval difference(s) of $ncharfile)"
#   else
#     echo no "($retval difference(s) remaining of $ncharfile)"
#     unset pass
#   endif
# endif
chk256e:

# Compare dos-mull to reference
if (! -e "$testdir/dos-mull.$ext") goto chk25x
if ! ($?dosmulltol) set dosmulltol = 1e-4
# call zcmpmfiles chk266 "6 dos-mull.$ext $testdir/dos-mull.$ext"
# chk266:
# compare_res_0 chk26a "Max deviation in dos-mull from reference" $retval $dosmulltol pass
# chk26a:
if ($?add0) then
  echo -n "         ..." ; $add0 dos-mull.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer dos-mull.$ext
endif
# if ($?have_mc) then
# # set mcexcl = "~START"  mcterm = "~CPU"
# # zcmpmcx chk29 "Max deviation in dos-mull.$ext" 1e-4 -1 1e-4 pass dos-mull.$ext $testdir/dos-mull.$ext
# # zcmpmcx chk266a "" 5e-3 40 1e-4 pass dos-mull.$ext $testdir/dos-mull.$ext
# # zcmpmcx chk266a "" -1e-1 20 1e-4 null dos-mull.$ext $testdir/dos-mull.$ext
# zcmpmcx chk266a "... checking" -1 40 1e-4 pass dos-mull.$ext $testdir/dos-mull.$ext
# endif
call zcmpnfiles chk266 "6 dos-mull.$ext $testdir/dos-mull.$ext"
chk266:
echo -n "$space ... files dos-mull.$ext and $testdir/dos-mull.$ext equivalent to 6 digits? ... "
if ($retval == 0) then
  echo  yes
else
  call zcmpnfiles chk263 "3 dos-mull.$ext $testdir/dos-mull.$ext"
  chk263:
  echo -n "no ... to 3 digits? ... "
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
chk266a:

if (-e dos-mull-nosym.$ext) then
call zcmpmfiles chk267 "6 dos-mull-sym.fe dos-mull-nosym.fe"
chk267:
set dossymdif = 0.5
compare_res_0 chk267a "... Max diff dos-mull-sym.$ext dos-mull-nosym.$ext" $retval $dossymdif pass
chk267a:
endif

chk25x:
if ($ext == "gdn" && ! $?clean) then
  if ($?have_mc) then
     set maxdif = `mcx -qr tdos1.gdn -qr tdos.sym -- -inc 'i>1&i<nr' -abs -max:g -w:nohead . | awk '{print $3}'`
     compare_res_0 chk25x1 "Difference sum(chan)PDOS-total DOS" $maxdif .0001 pass
     chk25x1:
  endif

  set tdostol = 1e-2  dostol = 1e-3  nl = 100
  zcmpmfiles_res_tol chk25x2 "Max deviation in tdos.gdn from reference $testdir/tdos.$ext" $tdostol pass 8 tdos.$ext $testdir/tdos.$ext 0
  chk25x2:

  zcmpmfiles_res_tol chk25x3 "Max deviation ($nl lines) in dos.gdn from reference $testdir/dos.sym" $dostol pass 8 dos.gdn $testdir/dos.sym $nl
  chk25x3:
endif

# --- check partial DOS ---
chk25:
egrep '^TSTPDOS' ctrl.$ext >/dev/null
set retval = $status
if ($retval != 0) goto chk29
if ($?MPIK && ! $?clean) then
  egrep '^TSMPDOS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk29
  runrdcmd chk26 %11f $testout "-cat:TSMPDOS --noerr ctrl.$ext"
else if (! $?clean) then
  runrdcmd chk26 %11f $testout "-cat:TSTPDOS --noerr ctrl.$ext"
else
  if (-e ctrl.$ext) then
  runrdcmd chk2e %11f $testout "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk2e
endif
chk26:

echo ' '
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif
call zdiffiles chk27 "CPU -1 $testout $refout"
chk27:
if ($?have_pldos && $?have_fplot && $?slow) then
  egrep '^PLOTDOS' ctrl.$ext >/dev/null
  set retval = $status
  if ($retval != 0) goto chk27c
  echo "$space ... using $plot for fplot"
  echo "$space ... using $pldos for pldos"
  query chk27a chk27c "generate a picture for DOS"
  chk27a:
  runrdcmd chk27b %11f . "-cat:PLOTDOS --noerr ctrl.$ext"
  chk27b:
   echo "$space $plot -disp -pr10 -f plot.dos"
                $plot -disp -pr10 -f plot.dos
  chk27c:
endif

# Compare pdos to reference
if ! ($?pdostol) set pdostol = 1e-4
if ($?add0) then
  echo -n "         ..." ; $add0 pdos.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer pdos.$ext
endif
call zcmpmfiles chk274 "4 pdos.$ext $testdir/pdos.$ext"
call zcmpnfiles chk274 "4 pdos.$ext $testdir/pdos.$ext"
chk274:
compare_res_0 chk274a "... Max deviation in pdos.$ext from reference $testdir/pdos.$ext" $retval $pdostol pass
chk274a:

# SrTiO3 checks
chk28:

if ($ext == "srtio3") then
set dq1 = `extract-lines --n=1  'RMS DQ' 1 $testout | awk '{print $NF}' | sed s/DQ=//`
set dq2 = `extract-lines --n=1  'RMS DQ' 2 $testout | awk '{print $NF}' | sed s/DQ=//`
set dq3 = `extract-lines --n=1  'RMS DQ' 3 $testout | awk '{print $NF}' | sed s/DQ=//`
set dq4 = `extract-lines --n=1  'RMS DQ' 4 $testout | awk '{print $NF}' | sed s/DQ=//`

set f1 = -`extract-lines --n=2  'estatic                  eigval' 1 $testout | tail -1 | awk '{print $10}'`
set f2 = `extract-lines --n=2  'estatic                  eigval' 2 $testout | tail -1 | awk '{print $8}'`

compare_res_0 chk2c8a "RMS DQ original lattice, unrotated density" $dq1 $drmsqtol1 pass
chk2c8a:

compare_res_0 chk2c8b "RMS DQ  rotated lattice, unrotated density" $dq2 $drmsqtol1 pass
chk2c8b:

compare_res_0 chk2c8c "RMS DQ original lattice,   rotated density" $dq3 $drmsqtol1 pass
chk2c8c:

compare_res_0 chk2c8d "RMS DQ  rotated lattice,   rotated density" $dq4 $drmsqtol1 pass
chk2c8d:

set ftol = .01
compare_res chk2c8e "x component of rotated force match -z component of unrotated force" $f1 $f2 $ftol pass
chk2c8e:


endif


chk29:
if ($?clean) then
else if ($?pass) then
    echo "$space test 2 PASSED ($ext)"
    set jobspassed = ($jobspassed 2)
else
    echo "$space test 2 FAILED ($ext)"
    set failed = ($failed 2)
endif

chk2e:

echo $joblist | egrep '\b3\b' >/dev/null
if ($status) goto chk3e
set jobid = 3
cat <<EOF

         --- Test 3.  Check of miscellaneous special features, programs lmfa,lmf ---

EOF
if ($?quiet) then
else if ($ext == "srtio3") then
cat <<EOF

         The SrTiO3 test checks the --dos and --quit=dos switches as implemented in lmf.
         The test makes the total DOS from a potential generated by the Mattheis construction,
         invoking it by the switch
           --dos:npts=1001:rdm:ef0:window=-.5,.5
         This creates the DOS in a Questaal 2D array format (:rdm), with the VBM at 0 (:ef0)

EOF
else if ($ext == "bzt") then
cat <<EOF

         Writes the density of Ba3ZnTa2O9 to disk in the plane, using the --wden switch.
         Begins with a self-consistent calculation, followed by a pass
         where the charge density is written for a mesh of points in the xy plane.
         A second density (superposition of atomic densities) is also created
         If plotting packages are available the second is subtracted from the first
         and a contour plot drawn of the (self-consistent) - (free atomic) densities.

EOF
else if ($ext == "te") then
cat <<EOF
         Compares the PBE functional internal functionals comparing internally compiled Kieron Burke's
         easypbe.f to the PBE functional in libxc. The two construct some derivatives in different ways.

         This test should establish that a self-consistent density generated by one method
         is nearly self-consistent for the other, and that the total energy is nearly identical.
         Forces and Fermi levels are also compared.

         Note that a fixed core is used for both tests.  If the cores are calculated independently,
         slight differences in  differences in the total energy appear (see test 1).

EOF
else if ($ext == "au") then
cat <<EOF
         Sets up a test similar to that compared against the WIEN code in
         http://molmod.ugent.be/sites/default/files/deltadftcodes/supplmat/SupplMat-WIEN2k.pdf
         The lattice constant and bulk modulus are nearly identical to the WEIN code.
         To achieve a high level of agreement it is necessary to include local orbitals on the Au 6d,
         and treat the core fully relativistically.

         This test :
         (1) calculates the free atom with relativistic cores

         (2) demonstrates some features of controlling parameters in the basp file with
               lmfa --basp~incrlmx ctrl.au
             and
               lmf au -vnkabc=-1000 -vpwmode=0 --optbas:wbas

         (3) Compares an all smooth-Hankel basis to one with APWs added.

         Optimizing the basis set with repect to RSMH has little effect here, and is skipped.
         Switch --optbas:wbas merely rewrites the basp file, taking potential functions P
         from the rst file.

         If you carry out the following for the all LMTO basis set
           foreach dalat ( -0.15 -0.125 -0.1 -0.075 -0.05 -0.025 0 0.025 0.05 0.075 0.1 0.125 0.15 )
             rm -f mixm.au log.au
             lmf -vdalat=\$dalat -vpwmode=0 --rs=0 ctrl.au
           end
         or the following for the PMT basis set
           foreach dalat ( -0.15 -0.125 -0.1 -0.075 -0.05 -0.025 0 0.025 0.05 0.075 0.1 0.125 0.15 )
             rm -f mixm.au log.au
             lmf -vdalat=\$dalat -vpwmode=11 -vpwemax=34 --rs=0 ctrl.au
           end
         and determine the equilibrium lattice constant, you should find in both cases it
         deviates from the WIEN2K lattice constant (7.8579) by approximately 0.0005.

EOF
else if ($ext == "zbgan") then
cat <<EOF

         Writes the pseudodensity for the highest lying valence band in Zincblende GaN.
         at the Gamma point.

         Begins by making a self-consistent density.

         Next the charge density is written on a square lattice of points normal to
         the z axis, passing through (0,0,0.25).  This plane passes through
         (0.25,0.25,0.25) where the  N atom is located.

         The smoothed (pseudo) density is created for a single k-point very close to
         Gamma, for the 9th band: it is mostly a N py orbital as can be seen from
         a contour plot. If plotting packages are in your path contour plot is drawn.

         You might want to view some other bands.  Consider the following:
          *Make the square of the 1st conduction band w.f., in the xy plane passing through Ga:
           lmf ctrl.zbgan -vgmax=24 --wden~q=0,0,0.001~lst=10~core=2~ro=0,0,.00~l1=-1,1,1,41~l2=1,-1,1,41
          *Make the square of bands 2-6 (essentially the Ga d orbitals)
            lmf ctrl.zbgan -vgmax=24 --wden~q=0,0,0.001~lst=2:6~core=2~ro=0,0,.00~l1=-1,1,1,41~l2=1,-1,1,41
          *Make the square of band 8, passing through N (mostly a N px orbital)
           lmf ctrl.zbgan -vgmax=24 --wden~q=0,0,0.001~lst=8~core=2~ro=0,0,.25~l1=-1,1,1,41~l2=1,-1,1,41
          *Make the square of the (pseudo) 1st conduction band, passing through Cd at q=(1,0,0)
            lmf ctrl.zbgan -vgmax=24 --wden~q=1,0,0~lst=10~core=-1~ro=.0,.0,.25~l1=-1,1,1,41~l2=1,-1,1,41

          *You can make pictures of any of these with the following:
            cp smrho.zbgan  rho001 ; fplot -f fp/test/plot.rho.zbgan; gs fplot.ps

EOF
else if ($ext == "c") then
cat <<EOF

         The C test tests the code's implementation of homogeneous background mode
         It checks that program lmf generates the correct total energy and
         ionization potential for a single C atom.

         Note that:

        *The total energy of the neutral atom computed by lmf (-74.996 Ry) is very close
         to the free-atom energy computed by lmfa (-74.995 Ry), and that the free-atom
         potential is approximately self-consistent;

        *Also the SUM of total energy of the ionized system (-74.443) and the
         interaction of a point charge with a homogeneous background E^I,
         estimated from E^I=9/5R, with 4*pi*R^3/3=vol => R=6.20 and E^I = 0.290 Ry
         evaluates to -74.443+0.290 = -74.153 Ry, is close to the total energy of
         the positively charged free ion as computed by lmfa (-74.171 Ry).

EOF
else if ($ext == "fe") then
cat <<EOF

         The fe test tests scaling magnetic part of the exchange-correlation field.

EOF
else if ($ext == "kfese") then
cat <<EOF
         This job tests scaling the magnetic part of the XC functional for KFeSe
         (an Fe based superconductor) in a striped antiferromagnetic phase.

EOF
else if ($ext == "felz") then
cat <<EOF

         The felz test tests the code's implementation of fixed-spin moment method
         with and without spin-orbit coupling.

EOF
else if ($ext == "gas") then
cat <<EOF
         Test of the automatic band minimum/maximum finder, starting from a given k-point.

         This test makes a quadratic fit to the valence band Evbm at q, estimates
         the extremal q from the fit, and iterates until grad_q Evbm is near zero.

         Note: the potential is taken from a Mattheis construction of the density,
         so the bandgap is larger than self-consistent LDA gap.

         Other tests to consider:
         The following minimizes grad_q Evbm using a Broyden algorithm:
           band-edge -maxit=20 --bin -edge -r=.02 -dqmx=.05 -band=9 -q0=.1,.2,.3 gas
         Iterations are not particularly stable, but the solution eventually converges
         approximately to the Gamma point.

         The following finds a minimum in the conduction band at the L point, using a Broyden algorithm:
           band-edge -maxit=20 --bin -edge -r=.02 -dqmx=.05 -band=10 -q0=.1,.2,.3 gas
         Using instead the quadratic method a maximum along the Gamma-L line is found:
           band-edge -maxit=12 --bin -edge2=.7 -r=.02 -dqmx=.05 -band=10 -q0=.1,.2,.3 gas

         The following floats to local valence band maximum (or to CBM), finding it within the resolution of the cluster radius
           band-edge -cmd=lmf -maxit=30 --bin -floatmx -r=.05 -dqmx=.05 -band=9 -q0=.1,.2,.03 gas
           band-edge -cmd=lmf -maxit=30 --bin -floatmn -r=.05 -dqmx=.05 -band=10 -q0=.1,.2,.03 gas

         The following test evaluates the conduction band mass at Gamma.
         The result (m*/m = 0.036) is larger than the self-consistent LDA result. (m*/m = 0.020).
         Experimentally m*/m = 0.067.  The mass is too small because LDA underestamates the gap.
           band-edge -cmd=lmf --bin -mass -alat=10.66 -r=.001,.002,.003,.004 -band=10 -q0=0,0,0 gas

         No check is made to see whether the test executes correctly,
         but it compares the output against $testdir/out.lmf.vbm.$ext .

EOF

else if ($ext == "ni") then
cat <<EOF
         The Ni test demonstrates the fixed-spin moment method in two ways:

         First 3 iterations are run with an global external field Beff=8 mRy,
         but no constraint on the moment. (This is not the fixed-spin moment method,
         but closely related).  Were this calculation taken to self-consistency
         a magnetic moment of 0.702 would result.
         Compare it to the self-consistent moment in the absence of Beff, 0.603.

         Next a constraint is imposed, where Beff is determined each iteration to
         satisfy the constraint that the total moment is 0.7.
         This calculation is converged to self-consistency.
         Beff is found to be 7.7 mRy.

EOF
else if ($ext == "copt") then
cat <<EOF
         The copt test demonstrates the restart file editor.
         It splits the rst file into a charge-only part, a spin-only part,
         saves them (files rsta1,rsta2); then reads and adds them,
         and checks the output from this file against the original one.
         Note: run this test only after completing test 1.

EOF
else if ($ext == "coptso" || $ext == "co") then
cat <<EOF
         This test uses both the supercell maker and restart file editor to generate a supercell
         with site and restart files generated from the original cell.
         Note: run this test only after completing test 1.

EOF
else if ($ext == "mgo") then
cat <<EOF
         This test demonstrates the application of an external potential.

         Volume of augmentation spheres matches the cell volume.
         This test confirms that the following local smooth quantities
           rho*ves, rho*eps, rho*mu and rho*vext
         approximately match interstitial ones.

        *After self-consistency, three passes are made:

           1. No external potential

           2. External potential = constant 0.1 Ry

           3. External potential V_G = 0.1 Ry for G = (-1,1,1)2*pi/alat

           4. External potential V_G = 0.1 Ry for G = (-1,1,1)2*pi/alat

        *Passes 1 and 2 should be equivalent apart from:
         - (0.1 Ry shift) : Energy bands, Fermi level, band parameters enu and C
         - (0.2 Ry shift) : sum of core energies for Oxygen (it alone has foca=0)
         - (0.8 Ry shift) : single particle sum for 8 occupied electrons

         The test should confirm that the local true and smooth rho*vext are identical for
         a constant potential.  They cancel, resulting in a net rho*vexti = 0.8 Ry.
         For pass 3 these three terms are very similar not identical.

         It also should confirm that the change in single-particle sum less the
         interaction  rho*vext equal the change in the non self-consistent
         Harris-Foulkes energy ehf.

        *Pass 3 adds a potential with one Fourier component, and makes the
         density self-consistent.  Because of symmetry this Fourier component
         gets divided into 8 components.
         A check is made that the Harris-Foulkes and Kohn-Sham total energies match.

        *Pass 4 constructs a supercell of 5 MgO units along the [11-1] direction.
         A constant potential shift is added to half the cell and the
         density made self-consistent.

         The response and dielectric functions are analyzed in the output.

         The following test, which can be run after this test completes, performs a
         similar function except that the potential is made for a single Fourier component.
         The first command is for serial mode, the second for MPIK mode
           rdcmd -vncell=-5 -cat:TSTEPS --noerr ctrl.mgo
           rdcmd -vncell=-5 -cat:TSTEPM --noerr ctrl.mgo

EOF
else if ($ext == "al") then
cat <<EOF
         This test on elemental demonstrates how to compute band- and k- resolved DOS and ballistic conductivity.


EOF
else if ($ext == "fept") then
cat <<EOF
         This spin-polarized test demonstrates the following features:

         * Self-consistent calculations are made for three conditions.
           The same input file is used; the different cases are distinguished mostly through site files.
           1.  FM case, in the usual CuAu structure -- data taken from Pearson, 1817275 (file site.fept)
           2.  FM case, with the cell doubled along z (file site2.fept)
           3. AFM case, with the two Fe along z antiferromagnetically aligned

         * Site data are read from site files, along with species labels.
           This ability (starting with version 7.12) simplifies the assignment and counting of species.

         * The Fe 3p and Pt 5p are treated with extended local orbitals.
           Note: versions prior to 7.12 calculated the envelope shape from the spin-1 potential,
                 This results in some inconsistencies.
                 The old behavior can be mimicked by adding 40 to HAM_AUTOBAS_PFLOAT

         This test carries out the following steps.

         * Ferromagnetic FePt is calculated for the usual 2-atom cell.
           The restart file editor is used to yield a good starting point for a cell doubled along z.
           The self-consistent result shows moment and total energy almost exactly double the two-atom case

         * The restart file editor is used to globally exchange up- and down- spins.
           Band passes with the two spin configurations are compared and yield the same results, mutatis mutandis.

         * An extra species is inserted into the restart file using the restart file editor, to prepare for the AFM case.
           A band pass with a third species present yield the same results as the two-species case

         * A trial AFM density is constructed from the self-consistent FM density, with the restart file editor.
           The AFM case is iterated to self-consistency.
           The AFM state is calculated to be more stable than the FM one, in contradiction to experiment.

EOF
else
endif # materials
endif # quiet


set refout=$testdir/out.lmf.neutral.$ext testout=out.lmf.$ext refoutdos=$testdir/out.lmdos.$ext

if ($ext == "c") then
set testhom

else if ($ext == "au") then
set testfr
set refout=$testdir/out.lmf.fr.$ext testout=out.lmf.$ext

else if ($ext == "bzt" || $ext == "zbgan") then
set testwden
set refout=$testdir/out.lmf.wden.$ext testout=out.lmf.$ext

else if ($ext == "te") then
set testxc2
set refout=$testdir/out.lmf.pbex.$ext testout=out.lmf.$ext

else if ($ext == "ni" || $ext == "felz") then
set testfsm
set refout=$testdir/out.lmf.fsmom.$ext testout=out.lmf.$ext

else if ($ext == "mgo") then
set testvext
set refout=$testdir/out.lmf.vext.$ext testout=out.lmf.$ext testout2=out.lmf.scell.$ext refout2=$testdir/out.lmf.scell.$ext

else if ($ext == "gas") then
set testvbm
set refout=$testdir/out.lmf.vbm.$ext testout=out.lmf.$ext

else if ($ext == "fe" || $ext == "kfese") then
set testbxc
set refout=$testdir/out.lmf.bxc.$ext testout=out.lmf.$ext

else if ($ext == "srtio3") then
set testdos
set refout=$testdir/out.lmf.dos.$ext testout=out.lmf.$ext

else if ($ext == "al") then
set testdos
set refout=$testdir/out.lmf.dos.$ext testout=out.lmf.$ext

else if ($ext == "fept") then
set testfm
if ($?MPIK) set testmfm
set refout=$testdir/out.lmf.$ext testout=out.lmf.$ext

else if ($ext == "copt" || $ext == "coptso" || $ext == "co") then
set refout=$testdir/out.lmf.rsedit.$ext testout=out.lmf.$ext
set testrsedit
endif

if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk3e
endif
if (! $?testbxc && ! ($ext == "fept")) then
if (($?MPIK || $?MPI) && ($ext != "bzt") && ($ext != "mgo") && ($ext != "au") && ($ext != "c")) then
  echo "$space ... skipping this test ... no MPI-specific checks"
  goto chk3e
endif
endif
set pass

query chk31 chk3e 'run this test'
chk31:
# ... Look for executables
findcmd chk31a rdcmd "$path" "optional"
chk31a:
findcmd chk31b lmf "$path" "$topdir"
chk31b:
findcmd chk31c lmfa "$path" "optional"
chk31c:
findcmd chk31d bc "$path" "no"
chk31d:
findcmd chk31e vextract "$path" "no"
chk31e:

# If output already available, just run through checks
if ($?haveout) goto chk32

if ($?clean) then
  echo "$space cp $cplst ."
               cp $cplst .
  echo "$space rm -f {mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds,rsta1,rsta2,rsta3,sigm}.$ext $testout"
               rm -f {mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds,rsta1,rsta2,rsta3,sigm}.$ext $testout
  if (-e ctrl.$ext) then
    runrdcmd chk3e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk3e
endif

# ... Setup: remove existing files and copy new ones
if ($ext == "copt") then
  if (! -e ctrl.copt && ! -e rst.copt) then
     echo 'test.fp: run copt test 1 before running test 3'
     exit -1
  endif
else if ($ext == "coptso") then
  if (! -e ctrl.coptso || ! -e rsta.coptso) then
     echo 'test.fp: run coptso test 1 before running test 3'
     exit -1
  endif
else
echo "$space rm -f {mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds,rsta1,rsta2,rsta3,sigm}.$ext"
             rm -f {mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds,rsta1,rsta2,rsta3,sigm}.$ext
echo "$space cp $cplst ."
             cp $cplst .
if ($ext == "fept") then
echo "$space rm -f specialspec1"
             rm -f specialspec1
endif
endif

# ... Run test
if (! $?clean  && $?testfsm) then
  runrdcmd chk32 %11f $testout "-cat:TESTFSM --noerr ctrl.$ext"
else if (! $?clean  && $?testxc2) then
  runrdcmd chk32 %11f $testout "-cat:TESTXC2 --noerr ctrl.$ext"
else if (! $?clean && $?testwden && $?MPIK) then
  runrdcmd chk32 %11f $testout "-cat:TWDENMK --noerr ctrl.$ext"
else if (! $?clean && $?testvext && $?MPIK) then
  runrdcmd chk31i %11f $testout "-cat:TESTVXM --noerr ctrl.$ext"
  chk31i:
  if ($ext == "mgo") then
  echo "$space begin test 4 (MPI) ..."
  runrdcmd chk32 %11f $testout2 "-vncell=5 -cat:TSTEPM --noerr ctrl.$ext"
  endif
else if (! $?clean && $?testvext) then
  runrdcmd chk31j %11f $testout "-cat:TESTVXT --noerr ctrl.$ext"
  chk31j:
  if ($ext == "mgo") then
  echo "$space begin test 4 ..."
  runrdcmd chk32 %11f $testout2 "-vncell=5 -cat:TSTEPS --noerr ctrl.$ext"
  endif
else if (! $?clean  && $?testwden) then
  runrdcmd chk32 %11f $testout "-cat:TSTWDEN --noerr ctrl.$ext"
else if (! $?clean  && $?testfr && $?MPIK) then
  runrdcmd chk32 %11f $testout "-cat:TSTFRM --noerr ctrl.$ext"
else if (! $?clean  && $?testfr) then
  runrdcmd chk32 %11f $testout "-cat:TSTFR --noerr ctrl.$ext"
else if (! $?clean  && $?testvbm && $ext == gas) then
  findcmd chk31f fmin "$path" "optional"
  chk31f:
  if (! $?fmin) then
    echo "$space fmin missing ... skipping test"
    goto chk3e
  endif

  findcmd chk31g pfit "$path" "optional"
  chk31g:
  if (! $?pfit) then
    echo "$space pfit missing ... skipping test"
    goto chk3e
  endif

  lmfa gas > $testout
  echo "$space band-edge -cmd=lmf -maxit=20 --bin -edge2=.7 -r=.02 -dqmx=.05 -band=9 -q0=.1,.2,.3 gas > $testout"
               band-edge -cmd=lmf -maxit=20 --bin -edge2=.7 -r=.02 -dqmx=.05 -band=9 -q0=.1,.2,.3 gas > $testout
  call zdiffiles chk31h "CPU -1 $testout $refout"
  chk31h:
  goto chk3e
else if (! $?clean  && $?testrsedit) then
  runrdcmd chk32 %11f $testout "-cat:TESTRSE --noerr ctrl.$ext"
else if (! $?clean && $?MPIK && $?testbxc) then
  runrdcmd chk32 %11f $testout "-cat:TSTMBXC --noerr ctrl.$ext"
else if (! $?clean  && $?testbxc) then
  runrdcmd chk32 %11f $testout "-cat:TESTBXC --noerr ctrl.$ext"
else if (! $?clean  && $?testmfm) then
  runrdcmd chk32 %11f $testout "-cat:TESTMFM --noerr ctrl.$ext"
else if (! $?clean  && $?testfm) then
  runrdcmd chk32 %11f $testout "-cat:TESTFM --noerr ctrl.$ext"
else if (! $?clean  && $?testdos) then
  runrdcmd chk32 %11f $testout "-cat:TSTDOS --noerr ctrl.$ext"
else if (! $?clean && $?MPIK) then
  runrdcmd chk32 %11f $testout "-cat:TMPKLMF --noerr ctrl.$ext"
else if (! $?clean) then
  runrdcmd chk32 %11f $testout "-cat:TESTLMF --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk3e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk3e
endif
chk32:

# ... Special handling of Au (external potential)
if ($ext == "au") then

  set dq1f  = `extract-lines 'START LMF ' Exit 3 $testout | grep 'RMS DQ=' | $vextract . DQ | tail -1`
  set dq2f  = `extract-lines 'START LMF ' Exit 4 $testout | grep 'RMS DQ=' | $vextract . DQ | tail -1`
  set ef1   = `extract-lines 'START LMF ' Exit 3 $testout | grep 'Fermi energy' | sed 's/;//' | tail -1 | awk '{print $4}'`
  set ef2   = `extract-lines 'START LMF ' Exit 4 $testout | grep 'Fermi energy' | sed 's/;//' | tail -1 | awk '{print $4}'`
  set ehf1  = `extract-lines 'START LMF ' Exit 3 $testout | $vextract c ehf`
  set ehf2  = `extract-lines 'START LMF ' Exit 4 $testout | $vextract c ehf`
  set ehf1r = `extract-lines 'START LMF ' Exit 3 $refout | $vextract c ehf`
  set ehf2r = `extract-lines 'START LMF ' Exit 4 $refout | $vextract c ehf`


  if (! $?quiet) then

    call zdiffiles chk3j2 "CPU -1 $testout $refout"
    chk3j2:

#     call zdiffiles chk3j3 "CPU -1 $testout $refout2"
#     chk3j3:


    echo
    echo "$space self-consistent RMS DQ, all sm-Hankel basis                   = $dq1f"
    echo "$space self-consistent RMS DQ, PMT basis (APW emax=3)                = $dq2f"

    echo
    echo "$space self-consistent Fermi level, all sm-Hankel basis              = $ef1"
    echo "$space self-consistent Fermi level, PMT basis (APW emax=3)           = $ef2"
    set diff = `echo $ef1 $ef2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                                                    =  $diff"

    echo
    set diff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space self-consistent Harris-Foulkes energy, all sm-Hankel basis    = $ehf1" "change relative to reference:" $diff
    set diff = `echo $ehf2 $ehf2r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space self-consistent Harris-Foulkes energy, PMT basis (APW emax=3) = $ehf2" "change relative to reference:" $diff
    set diff = `echo $ehf1 $ehf2   | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space PMT and sm-Hankel total energy difference                     =  $diff"

    echo
  endif

  call qprint chk3j1 "$space ... automatic pass checks :"
  chk3j1:

  cmp -l basp.$ext $testdir/basp.$ext >& /dev/null
  set retval = $status
  echo -n "$space file basp.$ext identical to reference $testdir/basp.$ext ? ... "
  if ($retval == 0) then
   echo yes
  else
    echo -n "no ... differ by at most 1 char? ... "
    set n = `cmp -l basp.$ext $testdir/basp.$ext | wc -l`
    if ($n > 1) then
      echo no; unset pass
    else
      echo ok
    endif
  endif

  if ($?drmsqtol3 == 0) set drmsqtol3 = 4e-6
  if ($?ehkhftol == 0) set dehkhftol = 1e-3

  compare_res_0 chk3ja "RMS DQ, sm Hankel basis " $dq1f $drmsqtol3 pass
  chk3ja:

  compare_res_0 chk3jb "RMS DQ, PMT basis " $dq2f $drmsqtol3 pass
  chk3jb:

  set diff = `echo $ehf1 $ehf2   | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  compare_res_0 chk3jc "PMT-smH energy difference " $diff $dehkhftol pass
  chk3jc:

  echo
  goto chk3e0

endif

# ... Special handling of MgO (external potential)
if ($ext == "mgo") then

  set smlocves0 = `extract-lines 'START LMF' Exit 3 $testout | egrep '\bval\*ves'  | tail -1 | awk '{print $3}'`
  set smlocves1 = `extract-lines 'START LMF' Exit 4 $testout | egrep '\bval\*ves'  | tail -1 | awk '{print $3}'`
  set smlocves2 = `extract-lines 'START LMF' Exit 5 $testout | egrep '\bval\*ves'  | tail -1 | awk '{print $3}'`

  set smistves0 = `extract-lines 'START LMF' Exit 3 $testout | egrep 'rhoval\*ves' | tail -1 | awk '{print $2}'`
  set smistves1 = `extract-lines 'START LMF' Exit 4 $testout | egrep 'rhoval\*ves' | tail -1 | awk '{print $2}'`
  set smistves2 = `extract-lines 'START LMF' Exit 5 $testout | egrep 'rhoval\*ves' | tail -1 | awk '{print $2}'`

  set smrhoeps0 = `extract-lines 'START LMF' Exit 3 $testout | egrep 'rhoeps:'     | tail -1 | awk '{print $3}'`
  set smrhoeps1 = `extract-lines 'START LMF' Exit 4 $testout | egrep 'rhoeps:'     | tail -1 | awk '{print $3}'`
  set smrhoeps2 = `extract-lines 'START LMF' Exit 5 $testout | egrep 'rhoeps:'     | tail -1 | awk '{print $3}'`

  set smisteps0 = `extract-lines 'START LMF' Exit 3 $testout | egrep 'rho\*exc'    | tail -1 | awk '{print $2}'`
  set smisteps1 = `extract-lines 'START LMF' Exit 4 $testout | egrep 'rho\*exc'    | tail -1 | awk '{print $2}'`
  set smisteps2 = `extract-lines 'START LMF' Exit 5 $testout | egrep 'rho\*exc'    | tail -1 | awk '{print $2}'`

  set smrhomu0  = `extract-lines 'START LMF' Exit 3 $testout | egrep 'rhomu:'      | tail -1 | awk '{print $3}'`
  set smrhomu1  = `extract-lines 'START LMF' Exit 4 $testout | egrep 'rhomu:'      | tail -1 | awk '{print $3}'`
  set smrhomu2  = `extract-lines 'START LMF' Exit 5 $testout | egrep 'rhomu:'      | tail -1 | awk '{print $3}'`

  set smistmu0  = `extract-lines 'START LMF' Exit 3 $testout | egrep 'rho\*vxc'    | tail -1 | awk '{print $2}'`
  set smistmu1  = `extract-lines 'START LMF' Exit 4 $testout | egrep 'rho\*vxc'    | tail -1 | awk '{print $2}'`
  set smistmu2  = `extract-lines 'START LMF' Exit 5 $testout | egrep 'rho\*vxc'    | tail -1 | awk '{print $2}'`

  set sumev0    = `extract-lines 'START LMF' Exit 3 $testout | egrep 'sumev='      | egrep 'srhov:' | head -1 | awk '{print $6}'`
  set sumev1    = `extract-lines 'START LMF' Exit 4 $testout | egrep 'sumev='      | egrep 'srhov:' | head -1 | awk '{print $6}'`
  set sumev2    = `extract-lines 'START LMF' Exit 5 $testout | egrep 'sumev='      | egrep 'srhov:' | head -1 | awk '{print $6}'`

  set sumec0    = `extract-lines 'START LMF' Exit 3 $testout | egrep 'sumec='      | tail -1 | awk '{print $2}'`
  set sumec1    = `extract-lines 'START LMF' Exit 4 $testout | egrep 'sumec='      | tail -1 | awk '{print $2}'`
  set sumec2    = `extract-lines 'START LMF' Exit 5 $testout | egrep 'sumec='      | tail -1 | awk '{print $2}'`

  set eferm0    = `extract-lines 'START LMF' Exit 3 $testout | egrep 'Fermi '      | tail -1 | awk '{print $4}' | sed 's/;//'`
  set eferm1    = `extract-lines 'START LMF' Exit 4 $testout | egrep 'Fermi '      | tail -1 | awk '{print $4}' | sed 's/;//'`
  set eferm2    = `extract-lines 'START LMF' Exit 5 $testout | egrep 'Fermi '      | tail -1 | awk '{print $4}' | sed 's/;//'`

  set ehf0      = `extract-lines 'START LMF' Exit 3 $testout | egrep 'nit=' | head -1 | $vextract . ehf `
  set ehf1      = `extract-lines 'START LMF' Exit 4 $testout | egrep 'nit=' | head -1 | $vextract . ehf `
  set ehf2      = `extract-lines 'START LMF' Exit 5 $testout | egrep 'nit=' | head -1 | $vextract . ehf `
  set ehf2x     = `extract-lines 'START LMF' Exit 5 $testout | egrep 'nit=' | tail -1 | $vextract . ehf `

  set ehk0      = `extract-lines 'START LMF' Exit 3 $testout | egrep 'nit=' | head -1 | $vextract . ehk `
  set ehk1      = `extract-lines 'START LMF' Exit 4 $testout | egrep 'nit=' | head -1 | $vextract . ehk `
  set ehk2      = `extract-lines 'START LMF' Exit 5 $testout | egrep 'nit=' | head -1 | $vextract . ehk `
  set ehk2x     = `extract-lines 'START LMF' Exit 5 $testout | egrep 'nit=' | tail -1 | $vextract . ehk `

  set smrhovxt1c = `extract-lines 'START LMF' Exit 4 $testout | egrep 'valence rho' | head -1 | awk '{print $3}'`
  set smrhovxt2c = `extract-lines 'START LMF' Exit 4 $testout | egrep 'valence rho' | head -1 | awk '{print $4}'`
  set smrhovxtic = `extract-lines 'START LMF' Exit 4 $testout | egrep 'valence rho' | head -1 | awk '{print $5}'`
  set smrhovxtc  = `extract-lines 'START LMF' Exit 4 $testout | egrep 'valence rho' | head -1 | awk '{print $6}'`
  set delsevc    = `echo $sumev1 - $smrhovxtc | bc`

  set smrhovxt1g = `extract-lines 'START LMF' Exit 5 $testout | egrep 'valence rho' | head -1 | awk '{print $3}'`
  set smrhovxt2g = `extract-lines 'START LMF' Exit 5 $testout | egrep 'valence rho' | head -1 | awk '{print $4}'`
  set smrhovxtig = `extract-lines 'START LMF' Exit 5 $testout | egrep 'valence rho' | head -1 | awk '{print $5}'`
  set smrhovxtg  = `extract-lines 'START LMF' Exit 5 $testout | egrep 'valence rho' | head -1 | awk '{print $6}'`
  set delsevg    = `echo "$sumev2 - $smrhovxtg" | bc`

  set ehfs      = `extract-lines 'START LMF' Exit 5 $testout2 | vextract c ehf`
  set ehks      = `extract-lines 'START LMF' Exit 5 $testout2 | vextract c ehk`
  set epsinv    = `extract-lines 'START LMF' Exit 5 $testout2 | egrep -A 1 eps.- | tail -1 | awk '{print $10}'`

  if (! $?quiet) then

  call zdiffiles chk3d2 "CPU -1 $testout $refout"
  chk3d2:
  if ($ext == "mgo") then
    call zdiffiles chk3d3 "CPU -1 $testout2 $refout2"
    chk3d3:
  endif

  echo "$space Integrated quantities                vext=0       vext=cnst   vext=e^iG.r"
  echo "$space smooth local int rho*ves          = $smlocves0   $smlocves1   $smlocves2"
  echo "$space interstitial int rho*ves          = $smistves0   $smistves1   $smistves2"
  echo

  echo "$space smooth local int rho*eps          =  $smrhoeps0    $smrhoeps1    $smrhoeps2"
  echo "$space interstitial int rho*eps          =  $smisteps0    $smisteps1    $smisteps2"
  echo

  echo "$space smooth local int rho*mu           =  $smrhomu0    $smrhomu1    $smrhomu2"
  echo "$space interstitial int rho*mu           =  $smistmu0    $smistmu1    $smistmu2"
  echo

  echo "$space (true  valence rho)*vext          =                $smrhovxt1c    $smrhovxt1g"
  echo "$space (sm-loc valence rho)*vext         =                $smrhovxt2c    $smrhovxt2g"
  echo "$space (istl   valence rho)*vext         =                $smrhovxtic    $smrhovxtig"
  echo

  echo "$space Fermi energy                      =  $eferm0    $eferm1    $eferm2"
  echo "$space foca=0 sum of core evals          = $sumec0   $sumec1   $sumec2"
  echo "$space single particle sum (1st iter)    =  $sumev0    $sumev1    $sumev2"
  echo "$space rhoval*vext                       =                $smrhovxtc    $smrhovxtg"

  echo "$space difference                        =  $sumev0    $delsevc    $delsevg"
  echo "$space relative to vext=0                =               " \
             `echo $delsevc $sumev0 | awk '{printf "%9.6f", $1-$2}'` "  " \
             `echo $delsevg $sumev0 | awk '{printf "%9.6f\n", $1-$2}'`

  echo "$space Harris-Foulkes energy             =  $ehf0    $ehf1    $ehf2  (1st it after s.c.)"
  echo "$space relative to vext=0                =               " \
             `echo $ehf1 $ehf0 | awk '{printf "%9.6f", $1-$2}'` "  " \
             `echo $ehf2 $ehf0 | awk '{printf "%9.6f\n", $1-$2}'`
  echo "$space Kohn-Sham energy                  =  $ehk0    $ehk1    $ehk2  (1st it after s.c.)"
  echo "$space HF-KS energy difference           = "\
             `echo $ehf0 $ehk0 | awk '{printf "%9.6f", $1-$2}'` "  "  \
             `echo $ehf1 $ehk1 | awk '{printf "%9.6f", $1-$2}'` "  "  \
             `echo $ehf2 $ehk2 | awk '{printf "%9.6f", $1-$2}'` " (1st it after s.c.)"
  echo "$space Harris-Foulkes energy, s.c.       =  $ehf0    $ehf1    $ehf2x"
  echo "$space HF-KS energy difference           = "\
             `echo $ehf0 $ehk0 | awk '{printf "%9.6f", $1-$2}'` "  "  \
             `echo $ehf1 $ehk1 | awk '{printf "%9.6f", $1-$2}'` "  "  \
             `echo $ehf2x $ehk2x | awk '{printf "%9.6f", $1-$2}'`

  echo
  echo "$space Supercell results:"
  echo "$space Harris-Foulkes energy, s.c.       = " `echo $ehfs | awk '{printf "%9.6f", $1}'` " per MgO atom" `echo $ehfs | awk '{printf "%9.6f", $1/5}'`
  echo "$space Kohn-Sham energy, s.c.            = " `echo $ehks | awk '{printf "%9.6f", $1}'` " per MgO atom" `echo $ehks | awk '{printf "%9.6f", $1/5}'`
  echo "$space HF-KS energy difference           = " `echo $ehfs $ehks | awk '{printf "%9.6f", $1-$2}'` "  "
  echo "$space eps^-1 for smallest k vector      =  " `echo $epsinv | awk '{printf "%9.6f", $1}'`

  endif

  echo
  echo "$space ... automatic pass checks :"
  if ($?detol == 0) set detol = 2e-6
  if ($?detols == 0) set detols = 1e-5

  set ref       = `extract-lines 'START LMF' Exit 3 $refout  | egrep 'Fermi '      | tail -1 | awk '{print $4}' | sed 's/;//'`
  compare_res chk3h2 "Fermi energy (vext=0)    " $eferm0 $ref $detol pass
  chk3h2:

  set ref       = `extract-lines 'START LMF' Exit 4 $refout  | egrep 'Fermi '      | tail -1 | awk '{print $4}' | sed 's/;//'`
  compare_res chk3h3 "Fermi energy (vext=const)" $eferm1 $ref $detol pass
  chk3h3:

  set ref       = `extract-lines 'START LMF' Exit 5 $refout  | egrep 'Fermi '      | tail -1 | awk '{print $4}' | sed 's/;//'`
  compare_res chk3h4 "Fermi energy (vext=e^iGr)" $eferm2 $ref $detol pass
  chk3h4:

  set ref       = `extract-lines 'START LMF' Exit 3 $refout  | $vextract x ehf `
  compare_res chk3h5 "Total energy (vext=0)    " $ehf0 $ref $detol pass
  chk3h5:

  set ref       = `extract-lines 'START LMF' Exit 4 $refout  | $vextract x ehf `
  compare_res chk3h6 "Total energy (vext=const)" $ehf1 $ref $detol pass
  chk3h6:

  set ref       = `extract-lines 'START LMF' Exit 5 $refout | egrep 'nit=' | tail -1 | $vextract . ehf `
  compare_res chk3h7 "Total energy (vext=e^iGr)" $ehf2x $ref $detol pass
  chk3h7:

  set ref       = `extract-lines 'START LMF' Exit 5 $refout | egrep 'nit=' | head -1 | $vextract . ehf `
  compare_res chk3h8 "Same, but first iteration"      $ehf2 $ref $detol pass
  chk3h8:

  set ref       = `echo $delsevg $sumev0 $ehf2 $ehf0 | awk '{printf "%9.6f\n", $1-$2}'`
  set ref2      = `echo $delsevg $sumev0 $ehf2 $ehf0 | awk '{printf "%9.6f\n", $3-$4}'`
  set ref3      = `echo $delsevg $sumev0 $ehf2 $ehf0 | awk '{printf "%9.6f\n", $1-$2-$3+$4}'`
  compare_res_0 chk3h9 "Difference ehf ($ref) - Difference sumev ($ref2) = " $ref3 $detol pass
  chk3h9:

  echo
  echo "$space Supercell results:"
  set ref       = `extract-lines 'START LMF' Exit 5 $refout2 | vextract c ehf`
  compare_res chk3ha "H-F total energy" $ehfs $ref $detol pass
  chk3ha:
  set ref       = `echo $ehfs $ehks | awk '{printf "%9.6f", $1-$2}'`
  compare_res_0 chk3hb "H-F H-K energy difference" $ref $detols pass
  chk3hb:


  goto chk3e0

# ... Special handling of BZT and zbgan (density makers)
else if ($ext == "bzt" || $ext == "zbgan") then

  if (! $?quiet) then
    call zdiffiles chk3g0 "CPU -1 $testout $refout"
    chk3g0:
  endif

  echo
  echo "$space ... automatic pass checks :"

# Compare smrho to reference
  set ndig = 4
  set retval = `cmp -l smrho.$ext $testdir/smrho.$ext |& grep -v EOF | wc | awk '{printf "%d", $1}'`
  echo -n "$space ... files smrho.$ext and $testdir/smrho.$ext identical? ... "
  if ($retval == 0) then
    echo yes
  else
    echo -n "no ... to 5 digits? ... "
    set ndig = 5
    call zcmpnfiles chk3g1 "$ndig smrho.$ext $testdir/smrho.$ext"
    chk3g1:
    if ($retval == 0) then
      echo yes
    else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
      echo ok "($retval difference(s) of $ncharfile)"
    else
      echo no "($retval difference(s) remaining of $ncharfile)"
      unset pass
    endif
  endif

  if ($?have_mc && $?have_fplot && $ext == "bzt" && $?slow) then
  echo "$space ... generate [self-consistent - free atom] density and shift origin"
  echo "$space mcx smrho.$ext atrho.$ext -- -roll:0,19 > rho001"
               mcx smrho.$ext atrho.$ext -- -roll:0,19 > rho001
  endif
  if ($?have_mc && $?have_fplot && $ext == "zbgan" && $?slow) then
  echo "$space ... generate sm density and shift origin"
#   echo "$space mcx smrho.$ext -roll:20,20 > rho001"
#                mcx smrho.$ext -roll:20,20 > rho001
  echo "$space cp smrho.$ext rho001"
               cp smrho.$ext rho001
  endif
  if ($?have_mc && $?have_fplot && $?slow) then
  query chk3g2 chk3e0 'make contour plot of density'
  chk3g2:
  echo "$space fplot -disp -pr10 -f $topdir/fp/test/plot.rho.$ext"
               fplot -disp -pr10 -f $topdir/fp/test/plot.rho.$ext
  endif
  goto chk3e0

# ... Special handling of Te (PBE functional)
else if ($ext == "te" ) then

  set dq1f  = `extract-lines 'START LMF' Exit 2 $testout | grep 'RMS DQ=' | $vextract . DQ | tail -1`
# set dqr1  = `extract-lines 'START LMF' Exit 2 $refout | grep 'RMS DQ=' | $vextract . DQ | tail -1`
  set dq2i  = `extract-lines 'START LMF' Exit 3 $testout | grep 'RMS DQ=' | $vextract . DQ | head -1`
  set dq2f  = `extract-lines 'START LMF' Exit 3 $testout | grep 'RMS DQ=' | $vextract . DQ | tail -1`
  set ehf1  = `extract-lines 'START LMF' Exit 2 $testout | $vextract c ehf`
  set ehk1  = `extract-lines 'START LMF' Exit 2 $testout | $vextract c ehk`
  set ehf2i = `extract-lines 'START LMF' Exit 3 $testout | $vextract i ehf | head -1`
  set ehf2  = `extract-lines 'START LMF' Exit 3 $testout | $vextract c ehf`
  set ehk2  = `extract-lines 'START LMF' Exit 3 $testout | $vextract c ehk`

  set ef1   = `extract-lines 'START LMF' Exit 2 $testout  | grep 'Fermi energy' | sed 's/;//' | tail -1 | awk '{print $3}'`
  set ef2   = `extract-lines 'START LMF' Exit 3 $testout  | grep 'Fermi energy' | sed 's/;//' | tail -1 | awk '{print $3}'`

  set fmax1  = `extract-lines 'START LMF' Exit 2 $testout  | grep 'Maximum Harris' | tail -1 | awk '{print $5}'`
  set fmax2i = `extract-lines 'START LMF' Exit 3 $testout  | grep 'Maximum Harris' | head -1 | awk '{print $5}'`
  set fmax2  = `extract-lines 'START LMF' Exit 3 $testout  | grep 'Maximum Harris' | tail -1 | awk '{print $5}'`


  if (! $?quiet) then

    call zdiffiles chk3i2 "CPU -1 $testout $refout"
    chk3i2:

    echo
    echo "$space RMS DQ, easypbe, self-consistent     = $dq1f"
    echo "$space RMS DQ, libxc, first iteration       = $dq2i"
    echo "$space RMS DQ, libxc, self-consistent       = $dq2f"

    echo
    echo "$space HF energy, easybe, self-consistent   = $ehf1"
    echo "$space HF energy, libxc, 1st iteration      = $ehf2i"
    set ediff = `echo $ehf1 $ehf2i  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                           = $ediff"
    echo "$space HF energy, libxc, self-consistent    = $ehf2"
    set ediff = `echo $ehf1 $ehf2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference relative to easypbe       = $ediff"

    echo
    echo "$space HK energy, easybe, self-consistent   = $ehk1"
    set ediff = `echo $ehf1 $ehk1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space HK-HF energy, easybe                 = $ediff"
    echo "$space HK energy, libxc, self-consistent    = $ehk2"
    set ediff = `echo $ehf2 $ehk2 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space HK-HF energy, libxc                  = $ediff"

    echo
    echo "$space Fermi level, easybe, self-consistent  = $ef1"
    echo "$space Fermi level, libxc, self-consistent   = $ef2"

    echo
    echo "$space Max force, easybe, self-consistent  = $fmax1"
    echo "$space Max force, libxc, 1st iteration     = $fmax2i"
    echo "$space Max force, libxc, self-consistent   = $fmax2"

    echo

  endif

  call qprint chk3i1 "$space ... automatic pass checks :"
  chk3i1:

  if ($?drmsqtol3i == 0) set drmsqtol3i = 1e-5
  if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-6
  if ($?ehkhftol == 0) set ehkhftol = 4e-6
  if ($?eftol == 0) set eftol = 2e-4

  compare_res_0 chk3ia "RMS DQ, easypbe " $dq1f 1d-6 pass
  chk3ia:

  compare_res_0 chk3ib "RMS DQ, libxc 1st iteration " $dq2i $drmsqtol3i pass
  chk3ib:

  compare_res_0 chk3ic "RMS DQ, libxc last iteration " $dq2f 1d-6 pass
  chk3ic:

  set ediff = `echo $ehf1 $ehf2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  compare_res_0 chk3id "HF energy, easypbe-libxc " $ediff $drmsqtol3i pass
  chk3id:

  set ediff = `echo $ehf2 $ehk2 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  compare_res_0 chk3ie "HK-HF energy, libxc" $ediff $ehkhftol pass
  chk3ie:

  set ediff = `echo $ef1 $ef2 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  compare_res_0 chk3if "Difference in Fermi levels" $ediff $eftol pass
  chk3if:

  goto chk3e0

# ... Special handling of SrTiO3
else if ($ext == "srtio3") then
  if (! $?quiet) then
    call zdiffiles chk3a0 "CPU -1 $testout $refout"
    chk3a0:
  endif

# Compare dos to reference
  set ndig = 4
  set retval = `cmp -l dos.$ext $testdir/dos.$ext |& grep -v EOF | wc | awk '{printf "%d", $1}'`
  echo -n "$space ... files dos.$ext and $testdir/dos.$ext equivalent to $ndig digits? ... "
  if ($retval == 0) then
    echo yes
  else
    echo -n "no ... to 3 digits? ... "
    set ndig = 3
    call zcmpnfiles chk3a2 "$ndig dos.$ext $testdir/dos.$ext"
    chk3a2:
    if ($retval == 0) then
      echo yes
    else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
      echo ok "($retval difference(s) of $ncharfile)"
    else
      echo no "($retval difference(s) remaining of $ncharfile)"
      unset pass
    endif
  endif

  goto chk3e0
# end of srtio3 checks

# ... Special handling of fept
else if ($ext == "fept") then
   set n1 = `grep -c 'FM, CuAu structure' $testout`
   set n2 = `grep -c 'FM, CuAu structure' $refout`

   set ehf1  = `extract-lines 'FM, CuAu structure' Exit $n1 $testout | $vextract c ehf `
   set mmom1 = `extract-lines 'FM, CuAu structure' Exit $n1 $testout | $vextract c mmom`
   set ehfr1  = `extract-lines 'FM, CuAu structure' Exit $n2 $refout | $vextract c ehf `
   set mmomr1 = `extract-lines 'FM, CuAu structure' Exit $n2 $refout | $vextract c mmom`

   set n1 = `grep -c 'cell doubled' $testout`
   set n2 = `grep -c 'cell doubled' $refout`

   set dq1  = `extract-lines 'FM, cell doubled along z' Exit $n1 $testout | grep DQ | sed -n 1,1p | $vextract . DQ`
   set ehf2  = `extract-lines 'FM, cell doubled along z' Exit $n1 $testout | $vextract c ehf `
   set mmom2 = `extract-lines 'FM, cell doubled along z' Exit $n1 $testout | $vextract c mmom`

   set dqr1 = `extract-lines 'FM, cell doubled along z' Exit $n2 $refout | grep DQ | sed -n 1,1p | $vextract . DQ`
   set ehfr2  = `extract-lines 'FM, cell doubled along z' Exit $n2 $refout | $vextract c ehf `
   set mmomr2 = `extract-lines 'FM, cell doubled along z' Exit $n2 $refout | $vextract c mmom`

   set n1 = `grep -c 'Band pass all spin up' $testout`
   set n2 = `grep -c 'Band pass all spin up' $refout`

   set dqfu = `extract-lines 'Band pass all spin up' '^x' $n1 $testout | grep DQ | sed -n 1,1p | $vextract . DQ`
   set momu = `extract-lines 'Band pass all spin up' '^x' $n1 $testout | grep ^x | sed -n 1,1p | $vextract . mmom`
   set ehfu = `extract-lines 'Band pass all spin up' '^x' $n1 $testout | grep ^x | sed -n 1,1p | $vextract . ehf`

   set dqfd = `extract-lines 'Band pass all spin down' '^x' $n1 $testout | grep DQ | sed -n 1,1p | $vextract . DQ`
   set momd = `extract-lines 'Band pass all spin down' '^x' $n1 $testout | grep ^x | sed -n 1,1p | $vextract . mmom`
   set ehfd = `extract-lines 'Band pass all spin down' '^x' $n1 $testout | grep ^x | sed -n 1,1p | $vextract . ehf`
   set ehfdr = `extract-lines 'Band pass all spin down' '^x' $n2 $refout | grep ^x | sed -n 1,1p | $vextract . ehf`

   set n1 = `grep -c 'FM band pass with extra species' $testout`
   set n2 = `grep -c 'FM band pass with extra species' $refout`

   set dqf3 = `extract-lines 'FM band pass with extra species' '^x' $n1 $testout | grep DQ | sed -n 1,1p | $vextract . DQ`
   set mom3 = `extract-lines 'FM band pass with extra species' '^x' $n1 $testout | grep ^x | sed -n 1,1p | $vextract . mmom`
   set ehf3 = `extract-lines 'FM band pass with extra species' '^x' $n1 $testout | grep ^x | sed -n 1,1p | $vextract . ehf`

   set n1 = `grep -c 'Iterate AFM case to self-consistency' $testout`
   set n2 = `grep -c 'Iterate AFM case to self-consistency' $refout`

   set dqa  = `extract-lines 'Iterate AFM case to self-consistency' Exit $n1 $testout | grep DQ | tail -1 | $vextract . DQ`
   set ehfa = `extract-lines 'Iterate AFM case to self-consistency' Exit $n1 $testout | grep ^c | $vextract c ehf`
   set ehfar = `extract-lines 'Iterate AFM case to self-consistency' Exit $n2 $refout | grep ^c | $vextract c ehf`

  if ($?add0) then
    echo -n "         ..." ; $add0 $testout
  else if ($?poszer) then
    echo -n "         ..." ; $poszer $testout
  endif

  if (! $?quiet) then

  call zdiffiles chk3d0 "CPU -1 $testout $refout"
  chk3d0:

  echo "$space Harris energy, original unit cell = $ehf1"
  echo "$space corresponding reference energy    = $ehfr1"
  set ediff = `echo $ehf1 $ehfr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space mag. moment,   original unit cell = $mmom1"
  echo "$space corresponding reference moment    = $mmomr1"
  set ediff = `echo $mmom1 $mmomr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space RMS DQ first iteration supercell  = $dq1"
  echo "$space corresponding reference           = $dqr1"
  set ediff = `echo $dq1 $dqr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"

  echo "$space Harris energy, supercell          = $ehf2"
  echo "$space corresponding reference energy    = $ehfr2"
  set ediff = `echo $ehf2 $ehfr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space mag. moment, supercell            = $mmom2"
  echo "$space corresponding reference moment    = $mmomr2"
  set ediff = `echo $mmom2 $mmomr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo -n "$space ehf(4 atoms)/2 - ehf(2 atoms)     = "; echo $ehf1 $ehf2   | awk '{printf "%9.6f\n", $2/2-$1}'
  echo -n "$space mom(4 atoms)/2 - mom(2 atoms)     = "; echo $mmom1 $mmom2 | awk '{printf "%9.6f\n", $2/2-$1}'
  echo ' '

  echo "$space RMS DQ FM spin up                 = $dqfu"
  echo "$space RMS DQ FM spin down               = $dqfd"
  echo "$space RMS DQ FM 3 species               = $dqf3"
  echo ' '

  echo "$space Harris energy spin up             = $ehfu"
  echo "$space Harris energy spin down           = $ehfd"
  echo "$space Harris energy 3 species           = $ehf3"

  set ediff = `echo $ehfu $ehfd  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space (up,down) difference              = $ediff"
  set ediff = `echo $ehfu $ehf3  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space (up,3) difference                 = $ediff"
  echo ' '

  echo "$space magnetic moment spin up           = $momu"
  echo "$space magnetic moment spin down         = $momd"
  echo "$space magnetic moment 3 species         = $mom3"
  set ediff = `echo -$momu $momd  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space (up,down) difference              = $ediff"
  set ediff = `echo $momu $mom3  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space (up,3) difference                 = $ediff"
  echo ' '

  echo "$space RMS DQ, self-consistent AFM       = $dqa"
  echo "$space ehf, self-consistent AFM          = $ehfa"
  echo "$space ehf, self-consistent FM           = $ehf2"
  set ediff = `echo $ehf2 $ehfa  | awk '{{k=($2-$1)} print k}'`
  echo "$space (AFM-FM)energy difference         = $ediff"
  echo ' '

  endif

  call qprint chk3f1 "$space ... automatic pass checks :"
  chk3f1:

if ($?defatol3 == 0) set defatol3 = 2e-6
if ($?dehf1tol3 == 0) set dehf1tol3 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol3 == 0) set dmom1tol3 = 1e-4
if ($?dmomntol3 == 0) set dmomntol3 = 1e-4
if ($?dfmax1tol3 == 0) set dfmax1tol3 = 0.1
if ($?dfmaxntol3 == 0) set dfmaxntol3 = 0.1
if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-4
echo " "

compare_res chk3f2 "ehf original cell" $ehf1 $ehfr1 $dehf1tol3 pass
chk3f2:

compare_res chk3f3 "ehf supercell    " $ehf2 $ehfr2 $dehf1tol3 pass
chk3f3:

compare_res chk3f4 "mom original cell" $mmom1 $mmomr1 $dmomntol3 pass
chk3f4:

compare_res chk3f5 "mom supercell    " $mmom2 $mmomr2 $dmomntol3 pass
chk3f5:

compare_res chk3f6 "ehf FM spin down " $ehfd $ehfdr $dehf1tol3 pass
chk3f6:

compare_res chk3f7 "ehf AFM          " $ehfa $ehfar $dehf1tol3 pass
chk3f7:

goto chk3e0
# end of fept checks
endif

# ... Special handling of al
if ($ext == "al" ) then

  set bbot    = `grep -A 1 'bndfp:  kpt 1 ' $testout | tail -1 | awk '{print $1}'`
  set bbotref = `grep -A 1 'bndfp:  kpt 1 ' $refout | tail -1 | awk '{print $1}'`

  set ef = `grep 'Fermi energy' $testout | tail -1 | awk '{print $4}' | sed 's/;//'`
  set efref = `grep 'Fermi energy' $refout | tail -1 | awk '{print $4}' | sed 's/;//'`

  set w = `echo $ef-"($bbot)" | bc`

  echo
  echo "$space ... calculate <|v|> at band bottom + 1/10 bandwidth"
  echo "$space" \
  "lmdos al -vnk=72 --dos~bands=1:3~rdm~window=-.75,-.65~npts=101~mode=3 --kres=$bbot+$w/10 > out.lmdos.al"
   lmdos al -vnk=72 --dos~bands=1:3~rdm~window=-.75,-.65~npts=101~mode=3 --kres=$bbot+$w/10 > out.lmdos.al

  set vbot = `grep 'at Ef' out.lmdos.al | head -1 | awk '{print $10}'`
  set vbotr = `grep 'at Ef' $refoutdos | head -1 | awk '{print $10}'`

  echo
  echo "$space ... calculate x component of ballistic conductivity Ef"
  echo "$space" \
  "lmdos al -vnk=72 --dos~bands=1:3~rdm~window=-.05,.05~npts=101~mode=1~vec=1,0,0  --kres=$ef >> out.lmdos.al"
   lmdos al -vnk=72 --dos~bands=1:3~rdm~window=-.05,.05~npts=101~mode=1~vec=1,0,0  --kres=$ef >> out.lmdos.al

  if (! $?quiet) then
  echo
  echo "$space Bottom of the band               = $bbot"
  echo "$space Reference value                  = $bbotref"
  set ediff = `echo $bbot $bbotref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  echo
  echo "$space Fermi energy                     = $ef"
  echo "$space Reference value                  = $efref"
  set ediff = `echo $ef $efref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  echo
  echo "$space bandwidth                        = $w"
  echo "$space Free electron value              = 0.866"

  echo
  echo "$space velocity at bot + bandwidth/10   = $vbot x 10^6 m/s"
  echo "$space reference                        = $vbotr x 10^6 m/s"
  echo "$space free electron value              = 0.64  x 10^6 m/s"
  echo


  call zdiffiles chk3k1 "CPU -1 $testout $refout"
  chk3k1:

  call zdiffilesx chk3k2 "CPU -1 out.lmdos.al $refoutdos"
  chk3k2:

# if (! $?mcx) then
#   echo "$space mcx not in path ... no checks on dos"
#   goto chk3e0
# endif
# echo
# echo "$space make |v| from dosq.al.1 dosq.al.2 dosq.al.3, separating bands 2,3"
# echo "$space" \
# "mcx [ k=1,2,3 dosq.al.'{k}' -inc x5==2 -a 'd{k}'  'd{k}' -coll 6 -p -xe ] -+ -+ -e1 'sqrt(x1)' d1 -coll 1,2,3 -tog -ccat > dos2.al"
# echo "$space" \
# "mcx [ k=1,2,3 dosq.al.'{k}' -inc x5==3 -a 'd{k}'  'd{k}' -coll 6 -p -xe ] -+ -+ -e1 'sqrt(x1)' d1 -coll 1,2,3 -tog -ccat > dos3.al"

  endif # ! quiet

  echo
  echo "$space ... automatic pass checks :"
  set detol = 1e-4
  compare_res chk3e1 "Band bottom  " $bbot $bbotref $detol pass
  chk3e1:
  set detol = 1e-6
  compare_res chk3e2 "Fermi energy  " $ef $efref $detol pass
  chk3e2:

  if (! $?gmaxdif) then
  else if (! $?mcx) then
    echo "$space ... mcx not installed ... no check on stdout"
  else
    if (! $?stdotol3) set stdotol3 = 1e-6
    echo -n "$space maximum numerical difference in stdout ($gmaxdif) <= tol ($stdotol3) ? ... "
    if (`echo ' ' | awk -v maxdif=$gmaxdif -v tol=$stdotol3 '{print (maxdif <= tol)}'` == 1) then
      echo yes
    else
      echo no "(max diff = $gmaxdif)"
      unset pass
    endif
  endif

  goto chk3e0
endif  # Al special treatment

# ... Special handling of coptso
if ($ext == "coptso" ) then
   set ehf1 = `cat $testout | grep '^x' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehf2 = `cat $testout | grep '^x' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set ehk1 = `cat $testout | grep '^x' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehk2 = `cat $testout | grep '^x' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set mmom1 = `cat $testout | grep '^x' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set mmom2 = `cat $testout | grep '^x' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`

   set ehfr1 = `zcat $refout | grep '^x' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehfr2 = `zcat $refout | grep '^x' | awk -v varlst="ehf" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set ehkr1 = `zcat $refout | grep '^x' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set ehkr2 = `zcat $refout | grep '^x' | awk -v varlst="ehk" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`
   set mmomr1 = `zcat $refout | grep '^x' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,1p | tail -1`
   set mmomr2 = `zcat $refout | grep '^x' | awk -v varlst="mmom" -v key=. 'BEGIN { nv=split(varlst,vars); } ; { if (! match($0,key)) next ; outstr = " " ; i = 0; while (i++ < nv) {vari = sprintf("%s%s",vars[i],"=[^ \t]*"); k=length(vars[i])+1; if (! match($0,vari)) next ; outstr = sprintf("%s %s",outstr,substr($0,RSTART+k,RLENGTH-k))} print outstr }' | sed -n 1,2p | tail -1`

  if ($?add0) then
    echo -n "         ..." ; $add0 $testout
  else if ($?poszer) then
    echo -n "         ..." ; $poszer $testout
  endif

  if (! $?quiet) then

  call zdiffiles chk3d1 "CPU -1 $testout $refout"
  chk3d1:

  echo "$space Harris energy, original unit cell = $ehf1"
  echo "$space corresponding reference energy    = $ehfr1"
  set ediff = `echo $ehf1 $ehfr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space H-K energy,    original unit cell = $ehk1"
  echo "$space corresponding reference energy    = $ehkr1"
  set ediff = `echo $ehk1 $ehkr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space mag. moment,   original unit cell = $mmom1"
  echo "$space corresponding reference moment    = $mmomr1"
  set ediff = `echo $mmom1 $mmomr1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space Harris energy, supercell          = $ehf2"
  echo "$space corresponding reference energy    = $ehfr2"
  set ediff = `echo $ehf2 $ehfr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space H-K energy,    supercell          = $ehk2"
  echo "$space corresponding reference energy    = $ehkr2"
  set ediff = `echo $ehk2 $ehkr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  echo "$space mag. moment, supercell            = $mmom2"
  echo "$space corresponding reference moment    = $mmomr2"
  set ediff = `echo $mmom2 $mmomr2  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                        = $ediff"
  echo ' '

  if ($?bc) then
    set ehf14 = `echo "4*$ehf1" | bc`
    set ediff = `echo $ehf2 $ehf14  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space Harris energy diff, scell - 4*cell = $ediff"
    set ehk14 = `echo "4*$ehk1" | bc`
    set ediff = `echo $ehk2 $ehk14  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space H-K    energy diff, scell - 4*cell = $ediff"
    set mom14 = `echo "4*$mmom1" | bc`
    set ediff = `echo $mmom2 $mom14 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space Magn. moment diff,  scell - 4*cell = $ediff"
  endif

  endif

if ($?defatol3 == 0) set defatol3 = 2e-6
if ($?dehf1tol3 == 0) set dehf1tol3 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol3 == 0) set dmom1tol3 = 1e-4
if ($?dmomntol3 == 0) set dmomntol3 = 1e-4
if ($?dfmax1tol3 == 0) set dfmax1tol3 = 0.1
if ($?dfmaxntol3 == 0) set dfmaxntol3 = 0.1
if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-4
echo " "

compare_res chk3da "ehf original cell" $ehf1 $ehfr1 $dehf1tol3 pass
chk3da:

compare_res chk3db "ehk original cell" $ehk1 $ehkr1 $dehf1tol3 pass
chk3db:

compare_res chk3dc "ehf supercell    " $ehf2 $ehfr2 $dehf1tol3 pass
chk3dc:

compare_res chk3dd "ehk supercell    " $ehk2 $ehkr2 $dehf1tol3 pass
chk3dd:

compare_res chk3de "mom original cell" $mmom1 $mmomr1 $dmomntol3 pass
chk3de:

compare_res chk3dg "mom supercell    " $mmom2 $mmomr2 $dmomntol3 pass
chk3dg:

goto chk3e0
# end of coptso checks
endif

# ... Special handling testbxc
if ($?testbxc) then
if (! $?vextract) then
  echo "$space missing vextract .. skipping checks"
  goto chk3e0
endif
set ehf1  =  `cat $testout | $vextract c ehf | sed -n 1,1p`
set ehf1r =  `zcat $refout | $vextract c ehf | sed -n 1,1p`

if ($ext == "fe") then
set mmom0  = `cat $testout | $vextract c mmom | sed -n 1,1p`
set mmom0r = `zcat $refout | $vextract c mmom | sed -n 1,1p`
set mmom1  = `cat $testout | $vextract c mmom | sed -n 2,2p`
set mmom1r = `zcat $refout | $vextract c mmom | sed -n 2,2p`
set mmom2  = `cat $testout | $vextract c mmom | sed -n 3,3p`
set mmom2r = `zcat $refout | $vextract c mmom | sed -n 3,3p`

set mmom3  = `cat $testout | $vextract c mmom | sed -n 4,4p`
set mmom3r = `zcat $refout | $vextract c mmom | sed -n 4,4p`
set mmom4  = `cat $testout | $vextract c mmom | sed -n 5,5p`
set mmom4r = `zcat $refout | $vextract c mmom | sed -n 5,5p`
endif

if ($ext == "kfese") then
  set mmom0  = `extract-lines 'START LMF' Exit 2 $testout | grep 'est. true mm' | sed s/=// | awk '{if ($(NF) > 1) print $(NF)}' | tail -1`
  set mmom0r = `extract-lines 'START LMF' Exit 2 $refout | grep 'est. true mm' | sed s/=// | awk '{if ($(NF) > 1) print $(NF)}' | tail -1`
  set mmom1  = `extract-lines 'START LMF' Exit 3 $testout | grep 'est. true mm' | sed s/=// | awk '{if ($(NF) > 1) print $(NF)}' | tail -1`
  set mmom1r = `extract-lines 'START LMF' Exit 3 $refout | grep 'est. true mm' | sed s/=// | awk '{if ($(NF) > 1) print $(NF)}' | tail -1`
endif

if (! $?quiet) then

  call zdiffiles chk3b0 "CPU -1 $testout $refout"
  chk3b0:

  echo "$space first iteration Harris energy    = $ehf1"
  echo "$space first iteration reference energy = $ehf1r"
  set ediff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo
  echo "$space magnetic moment M, unscaled Bxc  = $mmom0"
  echo "$space reference                        = $mmom0r"
  set ediff = `echo $mmom0 $mmom0r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"


  if ($ext == "fe") then
    echo
    echo "$space M, site Bxc scaled by 0.9        = $mmom1"
    echo "$space M, site+istl Bxc scaled by 0.9   = $mmom2"
    set ediff = `echo $mmom2 $mmom1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space change site -> full scaling      =  $ediff"
    set ediff = `echo $mmom0 $mmom1  | awk '{{k=($1/$2)} print k}'`
    echo "$space ratio Mscaled/Munscaled          = $ediff"

    echo
    echo "$space M, site Bxc scaled by 1.1        = $mmom3"
    echo "$space M, site+istl Bxc scaled by 1.1   = $mmom4"
    set ediff = `echo $mmom3 $mmom4  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space change site -> full scaling      =  $ediff"
    set ediff = `echo $mmom0 $mmom3  | awk '{{k=($1/$2)} print k}'`
    echo "$space ratio Mscaled/Munscaled          = $ediff"
    echo
  endif

  if ($ext == "kfese") then
    echo "$space M, site Bxc scaled by 0.9        = $mmom1"
    echo "$space reference                        = $mmom1r"
    set ediff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $ediff"
    set ediff = `echo $mmom1 $mmom0  | awk '{{k=($1/$2)} print k}'`
    echo "$space ratio Mscaled/Munscaled          = $ediff"
    echo
  endif

# End quiet
  endif

# pass checks
if ($?dehf1tol3 == 0) set dehf1tol3 = 2e-6
if ($?dmom1tol3 == 0) set dmom1tol3 = 1e-4

compare_res chk3bb "HF energy, unscaled Bxcf" $ehf1 $ehf1r $dehf1tol3 pass
chk3bb:

if (! $?mmom0) goto chk3bd
compare_res chk3bd "unscaled mag. moment" $mmom0 $mmom0r $dmom1tol3 pass
chk3bd:

if (! $?mmom0) goto chk3be
compare_res chk3be "mag. moment site Bxc scaled by 0.9" $mmom1 $mmom1r $dmom1tol3 pass
chk3be:

if ($ext == "kfese") goto chk3e0

if (! $?mmom0) goto chk3bf
compare_res chk3bf "mag. moment full Bxc scaled by 0.9" $mmom2 $mmom2r $dmom1tol3 pass
chk3bf:

if (! $?mmom0) goto chk3bg
compare_res chk3bg "mag. moment site Bxc scaled by 1.1" $mmom3 $mmom3r $dmom1tol3 pass
chk3bg:

if (! $?mmom0) goto chk3bh
compare_res chk3bh "mag. moment full Bxc scaled by 1.1" $mmom4 $mmom4r $dmom1tol3 pass
chk3bh:

goto chk3e0
# End special handling testbxc
endif

# ... Extract total energies, forces, magnetic moments 1st and last iter
extract_res_n chk32a efa erfa "etot=" 2 0 etot=
chk32a:

set ehf1  =  `cat $testout | grep ehf= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1  =  `cat $testout | grep ehk= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehf1r =  `zcat $refout | grep ehf= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1r =  `zcat $refout | grep ehk= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set dq1   =  `cat $testout | grep 'RMS DQ=' | sed -n 1,1p | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dq1r   = `zcat $refout | grep 'RMS DQ=' | sed -n 1,1p | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqn   =  `cat $testout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`


grep 'Maximum Harris force' $testout >/dev/null
if (! $status) then
  set fmax1  = `cat $testout | grep 'Maximum Harris force' | sed -n 1,1p | awk '{print $5}'`
  set fmax1r = `zcat $refout | grep 'Maximum Harris force' | sed -n 1,1p | awk '{print $5}'`
  set fmaxn  = `cat $testout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
  set fmaxnr = `zcat $refout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
endif

grep mmom= $testout >/dev/null
if (! $status) then
set mmom1  =  `cat $testout      | grep mmom= | sed -n 1,1p | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmom1r =  `zcat $refout | grep mmom= | sed -n 1,1p | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomn  =  `cat $testout      | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomnr =  `zcat $refout | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
if ($?testfsm) then
compare_resf chk32b mmom1u mmom1ur 'Mag. moment:' 3 1 zzz
chk32b:
compare_resf chk32c mmom1 mmom1r 'Mag. moment:' 3 2 zzz
chk32c:
endif
endif



set ediff = `echo $efa $erfa  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} printf "%10.2E", k}'`
if (! $?quiet) then
  echo " "
  if ($efa != "") then
  echo "$space Total energy last free atom      = $efa"
  echo "$space Total energy of reference        = $erfa"
  echo "$space                    difference    =  $ediff"
  echo ' '
  endif

  echo "$space first iteration Harris energy    = $ehf1"
  echo "$space first iteration reference energy = $ehf1r"
  set ediff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space first iteration K-Sham energy    = $eks1"
  echo "$space first iteration reference energy = $eks1r"
  set ediff = `echo $eks1 $eks1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehf1 $eks1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmax1) then
  echo "$space first iteration maximum force    = $fmax1"
  echo "$space first iteration reference force  = $fmax1r"
  endif
  if ($?mmom1) then
    echo "$space first iteration magnetic moment  = $mmom1"
    echo "$space first iteration reference moment = $mmom1r"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif
  if ($?testfsm) then
    echo "$space first iter unconstr. moment      = $mmom1u"
    echo "$space first iter unconstr. ref moment  = $mmom1ur"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmaxn) then
  echo "$space last iteration maximum force     = $fmaxn"
  echo "$space last iteration reference force   = $fmaxnr"
  endif
  if ($?mmom1) then
  echo "$space last iteration magnetic moment   = $mmomn"
  echo "$space last iteration reference moment  = $mmomnr"
  set mdiff = `echo $mmomn $mmomnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"
  endif

  echo "$space last iter RMS input-output drho  = $dqn"
  echo "$space last iter reference RMS drho     = $dqnr"
  echo " "

  zcat $refout | grep RELAX >/dev/null
  if ($status == 0) then
    call showout chk33 RELAX
chk33:
    echo ' '
  endif

  if ($?add0) then
    echo -n "         ..." ; $add0 $testout
  else if ($?poszer) then
    echo -n "         ..." ; $poszer $testout
  endif
  call zdiffiles chk34 "CPU -1 $testout $refout"
chk34:
endif

# ... Check that FA fit basis set is within tol of reference
#  echo "$fitbas" > tmp.$ext.ref
#  if ($?fitbas2) then
#    echo "$fitbas2" >> tmp.$ext.ref
#  endif
#  grep RSMH atm.$ext  > tmp.$ext.lmfa
#  cmp tmp.$ext.lmfa tmp.$ext.ref >/dev/null
#  set retval = $status
#  echo -n "$space lmfa fit basis identical to reference ? ..."
#  if ($retval == 0) then
#   echo yes
#  else
#    echo no
#    unset pass
#  endif

if ($?defatol3 == 0) set defatol3 = 2e-6
if ($?dehf1tol3 == 0) set dehf1tol3 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol3 == 0) set dmom1tol3 = 1e-4
if ($?dmomntol3 == 0) set dmomntol3 = 1e-4
if ($?dfmax1tol3 == 0) set dfmax1tol3 = 0.1
if ($?dfmaxntol3 == 0) set dfmaxntol3 = 0.1
if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-4

# pass checks
chk3c:

# ... Check that FA total energy is within tol of reference
compare_res chk3ca "FA etot (last species)" $efa $erfa $defatol3  pass
chk3ca:

compare_res chk3cb "1st  iter ehf" $ehf1 $ehf1r $dehf1tol3 pass
chk3cb:

if (! $?fmax1) goto chk3cc
compare_res chk3cc "1st  iter max force" $fmax1 $fmax1r $dfmax1tol3 pass
chk3cc:

if ($?testfsm) then
compare_res chk3ccc "1st  iter unconst. mmom" $mmom1u $mmom1ur $dmom1tol3 pass
chk3ccc:
endif

if (! $?mmom1) goto chk3cd
compare_res chk3cd "1st  iter mmom" $mmom1 $mmom1r $dmom1tol3 pass
chk3cd:

compare_res chk3ce "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chk3ce:

if ($?fmaxn) then
compare_res chk3cf "last iter max force" $fmaxn $fmaxnr $dfmaxntol3 pass
chk3cf:
endif

if ($?mmomn) then
compare_res chk3cg "last iter mmom" $mmomn $mmomnr $dmomntol3 pass
chk3cg:
endif

compare_res chk3ch "last iter RMS dq" $dqn $dqnr $drmsqtol3 pass
chk3ch:


# compare bnds to reference
if (-e bnds.$ext) then
  set ndig = 4
  call zcmpnfiles chk3ci "$ndig bnds.$ext $testdir/bnds.$ext"
  chk3ci:
  echo -n "$space ... files bnds.$ext and $testdir/bnds.$ext equivalent to $ndig digits? ... "
  if ($retval == 0) then
    echo  yes
  else
  #    set ndig = 4
  #    call zcmpnfiles chk3cj "$ndig bnds.$ext $testdir/bnds.$ext"
  #    chk3cj:
  #    echo -n "no ... to $ndig digits? ... "
    if ($retval == 0) then
      echo yes
    else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
      echo ok "($retval difference(s) of $ncharfile)"
    else
      echo no "($retval difference(s) remaining of $ncharfile)"
      unset pass
    endif
  endif
endif

set saveehfn = $ehfn

#  xxxx:
#  echo 'FIX'
#
#  set saveehfn = -74.997344
#  set ehfn = -74.424166
#  set refout=$testdir/out.lmf.ionized.$ext testout=out.lmf.$ext
#  zcat $testdir/out.lmf.ionized.$ext >out.lmf.$ext
#  goto chkx32

if ($?testfsm || $?testrsedit || $?testdos) goto chk3e0
echo " "
echo "$space ... repeat for ionized case"
echo " "
set refout=$testdir/out.lmf.ionized.$ext testout=out.lmf.$ext

#goto chkx32

# ... Setup: remove existing files and copy new ones
echo "$space rm -f {mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds}.$ext"
             rm -f {mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds}.$ext
echo "$space cp $cplst ."
             cp $cplst .

# ... Run lmf program
if (! $?clean && $?MPIK) then
  runrdcmd chkx32 %11f $testout "-cat:TMPLION --noerr ctrl.$ext"
else if (! $?clean) then
  runrdcmd chkx32 %11f $testout "-cat:TESTION --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk3e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk3e
endif
chkx32:

# ... Extract total energies, forces, magnetic moments 1st and last iter
extract_res_n chkx32a efa erfa "etot=" 2 0 etot=
chkx32a:
set ehf1  =  `cat $testout | grep ehf= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1  =  `cat $testout | grep ehk= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehf1r =  `zcat $refout | grep ehf= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eks1r =  `zcat $refout | grep ehk= | egrep -v '^   it' | sed -n 1,1p | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
set dq1   =  `cat $testout | grep 'RMS DQ=' | sed -n 1,1p | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dq1r   = `zcat $refout | grep 'RMS DQ=' | sed -n 1,1p | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqn   =  `cat $testout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`

grep 'Maximum Harris force' $testout >/dev/null
if (! $status) then
  set fmax1  = `cat $testout | grep 'Maximum Harris force' | sed -n 1,1p | awk '{print $5}'`
  set fmax1r = `zcat $refout | grep 'Maximum Harris force' | sed -n 1,1p | awk '{print $5}'`
  set fmaxn  = `cat $testout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
  set fmaxnr = `zcat $refout | grep 'Maximum Harris force' | tail -1 | awk '{print $5}'`
endif

grep mmom= $testout >/dev/null
if (! $status) then
set mmom1  =  `cat $testout      | grep mmom= | sed -n 1,1p | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmom1r =  `zcat $refout | grep mmom= | sed -n 1,1p | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomn  =  `cat $testout      | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
set mmomnr =  `zcat $refout | grep mmom= | tail -1 | awk '{sub(".*mmom=","");sub("ehf=.*",""); print $0}'`
endif

set ediff = `echo $efa $erfa  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} printf "%10.2E", k}'`

set ewald = `grep 'Energy for background' out.lmf.$ext | tail -1 | awk '{print $NF}'`

if (! $?quiet) then
  echo " "
  echo "$space Total energy last free atom      = $efa"
  echo "$space Total energy of reference        = $erfa"
  echo "$space                    difference    =  $ediff"
  echo ' '

  echo "$space first iteration Harris energy    = $ehf1"
  echo "$space first iteration reference energy = $ehf1r"
  set ediff = `echo $ehf1 $ehf1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space first iteration K-Sham energy    = $eks1"
  echo "$space first iteration reference energy = $eks1r"
  set ediff = `echo $eks1 $eks1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehf1 $eks1  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"
  if ($?fmax1) then
  echo "$space first iteration maximum force    = $fmax1"
  echo "$space first iteration reference force  = $fmax1r"
  endif
  if ($?mmom1) then
    echo "$space first iteration magnetic moment  = $mmom1"
    echo "$space first iteration reference moment = $mmom1r"
    set mdiff = `echo $mmom1 $mmom1r  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $mdiff"
  endif

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space Harris - Kohn-sham difference    = $ediff"

  if ($?fmaxn) then
  echo "$space last iteration maximum force     = $fmaxn"
  echo "$space last iteration reference force   = $fmaxnr"
  endif
  if ($?mmom1) then
  echo "$space last iteration magnetic moment   = $mmomn"
  echo "$space last iteration reference moment  = $mmomnr"
  set mdiff = `echo $mmomn $mmomnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"
  endif
  echo "$space last iter RMS input-output drho  = $dqn"
  echo "$space last iter reference RMS drho     = $dqnr"
  echo " "

  echo "$space Energy of charged system         = $ehfn"
  echo "$space Estat energy q*q/9/5/<r>         = $ewald"
  set ediff = `echo $ehfn $ewald | awk '{k=$1+$2; print k}'`
  echo "$space Corrected charged system energy  = $ediff"
  echo "$space Energy of neutral system         = $saveehfn"
  set ediff = `echo $ediff $saveehfn | awk '{k=$1-$2; print k}'`
  echo "$space difference                       = $ediff"

  zcat $refout | grep RELAX >/dev/null
  if ($status == 0) then
    call showout chkx33 RELAX
chkx33:
    echo ' '
  endif

  if ($?add0) then
    echo -n "         ..." ; $add0 $testout
  else if ($?poszer) then
    echo -n "         ..." ; $poszer $testout
  endif
  call zdiffiles chkx34 "CPU -1 $testout $refout"
chkx34:
endif

# ... Check that FA fit basis set is within tol of reference
#  echo "$fitbas" > tmp.$ext.ref
#  if ($?fitbas2) then
#    echo "$fitbas2" >> tmp.$ext.ref
#  endif
#  grep RSMH atm.$ext  > tmp.$ext.lmfa
#  cmp tmp.$ext.lmfa tmp.$ext.ref >/dev/null
#  set retval = $status
#  echo -n "$space lmfa fit basis identical to reference ? ..."
#  if ($retval == 0) then
#   echo yes
#  else
#    echo no
#    unset pass
#  endif

if ($?defatol3 == 0) set defatol3 = 2e-6
if ($?dehf1tol3 == 0) set dehf1tol3 = 2e-6
if ($?dehf1toln == 0) set dehf1toln = 2e-6
if ($?dmom1tol3 == 0) set dmom1tol3 = 1e-4
if ($?dmomntol3 == 0) set dmomntol3 = 1e-4
if ($?dfmax1tol3 == 0) set dfmax1tol3 = 0.1
if ($?dfmaxntol3 == 0) set dfmaxntol3 = 0.1
if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-4

# pass checks
chkx3c:

# ... Check that FA total energy is within tol of reference
if ($ext == "copt") then
else
compare_res chkx3ca "FA etot (last species)" $efa $erfa $defatol3  pass
chkx3ca:
endif

compare_res chkx3cb "1st  iter ehf" $ehf1 $ehf1r $dehf1tol3 pass
chkx3cb:

if (! $?fmax1) goto chkx3cc
compare_res chkx3cc "1st  iter max force" $fmax1 $fmax1r $dfmax1tol3 pass
chkx3cc:

if (! $?mmom1) goto chkx3cd
compare_res chkx3cd "1st  iter mmom" $mmom1 $mmom1r $dmom1tol3 pass
chkx3cd:

compare_res chkx3ce "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chkx3ce:

if ($?fmaxn) then
compare_res chkx3cf "last iter max force" $fmaxn $fmaxnr $dfmaxntol3 pass
chkx3cf:
endif

if ($?mmomn) then
compare_res chkx3cg "last iter mmom" $mmomn $mmomnr $dmomntol3 pass
chkx3cg:
endif

compare_res chkx3ch "last iter RMS dq" $dqn $dqnr $drmsqtol3 pass
chkx3ch:

# compare bnds to reference
if (-e bnds.$ext) then
set ndig = 4
call zcmpnfiles chkx3ci "$ndig bnds.$ext $testdir/bnds.$ext"
chkx3ci:
echo -n "$space ... files bnds.$ext and $testdir/bnds.$ext equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
#    set ndig = 4
#    call zcmpnfiles chkx3cj "$ndig bnds.$ext $testdir/bnds.$ext"
#    chkx3cj:
#    echo -n "no ... to $ndig digits? ... "
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) remaining of $ncharfile)"
    unset pass
  endif
endif
endif

chk3e0:
if ($?clean) then
else if ($?pass) then
    echo "$space test 3 PASSED ($ext)"
    set jobspassed = ($jobspassed 3)
else
    echo "$space test 3 FAILED ($ext)"
    set failed = ($failed 3)
endif

chk3e:

echo $joblist | egrep '\b4\b' >/dev/null
if ($status) goto chk4e
set jobid = 4

cat <<EOF

         --- Test 4:  Spin-orbit coupling ---

EOF

if ($?quiet) then
else if ($ext == "felz") then
cat <<EOF
         The felz test computes the orbital moment in Fe.
         lmf calculates the orbital magnetic moment.

           * In the first part of this test only LzSz is used.
             The APW basis with LzSz is also checked.

           * In the second part the FULL SPIN ORBIT is used.

           * In the third part L+S- is added perturbatively to LsSz
             To include the entire L.S perturbatively try
               fp/test/test.fp --so4 felz 4

           * Symmetry operations are reduced because of spin orbit coupling
             (note SOC token in SYMGRP)

           * Demonstrate local k-point shortening (SHORBZ=F)

           * Only 4x4x4 k points are used in this test.

EOF
else if ($ext == "gasls") then
cat <<EOF
         The GaAs test computes the energy bands at (0,0,0) (Gamma point),
         (1/4,1/4,1/4) and (1/2,1/2,1/2) (L point).
         The spin-orbit splitting of the valence states is tested.

         This test checks SO coupling in conjunction with conventional local orbitals.

         A second test includes the L+S- part as a perturbative corrction to LsSz
         To include the entire L.S perturbatively try
            fp/test/test.fp --so4 gasls 4

EOF
else if ($ext == "gaslc") then
cat <<EOF
         The GaAs test with local orbitals and GW self-energy computes
         the energy bands at (0,0,0) (Gamma point),(1/4,1/4,1/4)
         and (1/2,1/2,1/2) (L point).
         The spin-orbit splitting of the valence states is tested.

         This test also checks SO coupling in conjunction with
         extended local orbitals and floating orbitals.

EOF
#  goto chk4c
else if ($ext == "cdte") then
cat <<EOF
         The CdTe case tests fully relativistic local orbitals.

EOF
goto chk4c
endif
if ($?MPIK || $?MPI) then
  if ($ext == "felz") then
  else if ($ext == "cdte") then
    goto chk4c
  else
  echo "$space ... skipping this test ... no MPI-specific checks"
  goto chk4e
  endif
endif
set pass

set refout=$testdir/out.lmf.lzsz.$ext testout=out.lmf.lzsz.$ext
if (! -e $refout) then
  echo "$space ... skipping orbital moments test : missing reference file $refout"
  goto chk4c
endif
echo ' '
query chk41 chk4e 'run this test'
chk41:
set pass
if ($a == "s") goto chk4e
# ... Look for executables
findcmd chk41a rdcmd "$path" "$topdir"
chk41a:
findcmd chk41b lmf "$path" "$topdir"
chk41b:
findcmd chk41c lmfa "$path" "$topdir"
chk41c:

# ... remove related files
echo "$space rm -f {atm,ctrl,fs,moms,mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds,site,dos}.$ext specialspecc specialspeca"
             rm -f {atm,ctrl,fs,moms,mixm,rst,save,log,hssn,wkp,bsmv,syml,bnds,site,dos}.$ext specialspecc specialspeca
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  set testout=out.lmf.ls.$ext
  echo "$space rm -f $testout"
               rm -f $testout
  set testout=out.lmf.$ext
  echo "$space rm -f $testout"
               rm -f $testout
endif
# ... copy required files
echo "$space cp $cplst ."
             cp $cplst .

echo ' '
echo "$space ... Test with LzSz part of SO coupling only"

# ... Run lmf program
if (! $?clean && $?MPIK) then
  runrdcmd chk42 %11f $testout "-cat:MPISZ --noerr ctrl.$ext"
else if (! $?clean) then
  runrdcmd chk42 %11f $testout "-cat:TESTSZ --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk4e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  endif
  goto chk4e
endif
chk42:

  set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set eksn  =  `cat $testout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set eksnr =  `zcat $refout | grep ehk= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehk=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`

egrep ' pwmode=[^0]' $testout >/dev/null
if (! $status) then
  set epw  =  `cat $testout | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set epwr =  `zcat $refout | egrep ' pwmode=[^0]' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
endif

  set orbm  = `cat  $testout| grep "L+ - L-" |  tail -1 | awk '{print $5}'`
  set orbmr = `zcat $refout | grep "L+ - L-" |  tail -1 | awk '{print $5}'`

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  echo " "
  echo "$space last iteration K-Sham energy     = $eksn"
  echo "$space last iteration reference energy  = $eksnr"
  set ediff = `echo $eksn $eksnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"
  set ediff = `echo $ehfn $eksn  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

  if ($?epw) then
    echo " "
    echo "$space last iteration E(MTO + APW)      = $epw"
    echo "$space last iteration ref E(MTO + APW)  = $epwr"
    set ediff = `echo $epw $epwr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
    echo "$space difference                       = $ediff"
  endif

  echo " "
  echo "$space Orbital magnetic moment          = $orbm"
  echo "$space        reference moment          = $orbmr"
  set ediff = `echo $orbm $orbmr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

echo ' '
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif
call zdiffiles chk43 "CPU -1 $testout $refout"
chk43:

# pass checks
chk4achk:
if ($?dehf1toln == 0) set dehf1toln = 2e-6

compare_res chk4achka "Orbital moment" $orbm $orbmr $dorbmtol  pass
chk4achka:
compare_res chk4achkb "last iter ehf" $ehfn $ehfnr $dehf1toln pass
chk4achkb:
compare_res chk4achkc "last iter ehk" $eksn $eksnr $dehf1toln pass
chk4achkc:
if ($?epw) then
compare_res chk4achkd "last iter E(MTO+PW)" $epw $epwr $dehf1toln pass
chk4achkd:
endif

if ($?clean) then
else if ($?pass) then
    echo "$space test 4a PASSED ($ext)"
    set jobspassed = ($jobspassed 4a)
else
    echo "$space test 4a FAILED ($ext)"
    set failed = ($failed 4)
endif

set refout=$testdir/out.lmf.ls.$ext testout=out.lmf.ls.$ext
# ... Run lmf program for full L\dotS
echo ' '
echo "$space ... Calculate orbital moment with FULL SPIN ORBIT"
if ($ext == "felz") then
echo "$space A second band pass is made to test the --ef switch (Fermi level preset)"
echo ' '
endif

if (! $?clean && $?MPIK) then
  runrdcmd chk44a %11f $testout "-cat:MPISO --noerr ctrl.$ext"
else
  runrdcmd chk44a %11f $testout "-cat:TESTSO --noerr ctrl.$ext"
endif
chk44a:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  if ($ext == "felz") then
    set orbm1  = `cat $testout      | grep "L+ - L-" |  sed -n 1,1p | awk '{print $5}'`
    set orbm1r = `zcat $refout | grep "L+ - L-" |  sed -n 1,1p | awk '{print $5}'`
    set orbm2  = `cat $testout      | grep "L+ - L-" |  tail -1 | awk '{print $5}'`
    set orbm2r = `zcat $refout | grep "L+ - L-" |  tail -1 | awk '{print $5}'`
  else
    set orbm1  = `cat $testout      | grep "L+ - L-" |  tail -1 | awk '{print $5}'`
    set orbm1r = `zcat $refout | grep "L+ - L-" |  tail -1 | awk '{print $5}'`
  endif

  echo " "
  echo "$space Orbital magnetic moment WITH FULL SPIN ORBIT= $orbm1"
  echo "$space         reference moment                    = $orbm1r"

  if ($ext == "felz") then
  echo " "
  echo "$space Orbital magnetic moment preset Fermi energy = $orbm2"
  echo "$space         reference moment                    = $orbm2r"
  endif

echo ' '
call zdiffiles chk45 "CPU -1 $testout $refout"
chk45:

# ... Check that orbital moment is within tol of reference
compare_res chk4bchk "Orbital moment" $orbm1 $orbm1r $dorbmtol  pass
chk4bchk:

if ($ext == "felz") then
compare_res chk4cchk "Orbital moment (external Ef)" $orbm2 $orbm2r $dorbmtol  pass
chk4cchk:
endif

if ($?clean) then
else if ($?pass) then
    echo "$space test 4b PASSED ($ext)"
    set jobspassed = ($jobspassed 4b)
else
    echo "$space test 4b FAILED ($ext)"
    set failed = ($failed 4)
endif

#  third SO test : splitting of bands
chk4c:

set refout=$testdir/out.lmf.ls-bands.$ext testout=out.lmf.$ext
echo "$space ... Check band splitting at Gamma with FULL SPIN ORBIT"
if (! -e $refout) then
  echo ' '
  echo "$space ... skipping band splitting test :  missing reference file $refout"
  goto chk4d
endif

echo ' '
query chk4c1 chk4e 'run this test'
chk4c1:
set pass
if ($a == "s") goto chk4e
# ... Look for executables
findcmd chk4c1a rdcmd "$path" "$topdir"
chk4c1a:
findcmd chk4c1b lmf "$path" "$topdir"
chk4c1b:
findcmd chk4c1c lmfa "$path" "$topdir"
chk4c1c:

if ($?haveout) goto chk4c2

# ... remove related files
echo "$space rm -f {ctrl,rst,syml,wkp,bnds,atm,log,mixm}.$ext specialspecc specialspeca"
             rm -f {ctrl,rst,syml,wkp,bnds,atm,log,mixm}.$ext specialspecc specialspeca
# ... copy required files
echo "$space cp $cplst ."
             cp $cplst .

# ... Run lmf program
if (! $?clean && $?MPIK) then
  runrdcmd chk4c2 %11f $testout "-cat:MPISO --noerr ctrl.$ext"
else if (! $?clean) then
  runrdcmd chk4c2 %11f $testout "-cat:TESTSO --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk4e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  goto chk4e
endif
chk4c2:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

#  echo 'uncomment this line'
# set refout=$testdir/out.lmf.ls-bands.$ext testout=out.lmf.$ext

#  set strn    = `extract-lines --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $testout | tail -1`
#  set statesG = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
#  set delta   = `echo $strn | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
# #  set strn    = `extract-lines --gzip --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $refout | tail -1`
#  set strn    = `extract-lines --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $refout | tail -1`
#  set deltar  = `echo $strn | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
#  set strn    = `extract-lines --n=$lineeval  'k=  0.50000  0.50000  0.50000' 1 $testout | tail -1`
#  set statesL = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
#  set strn    = `extract-lines --n=$lineeval  'k=  0.25000  0.25000  0.25000' 1 $testout | tail -1`
#  set statesq = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `

#    echo " "
#    echo "$space The top valence states at:"
#    echo "$space Gamma point        : $statesG Ry"
#    echo "$space (1/4,1/4,1/4) point: $statesq Ry"
#    echo "$space L point            : $statesL Ry"
#    echo "$space Spin-Orbit splitting at the Gamma point = $delta eV"
#    echo "$space        reference splitting              = $deltar eV"

if (! $?quiet) then
echo " "

if ($ext == "cdte") then
   echo "$space" highest evals, SO splitting, gap at Gamma, high-lying p LO
   set gap = `extract-lines 'START LMF' Exit 2 $testout | grep gap | tail -1 | awk '{print $(NF-1)}'`
   set delta = `extract-lines 'START LMF' Exit 3 $testout | grep -A 2 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($7-$5)*13.605}'`
   set states = `extract-lines 'START LMF' Exit 3 $testout | grep -A 2 bndfp: | tail -1 | awk '{print $3, $5, $7}'`
   echo "$space" no LO : $states "  Delta" $delta eV "  gap" $gap eV

   set gap = `extract-lines 'START LMF' Exit 4 $testout | grep gap | tail -1 | awk '{print $(NF-1)}'`
   set delta = `extract-lines 'START LMF' Exit 5 $testout | grep -A 2 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($7-$5)*13.605}'`
   set states = `extract-lines 'START LMF' Exit 5 $testout | grep -A 2 bndfp: | tail -1 | awk '{print $3, $5, $7}'`
   echo "$space" sr LO : $states "  Delta" $delta eV "  gap" $gap eV

   set gap = `extract-lines 'START LMF' Exit 6 $testout | grep gap | tail -1 | awk '{print $(NF-1)}'`
   set delta = `extract-lines 'START LMF' Exit 7 $testout | grep -A 2 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($7-$5)*13.605}'`
   set states = `extract-lines 'START LMF' Exit 7 $testout | grep -A 2 bndfp: | tail -1 | awk '{print $3, $5, $7}'`
   echo "$space" fr LO : $states "  Delta" $delta eV "  gap" $gap eV

   set gapr = `extract-lines 'START LMF' Exit 6 $refout | grep gap | tail -1 | awk '{print $(NF-1)}'`
   set deltar = `extract-lines 'START LMF' Exit 7 $refout | grep -A 2 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($7-$5)*13.605}'`
   set states = `extract-lines 'START LMF' Exit 7 $refout | grep -A 2 bndfp: | tail -1 | awk '{print $3, $5, $7}'`
   echo "$space  " ref   : $states "  Delta" $delta eV "  gap" $gap eV

   echo "$space Corresponding results for low-lying p LO"
   set gap = `extract-lines 'START LMF' Exit 9 $testout | grep gap | tail -1 | awk '{print $(NF-1)}'`
   set delta = `extract-lines 'START LMF' Exit 10 $testout | grep -A 3 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($(NF)-$(NF-1))*13.605}'`
   set states = `extract-lines 'START LMF' Exit 10 $testout | grep -A 3 bndfp: | tail -1 | awk '{print $(NF-3), $(NF-1), $(NF)}'`
   echo "$space" sr LO : $states "  Delta" $delta eV "  gap" $gap eV

   set gap = `extract-lines 'START LMF' Exit 11 $testout | grep gap | tail -1 | awk '{print $(NF-1)}'`
   set delta = `extract-lines 'START LMF' Exit 12 $testout | grep -A 3 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($(NF)-$(NF-1))*13.605}'`
   set states = `extract-lines 'START LMF' Exit 12 $testout | grep -A 3 bndfp: | tail -1 | awk '{print $(NF-3), $(NF-1), $(NF)}'`
   echo "$space" fr LO : $states "  Delta" $delta eV "  gap" $gap eV

else  # Not CdTe
  echo "$space The top valence states at the calculated k-points"
  set i = 0
  set n = `grep -c ' bndfp:  kpt' $testout`
  while ($i < $n)
    @ i = $i + 1
    set strn = `extract-lines --n=$lineeval  'bndfp:  kpt' $i $testout | tail -1`
    set states = `echo $strn | awk -vn1=$eval1 -vn2=$eval2 '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
    extract-lines --n=1 'bndfp:  kpt' $i $testout | awk '{printf "     k = %s %s %s", $7,$8,$9}'
    echo " : $states Ry"
  end
endif  # CdTe or not
endif  # ! quiet

if ($ext == "cdte") then
   set delta = `extract-lines 'START LMF' Exit 7 $testout | grep -A 2 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($7-$5)*13.605}'`
   set strn = `extract-lines 'START LMF' Exit 7 $testout | grep -A 2 bndfp: | tail -1`
   set deltar = `extract-lines 'START LMF' Exit 7 $refout | grep -A 2 bndfp: | tail -1 | awk '{printf "%8.4f\n", ($7-$5)*13.605}'`
   set strnr = `extract-lines 'START LMF' Exit 7 $refout | grep -A 2 bndfp: | tail -1`
else  # Not CdTe

  set strn    = `extract-lines --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $testout | tail -1`
  set delta   = `echo $strn | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
  # set strnr   = `extract-lines --gzip --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $refout | tail -1`
  set strnr   = `extract-lines --n=$lineeval  'k=  0.00000  0.00000  0.00000' 1 $refout | tail -1`
  set deltar  = `echo $strnr | awk -vn=$evalso '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
    echo "$space Spin-Orbit splitting at the Gamma point = $delta eV ( states" `echo $strn | awk -vn=$evalso '{printf "%8.4f and %8.4f \n", $(n),$(n+1)}'` ")"
    echo "$space        reference splitting              = $deltar eV"

endif # CdTe or not

echo ' '
call zdiffiles chk4c3 "CPU -1 $testout $refout"
chk4c3:

echo -n "$space Levels at the Gamma point match reference? ... "
if ("$strn" == "$strnr") then
  echo  yes
else
  echo no
  unset pass
endif
compare_res chk4chkb "SO splitting at the Gamma point" $delta $deltar $gmtol  pass
chk4chkb:

if ($?clean) then
else if ($?pass) then
    echo "$space test 4c PASSED ($ext)"
    set jobspassed = ($jobspassed 4c)
else
    echo "$space test 4c FAILED ($ext)"
    set failed = ($failed 4c)
endif

# fourth SO test: L.S (pert)
chk4d:

set refout=$testdir/out.lmf.lspert-bands.$ext testout=out.lmf.$ext
if ($?so4) set refout=$testdir/out.lmf.lspert4-bands.$ext
echo ' '
echo "$space Check band splitting at Gamma with PERTURBATIVE SPIN ORBIT"
if (! -e $refout) then
  echo ' '
  echo "$space ... skipping L.S(pert) band splitting test :  missing reference file $refout"
  goto chk4de
endif
echo ' '
query chk4d1 chk4e 'run this test'
chk4d1:
set pass
if ($a == "s") goto chk4e
# ... Look for executables
findcmd chk4d1a rdcmd "$path" "$topdir"
chk4d1a:
findcmd chk4d1b lmf "$path" "$topdir"
chk4d1b:
findcmd chk4d1c lmfa "$path" "$topdir"
chk4d1c:

# ... remove related files
echo "$space rm -f {ctrl,rst,syml,wkp,bnds}.$ext"
             rm -f {ctrl,rst,syml,wkp,bnds}.$ext
# ... copy required files
echo "$space cp $cplst ."
             cp $cplst .

# ... Run lmf program
if (! $?clean) then
  set cat = TESTSO3
  if ($?so4) set cat = TESTSO4
  if ($?MPIK) set cat = MPISO3
  runrdcmd chk4d2 %11f $testout "-cat:$cat --noerr ctrl.$ext"
else
  echo "$space rm -f $testout"
               rm -f $testout
  if (-e ctrl.$ext) then
  runrdcmd chk4e %11f . "-cat:CLEAN --noerr ctrl.$ext"
  goto chk4e
endif
chk4d2:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

# Last qp is one that is checked
# if (! $?quiet) then
echo " "
echo "$space The top valence states at the calculated k-points"
set i = 0
set n = `grep -c ' bndfp:  kpt' $testout`
while ($i < $n)
 @ i = $i + 1
 set strn = `extract-lines --quiet --n=$lineevald  'bndfp:  kpt' $i $testout`
 set states = `echo $strn | awk -vn1=$eval1d -vn2=$eval2d '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
 extract-lines --n=1 'bndfp:  kpt' $i $testout | awk '{printf "     k = %s %s %s", $7,$8,$9}'
 echo " : $states Ry"
end
# endif

set delta   = `echo $states | awk -vn=$evalsod '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`
# set strnr = `extract-lines --gzip --quiet --n=$lineevald  'bndfp:  kpt' $i $refout`
set strnr = `extract-lines --quiet --n=$lineevald  'bndfp:  kpt' $i $refout`
set statesr = `echo $strnr | awk -vn1=$eval1d -vn2=$eval2d '{ j=0; while (j++ < NF) if (j >= n1 && j<=n2) printf " %s", $j} END {printf "\n"}' `
set deltar  = `echo $statesr | awk -vn=$evalsod '{printf "%8.4f\n", ($(n+1)-$(n))*13.605}'`

echo "$space Spin-Orbit splitting at the last k-point = $delta eV ( states" `echo $states | awk -vn=$evalsod '{printf "%8.4f and %8.4f \n", $(n),$(n+1)}'` ")"
echo "$space         reference splitting              = $deltar eV"

echo ' '
call zdiffiles chk4d3 "CPU -1 $testout $refout"
chk4d3:

echo -n "$space Levels at the final k-point match reference? ... "
if ("$strn" == "$strnr") then
  echo  yes
else
  echo no
  unset pass
endif
compare_res chk4dhkb "Splitting at the final point" $delta $deltar $gmtol  pass
chk4dhkb:

if ($?clean) then
else if ($?pass) then
    echo "$space test 4d PASSED ($ext)"
    set jobspassed = ($jobspassed 4d)
else
    echo "$space test 4d FAILED ($ext)"
    set failed = ($failed 4)
endif

chk4de:

# fifth SO test : orbital moment with L.S (pert)
set refout=$testdir/out.lmf.lspert.$ext testout=out.lmf.ls.$ext
if ($?so4) set refout=$testdir/out.lmf.lspert4.$ext
if (! -e $refout) then
  echo ' '
  echo "$space ... skipping L.S(pert) orbital moments test :  missing reference file $refout"
  goto chk4e
endif

# ... Run lmf program L\dotS, L+S- perturbation
echo ' '
echo "$space Calculate orbital moment with PERTURBATIVE SPIN ORBIT"
if (! $?clean) then
  set cat = TESTSO3
  if ($?so4) set cat = TESTSO4
  if ($?MPIK) set cat = MPISO3
  runrdcmd chk44e %11f $testout "-cat:$cat --noerr ctrl.$ext"
endif
chk44e:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

  set ehfn  =  `cat $testout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`
  set ehfnr =  `zcat $refout | grep ehf= | egrep -v '^   it' | tail -1 | awk '{match($0,"ehf=[^ ]+"); print substr($0,RSTART+4,RLENGTH-4) }'`

  echo " "
  echo "$space last iteration Harris energy     = $ehfn"
  echo "$space last iteration reference energy  = $ehfnr"
  set ediff = `echo $ehfn $ehfnr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $ediff"

  set orbm1  = `cat $testout      | grep "L+ - L-" |  tail -1 | awk '{print $5}'`
  set orbm1r = `zcat $refout | grep "L+ - L-" |  tail -1 | awk '{print $5}'`

  echo " "
  echo "$space Orbital magnetic moment WITH FULL SPIN ORBIT= $orbm1"
  echo "$space         reference moment                    = $orbm1r"

echo ' '
call zdiffiles chk45e "CPU -1 $testout $refout"
chk45e:

# ... Check that orbital moment is within tol of reference
compare_res chk4echk "Orbital moment" $orbm1 $orbm1r $dorbmtol  pass
chk4echk:

if ($?clean) then
else if ($?pass) then
    echo "$space test 4e PASSED ($ext)"
    set jobspassed = ($jobspassed 4e)
else
    echo "$space test 4e FAILED"
    set failed = ($failed 4)
endif

chk4e:

# --- Summary ---
echo ' '
if ($#failed <= 1) then
    if ($?clean) exit 0
    if ($jobspassed[1] != "") then
    echo "$space $testfile : all tests ($jobspassed) PASSED ($ext)"
    echo " "
    endif
    exit 0
else
    shift failed
    echo "$space $testfile : These tests FAILED:" $failed
    echo " "
    exit -1
endif

# ---------------- runjob --------------
exit
runjob:
  set quitjob=$retcall
  if ($outfile == ".") then
    echo "$space $callarg"
    echo " "
    $callarg
    set retval = $status
    if ($retval != 0) goto cleanup
    goto $quitjob
  endif

  if (`echo $outfile | awk '{print substr($1,1,2)}'` == '>>') then
    set appfile = `echo $outfile | awk '{print substr($1,3)}'`
    echo "$space $callarg  >> $appfile"
    $callarg >> $appfile
    set retval = $status
  else
    echo "$space $callarg  > $outfile"
    $callarg > $outfile
    set retval = $status
  endif
  if ($retval != 0) goto cleanup
  goto $quitjob

# ---------------- compare_resf --------------
# Extracts one element of a line in files $testout and $refout containing a keyword.
# Variables testout and refout point to file names and must be set beforehand ($refout is gzipped file)
# usage: compare_resf retcall testvar refvar keyword arg_number occur_number sed_strn
#   Variables testout and refout referring to file names must be set
#   testvar      : put result from file $testout into this variable
#   refvar       : put result from file $refout (compressed) into this variable
#   keyword      : string line must contain
#   arg_number   : extracts $arg_number'th entry in line, as defined by awk
#   occur_number : argument from $occur_number'th line; if zero, use last line
#   sed_strn     : purge this string from result before assigning
exit
compare_resf:
  set quitjob=$retcall
# echo $retcall $testvar $refvar $keyword $arg_number $occur_number $sed_strn
  set $testvar = `grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  set $refvar = `zcat $refout | grep "$keyword" | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  goto $quitjob

# ---------------- extract_res_n --------------
# Extracts nth token in a line containing a keyword
# usage: extract_res_n retcall testvar refvar keyword arg_number occur_number sed_strn
#   Variables testout and refout referring to file names must be set ($refout is gzipped file)
#   keyword      : string line must contain
#   testvar      : put result from file $testout into this variable
#   refvar       : put result from file $refout (compressed) into this variable
#   arg_number   : extracts $arg_number'th entry in line, as defined by awk
#   occur_number : argument from $occur_number'th line; if zero, use last line
#   sed_strn     : delete this string with from result before assigning
exit
extract_res_n:
  set quitjob=$retcall
#   echo $retcall $testvar $refvar keyword=$keyword $arg_number $occur_number $sed_strn
  set $testvar = `grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  set $refvar = `zcat $refout | grep "$keyword" | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  goto $quitjob

# ---------------- runrdcmd --------------
exit
runrdcmd:
  set quitjob=$retcall
  if ($outfile == ".") then
    $rdcmd -f:$rdcmdfmt $callarg
    set retval = $status
    echo ' '
    if ($retval == 0) then
      echo "$space Job(s) completed successfully"
      goto $quitjob
    endif
  else
    if (`echo $outfile | awk '{print substr($1,1,2)}'` == '>>') then
      set appfile = `echo $outfile | awk '{print substr($1,3)}'`
      echo "$space $callarg  >> $appfile"
      exit
#      $callarg >> $appfile
      set retval = $status
    else
      echo "$space ... the following job(s) will be executed by invoking "\""rdcmd $callarg"\" "(> $outfile)"
      $rdcmd -f:$rdcmdfmt --n $callarg
      echo "$space ... starting invocation of rdcmd:"
      echo "$space $rdcmd '-f:#rdcmd:%2f' $callarg  >& $outfile"
      $rdcmd '-f:rdcmd:%2f' $callarg >& $outfile
      set retval = $status
    endif
  endif

  if ($retval == 0) then
    echo "$space Job(s) completed successfully; output in $outfile"
    if ($?poszer) then
      echo -n "         ..." ; $poszer $outfile
    else if ($?add0) then
      echo -n "         ..." ; $add0 $outfile
    endif
    goto $quitjob
  else
    echo "$space ...oops... the following command returned with nonzero exit status:"
    echo -n "$space   "
    grep $rdcmd:t{:} $outfile | tail -1 | sed 's/rdcmd:  //'
    goto cleanup
  endif

# ---------------- cleanup --------------
exit
cleanup:
  if ($retval != 0) echo "$space job returned with error status $retval"
  if ($retval != 0) echo "$space ... $testfile aborting, test $ext job $jobid FAILED ($ext)"
  exit $retval

# ---------------- diffiles --------------
# calling argument should consist of four strings:
# 1st string = string that terminates diff
# 2nd string = integer that counts how many times terminator should occur before terminating
# 3nd string = first file name
# 4th string = second file name
# example: call diffiles chk69 "CPU 3 $testout $refout"
exit
diffiles:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  set files = ($callarg)
  set endstr = $files[1]
  shift files
  set nend = $files[1]
  shift files
  if ($nend == "-1") then
    set nend = `grep "$endstr" $files[1] | wc | awk '{print $1}'`
  endif

#    echo difffiles : $quitjob $nend
#    grep $endstr $files[1]

  query diff11 $quitjob "compare $files"
diff11:
  diff $files | awk -v endstr=$endstr -v nend=$nend -v endl=0 -v endr=0 '{if ($1 == "<" && endl < nend) print ; if ($1 == ">" && endr < nend) print ; if ($1 == ">" || $1 == "<" || endl >= nend && endr >= nend) ; else {print} ; if ($1 == "<" && $2 == endstr) {endl+=1}; if ($1 == ">" && $2 == endstr) {endr+=1};}' | sed -n 1,50p
  goto $quitjob

# ---------------- zdiffiles --------------
# File differences, with additional check for numerical differences
# callarg should consist of four strings; there is an optional fifth and sixth
# 1st word = string that terminates diff
# 2nd word = counts how many times terminator should occur before terminating
#            -1 -> last occurence
# 3nd word = first file name
# 4th word = second file name
# 5th word = (optional) tolerance.  Numerical differences < tolerance are counted as 0
#            If present, and not "-", passed to mcx as the argument to ~tol=
# 6th word = (optional) if present, it is used instead of mcexcl
#
# Returns ndif = number of differences, and maxdif = difference (if mcx is available)
# Example: call zdiffiles chk69 "CPU 3 $testout $refout"
exit
zdiffiles:

  set quitjob=$retcall

  set noglob
  set files = ($callarg)
  unset noglob

  set endstr = $files[1]
  shift files
  set nend = $files[1]
  shift files

  if ($nend == "-1") then
    set nend = `grep "$endstr" $files[1] | wc | awk '{print $1}'`
  endif

  if ($?quiet) goto zdiffiles2

  if ( $?slow == 0 ) echo "$space ... compare $files[1] $files[2]"
  query zdiff1 $quitjob "compare $files[1] $files[2]"
zdiff1:
  zdiff -Icpudel -Iwritten $files[1] $files[2] | awk -v endstr="$endstr" -v nend=$nend -v endl=0 -v endr=0 '{if ($1 == "<" && endl < nend) print ; if ($1 == ">" && endr < nend) print ; if ($1 == ">" || $1 == "<" || endl >= nend && endr >= nend) ; else {print} ; if ($1 == "<" && $2 == endstr) {endl+=1}; if ($1 == ">" && $2 == endstr) {endr+=1};}' | sed -n 1,50p
  echo " "

zdiffiles2:
  if (! $?mcx) goto $quitjob

  if (! $?mcexcll) set mcexcll
  if (! $?mcexcl) set mcexcl
  if ($?mcexcl) set mcexcll = "$mcexcl"
  if ($#files > 3) then
    set mcexcll = "$files[4]"
  endif

  if (! $?mcterm) set mcterm
# Don't do this step ... check that each test initializes its own gmaxdif
# if (! $?gmaxdif) set gmaxdif = 0
  set toldif
  if ($#files > 2) then
    if ("$files[3]" != "-") set toldif = "~tol=$files[3]"
  endif

  set maxdif = `$mcx -cmpf~fn1=$files[1]~fn2=$files[2]~max$toldif$mcterm$mcexcll`
  set ndif = `$mcx -cmpf~fn1=$files[1]~fn2=$files[2]~ndiff$toldif$mcterm$mcexcll`
#  set gmaxdif = `echo $gmaxdif $maxdif  | awk '{print ($1>$2)?$1:$2}'`

  echo "$space $ndif numerical differences in $files[1] compared to ref, max diff = $maxdif"
  if ($?slow > 0 && $?verb) then
    echo
    echo "$space The following make a detailed comparison of numerical differences:"
    echo "$space $mcx -cmpf~fn1=$files[1]~fn2=$files[2]~vverb$toldif$mcterm$mcexcll"

    query zdiff2  $quitjob "show comparison"
zdiff2:
    $mcx -cmpf~fn1=$files[1]~fn2=$files[2]~vverb$toldif$mcterm$mcexcll
    echo
    echo "$space"'*'"hit <return> to continue"
    set a = ($<)
  endif
  if (! $?quiet) echo
  goto $quitjob

# ---------------- zdiffilesx --------------
# Identical to zdiffiles, but gmaxdif is accumulated
# callarg should consist of four strings; there is an optional fifth and sixth
# 1st word = string that terminates diff
# 2nd word = counts how many times terminator should occur before terminating
#            -1 -> last occurence
# 3nd word = first file name
# 4th word = second file name
# 5th word = (optional) tolerance.  Numerical differences < tolerance are counted as 0
#            If present, and not "-", passed to mcx as the argument to ~tol=
# 6th word = (optional) if present, it is used instead of mcexcl
#
# Returns ndif = number of differences, and maxdif = difference (if mcx is available)
exit
zdiffilesx:

  set quitjob=$retcall

  set noglob
  set files = ($callarg)
  unset noglob

  set endstr = $files[1]
  shift files
  set nend = $files[1]
  shift files

  if ($nend == "-1") then
    set nend = `grep "$endstr" $files[1] | wc | awk '{print $1}'`
  endif

  if ($?quiet) goto zdiffilesx2

  if ( $?slow == 0 ) echo "$space ... compare $files[1] $files[2]"
  query zdiffx1 $quitjob "compare $files[1] $files[2]"
zdiffx1:
  zdiff -Icpudel -Iwritten $files[1] $files[2] | awk -v endstr="$endstr" -v nend=$nend -v endl=0 -v endr=0 '{if ($1 == "<" && endl < nend) print ; if ($1 == ">" && endr < nend) print ; if ($1 == ">" || $1 == "<" || endl >= nend && endr >= nend) ; else {print} ; if ($1 == "<" && $2 == endstr) {endl+=1}; if ($1 == ">" && $2 == endstr) {endr+=1};}' | sed -n 1,50p
  echo " "

zdiffilesx2:
  if (! $?mcx) goto $quitjob

  if (! $?mcexcll) set mcexcll
  if (! $?mcexcl) set mcexcl
  if ($?mcexcl) set mcexcll = "$mcexcl"
  if ($#files > 3) then
    set mcexcll = "$files[4]"
  endif

  if (! $?mcterm) set mcterm
# Don't do this step ... check that each test initializes its own gmaxdif
# if (! $?gmaxdif) set gmaxdif = 0
  set toldif
  if ($#files > 2) then
    if ("$files[3]" != "-") set toldif = "~tol=$files[3]"
  endif

  set maxdif = `$mcx -cmpf~fn1=$files[1]~fn2=$files[2]~max$toldif$mcterm$mcexcll`
  set ndif = `$mcx -cmpf~fn1=$files[1]~fn2=$files[2]~ndiff$toldif$mcterm$mcexcll`
  if (! $?gmaxdif) set gmaxdif = 0
  set gmaxdif = `echo $gmaxdif $maxdif  | awk '{print ($1>$2)?$1:$2}'`

  echo "$space $ndif numerical differences in $files[1] compared to ref, max diff = $maxdif"
  if ($?slow > 0 && $?verb) then
    echo
    echo "$space The following make a detailed comparison of numerical differences:"
    echo "$space $mcx -cmpf~fn1=$files[1]~fn2=$files[2]~vverb$toldif$mcterm$mcexcll"

    query zdiffx2  $quitjob "show comparison"
zdiffx2:
    $mcx -cmpf~fn1=$files[1]~fn2=$files[2]~vverb$toldif$mcterm$mcexcll
    echo
    echo "$space"'*'"hit <return> to continue"
    set a = ($<)
  endif
  if (! $?quiet) echo
  goto $quitjob

# ---------------- compare_res --------------
# Compares two numbers $testvar-$refvar and unsets $passvar if |testvar-refvar|<tol
# usage: compares_res retcall keyword testvar refvar tol passvar
#   keyword      : label (for printout)
#   testvar      : first number
#   refvar       : second number
#   tol          : tolerance
#   passvar      : $passvar is unset if |testvar-refvar|<tol
exit
compare_res:
  set quitjob=$retcall
# echo $retcall $keyword $testvar $refvar $tol $passvar
  set toll = `echo $tol | sed s/d/e/`
  echo -n "$space $keyword ($testvar) within tol ($toll) of reference ($refvar)? ... "
  if (`echo $testvar $refvar | awk -v tol=$toll '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjob

# ---------------- compare_res_0 --------------
# Compares a number $testvar and unsets $passvar if |testvar|<tol
# usage: compares_res_0 retcall keyword testvar tol passvar
# Example:
# compare_res_0 chk274a "Max deviation in pdos from reference" $retval $pdostol pass
#   keyword      : label (for printout)
#   testvar      : first number
#   tol          : tolerance
#   passvar      : $passvar is unset if |testvar|<tol
exit
compare_res_0:
  set quitjob=$retcall
#  echo $retcall $keyword $testvar $tol $passvar
  set toll = `echo $tol | sed s/d/e/`
  echo -n "$space $keyword ($testvar) within tol ($toll)? ... "
  if (`echo $testvar 0 | awk -v tol=$toll '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjob

# ---------------- zcmpmfiles_res_0 --------------
# Compares two files, stripping all but numerical fields.
# Checks for max absolute difference and unsets $passvar if difference<$tol
# Files with .gz or .Z extensions are assumed to be gzipped.
# usage: zcmpmfiles_res_0 retcall keyword testvar tol passvar ndig srcfile reffile
# See also zcmpmfiles_res_tol, which accomplished the same thing but extra argument nlines
#   retcall      : return to this point in script on exit
#   keyword      : label (for printout)
#   tol          : tolerance in maximum allowed deviation
#   passvar      : $passvar is unset if |testvar|<tol
#   ndig         : number of digits numbers in file are stripped to
#   srcfile      : first file to compare
#   reffile      : second file to compare
#   nlines       : number of lines to compare. Use 0 for all lines.  Inoperative if either file is a zipped file.
# Example:
# zcmpmfiles_res_0 chk1ck "Max deviation in bnds.$ext from reference" $bndstol pass 4 bnds.$ext $testdir/bnds.$ext.gz
exit
zcmpmfiles_res_0:
  set quitjobl=$retcall
# echo $retcall $keyword $tol $?passvar $ndig $srcfile $reffile

  unset retval
  call zcmpmfiles zcmpmfilesx "$ndig $srcfile $reffile"
zcmpmfilesx:
  echo -n "$space $keyword ($retval) within tol ($tol)? ... "
  if (`echo $retval 0 | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjobl

# ---------------- zcmpmfiles_res_tol --------------
# Compares two files, stripping all but numerical fields.
# Checks for max absolute difference and unsets $passvar if difference<$tol
# Files with .gz or .Z extensions are assumed to be gzipped.
# usage: zcmpnfiles_res_tol retcall keyword testvar tol passvar ndig srcfile reffile nlines
# See also zcmpmfiles_res_0, which accomplished the same thing but without nlines
#   retcall      : return to this point in script on exit
#   keyword      : label (for printout)
#   tol          : tolerance in maximum allowed deviation
#   passvar      : $passvar is unset if |testvar|<tol
#   ndig         : number of digits numbers in file are stripped to
#   srcfile      : first file to compare
#   reffile      : second file to compare
#   nlines       : number of lines to compare. Use 0 for all lines.  Inoperative if either file is a zipped file.
# Example:
# zcmpmfiles_res_tol chk1ck "Max deviation in bnds.$ext from reference" $bndstol pass 4 bnds.$ext $testdir/bnds.$ext.gz 0
exit
zcmpmfiles_res_tol:
  set quitjobl=$retcall
# echo $retcall $keyword $tol $?passvar $ndig $srcfile $reffile $nlines

  unset retval
  if ($nlines == 0) then
    call zcmpmfiles zcmpmfilesy "$ndig $srcfile $reffile"
  else
    call zcmpmfiles zcmpmfilesy "nlines=$nlines $ndig $srcfile $reffile"
  endif
zcmpmfilesy:

  echo -n "$space $keyword ($retval) within tol ($tol)? ... "
  if (`echo $retval 0 | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjobl

# ---------------- zcmpnfiles --------------
# Compares two files, treating each field as a number.
# call arguments should contain 3 strings: n test-file reference-file
# |n| = number of digits which numbers are truncated to.
# If n<0, sort files before comparing them
# Alternatively call arguments can contain 4 strings : nlines=# n test-file reference-file
# nlines=# specifies that the check is made on the first # lines only
# Files with .gz or .Z extensions are assumed to be gzipped.
# (nlines doesn't work with gzipped files; sorry)
# Returns with retval = number of differences in reduced files.
# Example :  call zcmpnfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $tmpdir/tmp_compnfile_1 $tmpdir/tmp_compnfile_2
exit
zcmpnfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)

  set nlines ; unset nlines
  switch ($zcmpnargs[1])
    case "nlines=*":
      set nlines = `echo $zcmpnargs[1] | sed s/nlines=//`
      @ nlines = $nlines  # Checks to make sure this is an integer
      shift zcmpnargs

    default:
  endsw

  set lsort; unset lsort
  @ digits = $zcmpnargs[1]
  if ($digits < 0) then
    @ digits = - $digits
    set lsort
  endif
# set a = ' { for (i = NF; i > 0; --i) printf " %.'$digits'f", $i; printf "\n" }'
# set a = ' { for (i = 1; i <= NF; i++) { k = sprintf("%.'$digits'f",$i); if (k+k == 0) k = 0 ; printf "%s ", k}; printf "\n" }'
  set a = ' { for (i = 1; i <= NF; i++) { k = sprintf("%.'$digits'f",$i); if (k+k == 0) k = 0 ; printf "%s ", k}; print "" }'

  set fn1 = $tmpdir/tmp_compnfile_1
  set fn2 = $tmpdir/tmp_compnfile_2
  if ("$zcmpnargs[2]:e" == 'gz' || "$zcmpnargs[2]:e" == 'Z') then
    set cat1 = 'gunzip -c'
  else
    set cat1 = cat
    if ($?nlines) then
      set cat1 = "head -$nlines"
    endif
  endif
  if ("$zcmpnargs[3]:e" == 'gz' || "$zcmpnargs[3]:e" == 'Z') then
    set cat2 = 'gunzip -c'
  else
    set cat2 = cat
    if ($?nlines) then
      set cat2 = "head -$nlines"
    endif
  endif

  # if (! $?quiet) then
  #   if ($?lsort) then
  #   if ($?exclude) then
  #     echo "$cat1 $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk '"$a"' >" $fn1
  #     echo "$cat2 $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk '"$a"' >" $fn2
  #   else
  #     echo "$cat1 $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk '"$a"' >" $fn1
  #     echo "$cat2 $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk '"$a"' >" $fn2
  #   endif
  #   else
  #   if ($?exclude) then
  #     echo "$cat1 $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | awk '"$a"' >" $fn1
  #     echo "$cat2 $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | awk '"$a"' >" $fn2
  #   else
  #     echo "$cat1 $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | awk '"$a"' >" $fn1
  #     echo "$cat2 $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | awk '"$a"' >" $fn2
  #   endif
  #   endif
  # endif

  if ($?lsort) then
    if ($?exclude) then
      $cat1  $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk "$a" > $fn1
      $cat2  $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk "$a" > $fn2
  else
      $cat1  $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk "$a" > $fn1
      $cat2  $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk "$a" > $fn2
  endif
  else

    if ($?exclude) then
      $cat1  $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
      $cat2  $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2
    else
      $cat1  $zcmpnargs[2] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
      $cat2  $zcmpnargs[3] | sed 's:\([1-9]\)-:\1 -:g' | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2
    endif
  endif

  set ncharfile = `wc $fn1 | awk '{print $3}'`
  set nwordfile = `wc $fn1 | awk '{print $2}'`
  set nlinefile = `wc $fn1 | awk '{print $1}'`
  cmp $fn1 $fn2 >/dev/null
  set retval = $status

  if ($retval == 0) rm -f $fn1 $fn2
  if ($retval == 0) goto $quitjob

  set retval = `cmp -l $fn1 $fn2 |& grep -v EOF | wc | awk '{printf "%d", $1}'`
  if ($retval == 0) set retval = '-1'
  rm -f $fn1 $fn2
  goto $quitjob

# ---------------- zcmpmfiles --------------
# Compares two files, treating each field as a number.
# Call arguments should contain 3 strings: no-digits test-file reference-file
# Alternatively call arguments can contain 4 strings : nlines=# n test-file reference-file
# nlines=# specifies that the check is made on the first # lines only
# Files with .gz or .Z extensions are assumed to be gzipped.
# Returns with retval = max numerical difference
# Example :  call zcmpmfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $tmpdir/tmp_compnfile_1 $tmpdir/tmp_compnfile_2
exit
zcmpmfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)

  set nlines ; unset nlines
  switch ($zcmpnargs[1])
    case "nlines=*":
      set nlines = `echo $zcmpnargs[1] | sed s/nlines=//`
      @ nlines = $nlines  # Checks to make sure this is an integer
      shift zcmpnargs

    default:
  endsw

  set digits = $zcmpnargs[1]
# set a = ' { for (i = NF; i > 0; --i) printf " %.'$digits'f", $i; printf "\n" }'
  set a = ' { for (i = 1; i <= NF; i++) { k = sprintf("%.'$digits'f",$i); if (k+k == 0) k = 0 ; printf "%s ", k}; printf "\n" }'

  set fn1 = $tmpdir/tmp_compnfile_1
  set fn2 = $tmpdir/tmp_compnfile_2
  if ("$zcmpnargs[2]:e" == 'gz' || "$zcmpnargs[2]:e" == 'Z') then
    set cat1 = 'gunzip -c'
  else
    set cat1 = cat
    if ($?nlines) then
      set cat1 = "head -$nlines"
    endif
  endif
  if ("$zcmpnargs[3]:e" == 'gz' || "$zcmpnargs[3]:e" == 'Z') then
    set cat2 = 'gunzip -c'
  else
    set cat2 = cat
    if ($?nlines) then
      set cat2 = "head -$nlines"
    endif
  endif

  $cat1  $zcmpnargs[2] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
  $cat2  $zcmpnargs[3] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2

  set retval = `diff -y --width=300 $fn1 $fn2 | grep '|' | awk -v top=0 '{n=split($0,a,"|"); n1=split(a[1],b1); n2=split(a[2],b2); { j=0; while (j++ < n1) if (j <= n1 && j<=n2) {x = (b1[j]-b2[j])>0?(b1[j]-b2[j]):(b2[j]-b1[j]); top = (top-x)>0?top:x; }}} END {printf "%12.4e\n", top}'`
  rm -f $fn1 $fn2
  goto $quitjob

# ---------------- zcmpmcx --------------
# Compares numerical elements in two files read by mcx -cmpf
# Checks for max absolute difference and unsets $passvar if
# (1) difference>$toldif
# (2) fewer ndiff differences encountered
# usage: zcmpmcx retcall keyword testvar toldif ndiff prec passvar srcfile reffile
#   retcall      : return to this point in script on exit
#   keyword      : label (for printout).
#   toldif       : tolerance in maximum allowed deviation.  toldif<0 not used
#   ndiff        : number of differences allowed.  ndiff<0 not used
#   prec         : precision with with numbers are compared
#   passvar      : if $passvar is null, check prints out keyword
#                : otherwise, check prints out keyword only, if it is not blank
#                : and varable $passvar is unset if criteria (1) and (2) are not met
#   srcfile      : first file to compare
#   reffile      : second file to compare
# variable $mcx must point to the mcx calculator
# Returns maxdif = max abs diff  and ndif = number of differences > prec
# Returns retval = 0 if both criteria satisfied, retval > 0 if one is not satisfied
#
# Case passvar is "null", what is printed out:
# If keyword is empty zcmpmcx is silent.  Example:
#   zcmpmcx chk266 "" -1e-1 20 1e-4 null dos-mull.$ext $testdir/dos-mull.$ext
# If keyword is not empty and passvar is null, zcmpmcx prints out keyword only.
#   zcmpmcx chk266 "" -1e-1 20 1e-4 null dos-mull.$ext $testdir/dos-mull.$ext
#
# Case passvar is not null, $passvar may be unset.  zcmpmcx prints out an internal string
# If keyword is empty it zcmpmcx makes its own header, "check $srcfile".  Example
#   zcmpmcx chk266a "" -1 40 1e-4 pass dos-mull.$ext $testdir/dos-mull.$ext
# Prints out something like:
#   check dos-mull.gdn : ndiff (33) within tol (40)? ... yes
# Otherwise zcmpmcx uses keyword for the header.  Example:
#   zcmpmcx chk266a "... checking" -1 40 1e-4 pass dos-mull.$ext $testdir/dos-mull.$ext
#   prints out something like
#   ... checking ndiff (33) within tol (40)? ... yes
exit
zcmpmcx:
  set quitjobl=$retcall
# echo $retcall $keyword $toldif $ndiff $prec $passvar $srcfile $reffile

  if (! $?mcx) then
    echo "no mcx in path ... no check made"
    goto $retcall
  endif

  set retval = 0

  set mcxarg = "-cmpf~fn1=$srcfile~fn2=$reffile"
  if (! $?mcincl) set mcincl
  if (! $?mcexcl) set mcexcl
  if (! $?mcterm) set mcterm
  set mcxarg = "$mcxarg$mcexcl$mcincl$mcterm"
  if ($?prec) then
    set mcxarg = "$mcxarg~tol=$prec"
  endif
  set maxdif = `$mcx $mcxarg~max`
  set ndif   = `$mcx $mcxarg~ndiff`

  if ("$keyword" != "") then
    echo -n "$space $keyword"
  endif

  set res = "yes"
  if ($ndiff > 0 && `echo $toldif | awk '{print  ($1 > 0)}'`) then
     if (`echo $maxdif $toldif | awk '{print  ($1 <= $2)}'` == 0) set res = "no"
     if ($res == "no") set retval = 1
     if ($passvar != "null") then
       if ("$keyword" == "") then
         echo -n "$space check $srcfile : "
       endif
       echo -n " max deviation ($maxdif) within tol ($toldif)? ..." $res";"
     endif
     set res = "yes"
     if ($ndif >= $ndiff) set res = "no"
     if ($passvar != "null") then
       echo "  ndiff ($ndif) within tol ($ndiff)? ..." $res
     endif
     if ($res == "no") set retval = 1
  else if (`echo $toldif | awk '{print  ($1 > 0)}'`) then
     if (`echo $maxdif $toldif | awk '{print  ($1 <= $2)}'` == 0) set res = "no"
     if ($passvar != "null") then
       if ("$keyword" == "") then
         echo -n "$space check $srcfile : "
       endif
       echo " max deviation ($maxdif) within tol ($toldif)? ..." $res
     endif
     if ($res == "no") set retval = 1
  else if ($ndiff > 0) then
     if ($ndif >= $ndiff) set res = "no"
     if ($res == "no") set retval = 1
     if ($passvar != "null") then
       if ("$keyword" == "") then
         echo -n "$space check $srcfile : "
       endif
       echo " ndiff ($ndif) within tol ($ndiff)? ..." $res
     endif
#  else  No checks made

  endif

  if ($passvar != "null") then
    if ($retval > 0) unset $passvar
  endif
  goto $retcall

# ---------------- qprint (print only quiet not set) --------------
exit
qprint:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  echo "$callarg"
  goto $quitjob

# ---------------- showout --------------
exit
showout:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  echo ' '
  echo "$space ... Compare $callarg to line(s) in file $refout":
  grep "$callarg" $testout
  if (`cat $testout | grep "$callarg" | wc | awk '{print $1}'` > 1) echo ' ---'
  zcat $refout | grep "$callarg"
  goto $quitjob

# ---------------- findcmd --------------
# Finds an executable program within the supplied path
# Usage: findcmd return_label executable_command path_name make_path
# If $executable_command is not found, findcmd does one of the following:
# If make_path = 'no' : returns silently.
# Otherwise findcmd aborts with a message, which assumes
# $make_path is the path where $executable_command is made.
exit
findcmd:
set found = 'no'
foreach ac_dir ($path_name)
 if (-x $ac_dir/$prog_cmd) then
   set $prog_cmd = $ac_dir/$prog_cmd
   set found = 'yes'
   break
 endif
end
if (! $?quiet) then
  if ($found == 'yes') echo "$space ... using executable $ac_dir/$prog_cmd"
  if ($found == 'no')  echo "$space ... no executable $prog_cmd found in path"
endif

# echo make_path $make_path

if ($found == 'no' && $make_path != "optional") then
  echo "  "
  echo "  Sorry, $testfile cannot find program '"$prog_cmd"' it needs to execute."
  echo "  '"$prog_cmd"' was not found in supplied path, or in the following:"
  echo "        $topdir"
# echo "  ... This script ($testfile) requires binary "'"rdcmd"'" to run."
  echo "  You must create or put '"$prog_cmd"' in your path before invoking this script."
  echo "  Normally '"$prog_cmd"' is created as part of the installation process."
  echo "  Invoking '"make $prog_cmd"' in $make_path should create it."
  echo "  $testfile aborting ..."
  exit -1
endif
goto $retcall

# ---------------- query --------------
exit
query:
  unset skip
  if ($?slow != 0) then
    echo "$space"'*'"hit <return> to $callarg, s <return> to skip it."
    set a = ($<)
    if ($a == "") goto $retcall
    switch ($a)
      case "quit":
      case "q":
      case "a":
        exit
      case "i":
        unset slow
        breaksw
      case "s":
        set skip
        breaksw
      case "t":
        time
        goto query
      default:
        echo 'q to quit; i unsets slow; s skips this job, t shows time'
        goto query
    endsw
  endif
  if ($?skip) goto $retcall2
  goto $retcall

# ---------------- List tests --------------
showtests:
cat <<EOF
  Usage: invoke with:   $testfile [switches] material-name job-list

   Material:    tests
       copt:    sampling with Fermi distribution, forces, special two-kappa basis, perturbation treatment of core
                Also shows --cvK switch to print LDA estimate of specific heat.
                Job 3 demonstrates restart file editor (combining densities)
         te:    molecular statics in an open structure, floating orbitals, use of GGA
        zrt:    ZrO_2 fluorite in tetragonal setting, with tetragonal distortion.  Tests read of ASCII rst file
         co:    tetrahedron/sampling, branches of charge mixing, color weights for bands, color weights with SO coupling
                Job 1 demonstrates various features of the band drawing mode (--band~col=...)
                Job 2 demonstrates partial dos (--pdos:mode=2:sites=2)
                Job 3 demonstrates restart file editor (making a supercell)
     cr3si6:    insulator, forces with correction mode 1, pnu floated with energy-weighted density matrix
                Job 2 tests m-resolved partial DOS (--pdos~mode=2~group=1:3;4:9~lcut=2,1)
         fe:    spin-polarized calculation, joint DOS, optics editor (--popted)
                Job 1 Also demonstrates how to generate bands uniform mesh in a plane, to draw Fermi surfaces
                Job 2 tests Mulliken dos (--mull:mode=2) and core level spectroscopy (--cls)
         cu:    A test of high-lying local orbitals and compares a pure LMTO with a mixed LMTO, LAPW basis
                Job 1 generates energy bands in a mixed LMTO, LAPW basis
         au:    A test designed generate equilibrium volume to high accuracy. Compared against standard as described in
                http://molmod.ugent.be/sites/default/files/deltadftcodes/supplmat/SupplMat-WIEN2k.pdf
        gas:    Tests local 3d orbital, and joint DOS
                Job 1 generates energy bands with empty and space-filling spheres,
                and generates bands on a union coarse, uniform mesh over the BZ and a dense mesh near Gamma
                Job 2 tests partial DOS (--pdos~lcut=2,2,1,1~mode=1)
                Job 3 tests band edge finder
     srtio3:    Frozen augmented wave functions, APW addition to basis, extended local orbitals
                Job 1 tests renormalized free O atom (improved starting density), and re-use of a restart file with shift basis vectors
                Job 2 tests 90 degree rotation of lattice DOS (ROT=y:pi/2)
                Job 3 checks command-line switches --quit=dos --dos
       cdte:    Job 1 Tests LDA+U implementation in various modes, and shear (SHEAR=)
                Job 4 tests SO coupling, including Dirac local orbitals
       felz:    Different approximations to SO coupling, with APW basis
                Job 1 also tests fixed-spin moment
                Job 4 tests magnetic symmetry (SOC)
      gasls:    SO coupling in conjunction with conventional local orbitals.
      zbgan:    joint DOS and Im eps for GaN in the zincblende structure, including nonequilibrium cases
           :    Job 3 writes the pseudodensity for the highest lying valence band in Zincblende GaN.
        gdn:    LDA+U, and also LDA with partial core occupancy, application of external B-field
                Job 2 tests partial DOS (--mull:mode=1)
       eras:    LDA+U
          c:    homogeneous background
       fept:    Extended local orbitals in a spin polarized environment; supercells of CuAu structure
         al:    Demonstrate k-resolved partial DOS at Fermi energy
        mgo:    External potential, with augmentation spheres filling cell volume
        crn:    Job 2 tests core level spectroscopy with core hole (--cls:5,0,1)
         cs:    Input file with large MT radius and use of kmxv to stabilize overlapping FA density
            ... the following are not calculated in the checks fp/test/test.fp --all
       tio2:    high-lying local orbitals, lattice relaxation of relatively complex structure
         na:    compare conventional and extended Na 2p orbitals
     coptso:    test SO coupling, including site-resolved contributions and scaling of L.S
         er:    LDA+U for a metal, and in conjunction with SO coupling.  Also tests HAM_V0BETA
      kfese:    an Fe based superconductor in the striped AFM phase
         ni:    global external field, fixed-spin moment method
        bzt:    generates charge density of Ba3ZnTa2O9 on a mesh a plane, and draws contour plot of it

  jobs:   1: Basic check of programs lmfa,lmf
          2: Core-level spectroscopy (EELS), Mulliken analysis, partial DOS
          3: Check of miscellaneous special features
          4: Different implementations of spin-orbit coupling

EOF
exit

# ---------------- usage: --------------
usage:
cat <<EOF
 usage: test.fp [switches] [file-extension] [testcase-list]
        e.g., "test.fp copt 1"
        If file-extension is missing, test.fp uses copt
        Switches:
        --list       list available materials and tests (no tests are made)
        --no-iactive runs tests without prompting user
        --quiet      runs tests with minimal output and without prompting user
        --all        run through a default list of test cases
        --verb       verbose mode
        --noplot     skip any steps that generate a plot'
        --clean      clean up files generated by this script
        --add0       add suppressed or leading zeros in output for real numbers \`.nnn'
        --poszer     strips (-) sign from numbers represented as 0
        --libxc      Use PBE functional (requires libxc to be linked with code)
                     Applies to certain tests, job 1 only
        --whichexec  prints out which lmf executable it finds in path and exits
        --MPIK       Runs lmf in MPIK mode
                     mpirun -n # lmf ...
        --MPI=#      Runs lmf in MPI mode
                     Note: script mpix must be present, to run the MPI job as
                     mpirun -n # lmf-MPI ...

EOF
exit -1
