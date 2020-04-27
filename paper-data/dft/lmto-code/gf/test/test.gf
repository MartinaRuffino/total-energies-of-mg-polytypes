#!/bin/tcsh -f

# A shell script testing operation of Greens function lmgf
# set verbose

alias call 'set retcall = \!\!:2 ; set callarg = \!\!:3 ; goto \!\!:1'
alias runjob 'set retcall = \!\!:1; set outfile = \!\!:2 ; set callarg = \!\!:3 ; goto runjob'
alias runrdcmd 'set retcall = \!\!:1; set rdcmdfmt = \!\!:2 ; set outfile = \!\!:3 ; set callarg = \!\!:4 ; goto runrdcmd'
alias findcmd  'set retcall = \!\!:1 ; set prog_cmd = \!\!:2 ; set path_name = \!\!:3 ; set make_path = \!\!:4 ; goto findcmd'
alias compare_res 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set refvar = \!\!:4 ; set tol = \!\!:5 ; set passvar = \!\!:6 ; goto compare_res'
alias compare_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set tol = \!\!:4 ; set passvar = \!\!:5 ; goto compare_res_0'
alias get_resf 'set retcall = \!\!:1; set testvar = \!\!:2 ; set refvar = \!\!:3 ; set keyword = \!\!:4  ; set arg_number = \!\!:5 ; set occur_number = \!\!:6 ; set sed_strn = \!\!:7 ; goto get_resf'
alias cnvt_d_fmt  'set retcall = \!\!:1; set testvar = \!\!:2 ; set testval = \!\!:3 ; goto cnvt_d_fmt'
alias zcmpmfiles_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set tol = \!\!:3 ; set passvar = \!\!:4 ; set ndig = \!\!:5 ; set srcfile = \!\!:6 ; set reffile = \!\!:7 ; goto zcmpmfiles_res_0 '
alias zcmpmfiles_res_tol 'set retcall = \!\!:1; set keyword = \!\!:2 ; set tol = \!\!:3 ; set passvar = \!\!:4 ; set ndig = \!\!:5 ; set srcfile = \!\!:6 ; set reffile = \!\!:7 set nlines = \!\!:8; goto zcmpmfiles_res_tol '
alias zcmpmfiles_res_mc 'set retcall = \!\!:1; set keyword = \!\!:2 ; set tol = \!\!:3 ; set passvar = \!\!:4 ; set ndig = \!\!:5 ; set srcfile = \!\!:6 ; set reffile = \!\!:7; set nlines = \!\!:8;; set count = \!\!:9; goto zcmpmfiles_res_mc'
alias zcmpmcx 'set retcall = \!\!:1; set keyword = \!\!:2 ; set toldif = \!\!:3 ; set ndiff = \!\!:4 ; set prec = \!\!:5 ; set passvar = \!\!:6 ;set srcfile = \!\!:7 ; set reffile = \!\!:8 ; goto zcmpmcx '
alias query 'set retcall = \!\!:1 ; set retcall2 = \!\!:2 ; set callarg = \!\!:3 ; goto query'

set allargs = ($argv)

set a
set slow
set testfile = $0
set testdir = $testfile:h
set topdir  = `cd $testdir/../..; pwd`
set tmpdir = $cwd
set space = '        '
set failed = 0
alias zcat 'gunzip -c'
alias zcat 'cat'
set mcexcl = "~excl=CPU~excl=cpudel"

#alias mpix mpirun
set nkj = 12
set nmpi = 4
set jobspassed

# Prepend current working-directory, top-level and related dir to path
set path = ($cwd $topdir $topdir/utils $topdir/testing $path)

# Look for mcx matrix calculator program
set plot = `which fplot`
if (-x "$plot") then
  if `$plot -h | sed -n 1,1p | awk '{print ($1 == "fplot")}'` set have_fplot
endif
set mcx = `which mcx`
if (-x "$mcx") then
  if `$mcx --h |& grep mc | grep vsn | awk '{print ($7 == "(vsn" && ($8 * 1 >= 1.04))}'` set have_mc
endif
# see if ghostscript is available
set gs = `which gs`
if (-x "$gs") then
  if `$gs --help | sed -n 1,1p | awk '{print ($2 == "Ghostscript")}'` set have_ghostscript
endif

set allargs = "$argv"
if (`echo $argv | awk '{i=0; j=0; while (i++<NF) {if ($i == "--all") {j = 1}}; print j}'` == 1) goto all

# --- Pick off switches ---
while (`echo $1 | sed -e 's/\(.\).*/\1/' `  ==  "-")

  set arg1 = $1; shift
  if ($?verb) echo test.gf: parsing switch $arg1
  switch ($arg1)
    case "--quiet":
      set quiet
      unset slow
      breaksw
    case "--add0":
      set add0 = `which add0`
      if (! -x "$add0") then
        echo "test.gf (abort): missing add0"
        exit
      endif
      breaksw
    case "--poszer":
      set poszer = `which poszer`
      if (! -x "$poszer") then
        echo "test.gf (abort): missing poszer"
        exit
      endif
      breaksw
    case "--MPIK=*":
     set narg = `echo $arg1 | sed s/--MPIK=// | awk '{print split($0, a, "," )}'`
     if ($narg < 1 || $narg > 1) then
       echo 'lmgw (abort): bad argument list to MPI=..'
       exit -1
     endif
#    extract nmpi from argument
     set nmpi = `echo $arg1 | sed s/--MPIK=// | awk '{split($0, a, "," ); print a[1]}'`
#    check to ensure argument is a valid integer
     @ nmpi = $nmpi
    case "--MPIK":
      set MPIK
      breaksw
    case "--clean":
      set clean
      breaksw
    case "--noexec":
      set noexec
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
      findcmd chk00 lmgf "$path" "$topdir"
      chk00:
      exit 0
      breaksw
    case "--noplot":
      set noplot
      set have_pldos
      unset have_pldos
      set have_fplot
      unset have_fplot
      breaksw
    case "--gamma=*":
      set gamma=`echo $arg1 | sed s/--gamma=//`
      breaksw
    case "--idxdn":
      set idx
      breaksw
    case "--codecheck":
      set codecheck
      breaksw
    case "--fe2":
      set fe2
      breaksw
    case "--verb*":
    case "-verb*":
      set verb = 1
      breaksw
    case "--all":
      goto all
    default:
      echo unrecognized switch $arg1
      goto usage
  endsw

end

echo ' '
echo "         ---- test.gf: test ASA Green's function program lmgf ---"

# --- Use co as default in the absence of specific choice ---
if ($#argv == 0) then
  set ext = co
  echo "$space .... no file extension specified; use test case " $ext
else
  set ext = $argv[1]
  shift
endif

if (! -e $testdir/ctrl.$ext) then
   echo ' '
   echo " test.gf aborting ... missing file $testdir/ctrl.$ext"
   goto usage
endif

set nkj = 12

if ($ext == "mnpt") then
  echo "$space Case MnPt: tests f orbitals, downfolding, density matrix, spin-polarization, ineq. sphere sizes ..."
  set cplst = ( $testdir/{ctrl}.mnpt )
  set lmargs1 = (-vccor=f -vgfmod=1 -vnk=10 -vldmat=0)
  set gfargs2 = (-vccor=f -vgfmod=1 -vnk=10 -vldmat=0 -vgamrep=t)
  set lmargs3 = (-vnit=10 -vccor=f -vgfmod=1 -vnk=8)
  set lmargs4 = (-vnit=10 -vccor=f -vgfmod=1 -vnk=6)
  set gfargs4 = (-vnit=10 -vccor=f -vgfmod=1 -vnk=6 -vidx=2 -vidxp=3)
  set lmargs5 = "-vnk=$nkj -vbzj=0 -vccor=f -vnit=1 -vgamrep=t"
  set lmargs6 = "-vnit=12 -vscr=223 -vccor=f -vgfmod=1 -vnk=4 -vnl=3 -vgamrep=f"
  if ($?gamma) then   # this has been tested for gamma=0 and gamma=2
    set gfargs2 = (-vccor=f -vgfmod=1 -vnk=10 -vldmat=0 -vgamrep=$gamma)
    set gfargs4 = (-vnit=10 -vccor=f -vgfmod=1 -vnk=6 -vidx=2 -vidxp=3 -vgamrep=$gamma)
    set gfargs5 = "-vnk=$nkj -vbzj=0 -vccor=f -vnit=1 -vidx=2 -vidxp=1 --pr40 -vgamrep=2"
    set gfargs6 = (-vnit=12 -vscr=223 -vccor=f -vgfmod=1 -vnk=4 -vnl=3 -vgamrep=$gamma)
  else  # Gamma as the default
    set gfargs2 = (-vccor=f -vgfmod=1 -vnk=10 -vldmat=0 -vgamrep=t)
    set gfargs4 = (-vnit=10 -vccor=f -vgfmod=1 -vnk=6 -vidx=2 -vidxp=3)
    set gfargs5 = "-vnk=$nkj -vbzj=0 -vccor=f -vnit=1 -vidx=2 --pr40 -vgamrep=2"
    set gfargs5 = "-vnk=$nkj -vbzj=0 -vccor=f -vnit=1 -vidx=2 -vidxp=1 --pr40 -vgamrep=2"
    set gfargs6 = "-vnit=12 -vscr=223 -vccor=f -vgfmod=1 -vnk=4 -vnl=3 -vgamrep=0"
  endif
  if ($?idx) then
    set gfargs2 = ($gfargs2 -vidx=2)
    set gfargs6 = ($gfargs6 -vidx=2)
  endif
  set joblist = ($argv);
  if ( $#joblist == 0 ) set joblist = (2 3 4 5 6)
  set dqtol4  = 5e-4
  set dqtol2  = 2e-4
  set detol2  = 3e-4
  set rmlst2 = ({ctrl,log,mixm,mn,moms,pt,sdot,str,vshft,wkp}.mnpt)
  set rmlst3 = ({ctrl,dmat,log,mixm,mn,pt,sdot,str,vshft}.mnpt)
  set rmlst4 = ({ctrl,log,mixm,mn,pt,save,sdot,str,sv,vshft}.mnpt)
  set rmlst5 = ({ctrl,sstiff,str,sdot,psta,mn,pt,vshft,log,mixm,sv,save,jr,moms,dmat,wkp,syml,bnds}.$ext psta.mnpt.init out.lm.job1 out.lm.job10 out.lm.job11 ps.dat)
  set rmlst6 = ({ctrl,log,mixm,mn,psta,pt,save,sdot,str,sv,vshft}.mnpt psta.mnpt.init)
else if ($ext == "fe2b") then
  set cplst   = $testdir/fe2b/SOC/input-files.tar
  set cplstsf = $testdir/fe2b/SpecFun/input-files.tar
  set lmargsa = "-vso=1 -vnit=1 -vxco=0.4 -vnk=16 -vbeta=0.0 --pr35,25 --rs=1,0"
  set gfargsa = "$lmargsa"
  set rmlsta = "all"
  set rmlstb = "all"
  set lmargsb = "-vso=0 --band:fn=syml -vnit=1 -vxco=0.4 -vnk=32 -vbeta=0.0 fe2b --pr45,25"
  set gfargsb = "$lmargsb"
  set ckspfr = '   -0.300000    1.801371   19.116031   18.625342'
  if ($?codecheck) then
    set gfargsb = "-vcodecheck=t -vso=0 --band:fn=syml -vnit=1 -vxco=0.4 -vnk=8 -vbeta=0.0 fe2b --pr45,25"
    set ckspfr = '   -0.300000    1.801371   19.586419   18.848067'
  endif
else if ($ext == "gas") then
  echo '         Case GaAs: an insulator'
  set lmargs4 = (-vgamma=t)   # Always gamma
# set filel = ({ctrl,as,e1,e2,ga,log,mixm,moms,out,plot,save,sdot,str,sv,syml,bnds,vshft,wkp}.$ext)
  set cplst = ($testdir/ctrl.$ext)
  set dqtol4  = 1e-6
  set rmlst4 = ({as,ctrl,e1,e2,ga,log,mixm,save,sdot,str,sv,vshft}.gas)
else if ($ext == "mnn") then
  echo "$space Case mnn tests f orbitals, downfolding, gamma-rep, large-memory, hexagaonal symmetry..."
  set cplst = ($testdir/ctrl.$ext)
  set lmargs1 = (-vzb=f -vccor=f -vgfmod=1 -vnk=8 -vgamma=t)
  set gfargs1 = (-vzb=f -vccor=f -vgfmod=1 -vnk=8 -vidx=1 -vgamma=t)
  set gfargs2 = (-vzb=f -vccor=f -vgfmod=1 -vnk=8 -vidx=1 --mem -vgamma=t)
  set lmargs4 = (-vzb=f -vccor=f -vgfmod=1 -vnk=8 -vidx=2 -vgamma=t)
  set lmargs5 = (-vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vnk=8 -vidx=2 -vidx2=idx)
  set gfargs5 = (-vnk=12 -vbzj=0 -vccor=f -vnk=8 -vidx=2 -vidx2=2 -vgamma=1 --sites:pair:1,3,5,7)
  if ($?gamma) then
    set gfargs5 = (-vnk=12 -vbzj=0 -vccor=f -vnk=8 -vidx=2 -vidx2=2 -vgamma=$gamma --sites:pair:1,3,5,7)
  endif
# set lmargs5 = (-vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vnk=8 -vidx=2 -vidx2=3)
# set gfargs5 = (-vnk=12 -vbzj=0 -vccor=f -vnk=8 -vidx=2 -vidx2=3 -vgamma=t)
  set lmargs6 = (-vnit=15 -vscr=223 -vzb=f -vccor=f -vgfmod=1 -vnk=6 -vidx=2 -vbzj=0 -vgamma=t)
  set dqtol1  = 8e-4
  set detol1  = 2e-4
  set detol2  = 8e-3
  set dqtol2  = .003
  set dqtol4  = 3e-4
  set detol6  = 2e-5
#   set lmargs1 = (-vccor=f -vgfmod=1 -vnk=10 -vnk2=4 -vldmat=0)
#   set lmargs3 = (-vnit=10 -vccor=f -vgfmod=1 -vnk=6)
  set joblist = ($argv);
  if ( $#joblist == 0 ) set joblist = (1 2 4 5 6)
  set rmlst1 = {a1,c1,ctrl,ea1,ec1,log,mixm,moms,sdot,str,vshft,wkp}.$ext
  set rmlst2 = {a1,c1,ctrl,ea1,ec1,log,mixm,moms,sdot,str,vshft,wkp}.$ext
  set rmlst4 = ({a1,c1,ctrl,ea1,ec1,log,mixm,save,sdot,str,sv,vshft}.$ext)
  set rmlst5 = ({ctrl,sstiff,str,sdot,moms,wkp,c1,a1,ec1,ea1,vshft,log,mixm,sv,save,psta,jr}.$ext psta.mnn.init)
  set rmlst6 = ({a1,c1,ctrl,ea1,ec1,log,mixm,psta,save,sdot,str,sv,vshft}.mnn psta.mnn.init)

else if ($ext == "bi2te3") then
  set lmargs7; set twoc=0
  set rmlst7 = ({a1,c1,ctrl,ea1,ec1,log,mixm,psta,save,sdot,str,sv,vshft}.$ext)
  set cplst = ($testdir/{ctrl,rsta,syml,vshft}.$ext)
  set dqtol7 = .004
  set detol7 = .38

else if ($ext == "nife") then
  echo "$space Case nife tests noncollinearity, linear response, spin statics ..."
  set cplst = ($testdir/{ctrl,eula}.$ext $testdir/eref.dat)
  set touchl = a1.$ext
# set filel = ({ctrl,eula,str,sdot,a\*,log,mixm,ma,vshft,sv,save,psta,dos,evec,gfqp}.$ext psta.nife.init eref.dat)
  set lmargs6 = "-vnit=27 -vscr=882 --spinoronly -vsforce=1 -vsdmod=11 -vfscal=15 -vnk=2"
  set dptol6 = 1e-6
  set dqtol6  = 1e-3
  set dmtol6  = .05
  set detol6  = 6e-4
  set joblist = ($argv);
  if ( $#joblist == 0 ) set joblist = (6)
  set rmlst6 = ({a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a1,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a2,a30,a31,a32,a3,a4,a5,a6,a7,a8,a9,ctrl,dos,eula,gfqp,log,ma,mixm,psta,save,sdot,str,sv,vshft}.nife psta.nife.init eref.dat)
else if ($ext == "fepd") then
  echo "$space Case FePd tests disordered local moments ..."
  set cplst = ($testdir/fepd/*.$ext)
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (9)
  set lmargs9 = "-vnit=1 -vbeta=0.05 -vnk=12 --pr45,25"
  if ($?gamma) then
    set gfargs9 = "$lmargs9 -vgamma=$gamma"
  endif
  set rmlst9 = ({ctrl,shfac,fe\#1,fe\#2,fe\#3,fe\#4,fe\#5,fe\#6,fe\#7,fe\#8,fe\#9,fe\#10,fe\#11,fe\#12,fe\#13,fe\#14,fe\#15,fe\#16,fe\#17,fe\#18,fe\#19,fe\#20,fe\#21,fe\#22,fe\#23,fe\#24,fe\#25,fe\#26,fe\#27,fe\#28,fe\#29,fe\#30,fe\#31,fe\#32,fe,log,mixm,omega1,om-hr1,pd,save,sdot,str,sv,vshft}.$ext)
else if ($ext == "fev" || $ext == "fevg") then
  echo "$space Case FeV tests chemical CPA ..."
  set cplst = ($testdir/{ctrl,fe\#1,fe\#2,omega1,vshft}.$ext)
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (5 10)
  set lmargsa = "-vnit=1 -vxco=0.01 -vbeta=0.1 -vnk=16 --pr45,25"
  set gfargsa = "$lmargsa"
  set rmlsta = ({shfac,ctrl,fe\#1,fe\#2,fe,log,mixm,omega1,om-hr1,save,sdot,str,sv,vshft}.$ext)
  set gfargs5 = "-vxco=0.01 -vbeta=0.1 -vnk=16 --pr45,25"
  if ($?gamma) then
    set gfargs5 = "-vxco=0.01 -vbeta=0.1 -vnk=16 --pr45,25 -vgamma=$gamma"
    set gfargsa = "$lmargsa -vgamma=$gamma"
  endif
  set nkj=16 fermi=0
  set rmlst5 = ({ctrl,sstiff,syml,wkp,moms,str,sdot,a,vshft,log,mixm,sv,save,bnds,jr,eula,tcqres,shfac}.$ext)
 else if ($ext == "femnpt") then
  echo "$space Case FeMnPt tests chemical CPA ..."
  set cplst = ($testdir/{ctrl,emesh,fe\#1,fe\#2,pt2,pt3,pt4,pt,vshft}.$ext)
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (5 10)
  set nkj = 6
  set gfargs5 = "-vnit=1 --pr40 -vnk=6 --no-iactiv -vnz=21 -ef=0"
  set gfargs5x = "-vnit=1 --pr40,30 -vnk=6 --no-iactiv -vnz=21 -ef=0 -vgfmod=1 --quit=rho"
  set lmargsa = "-vnit=1 --pr40"
  set gfargsa = "$lmargsa -vnz=31"
  set rmlst5 = ({mixm,shfac,ctrl,emesh,fe\#1,fe\#2,fe,jr,log,omega1,omega2,omega5,omega6,pt2,pt3,pt4,pt,sdot,str,vshft}.$ext)
  set rmlsta = ({shfac,ctrl,emesh,fe\#1,fe\#2,fe,log,mixm,omega1,omega2,omega5,omega6,om-hr1,om-hr2,om-hr5,om-hr6,pt2,pt3,pt4,pt,save,sdot,str,sv,vshft.femnpt}.$ext)
  set fermi=0
 else if ($ext == "co") then
  echo '         Case co : a hexagonal environment with equivalent atoms.'
  set cplst = ($testdir/{ctrl.co,a.co,q.bulk,syml.co,occnum.co})
  set cpagain
# set filel = ({ctrl,moms,str,sdot,a,vshft,log,mixm,ma,sv,save,dmat,jr,wkp,syml,bnds,occnum,dmats}.$ext q.bulk)
  set lmargs1 = "-vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vldmat=0 -vcbya=1.625"
  set gfargs1 = "-vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vldmat=0 -vcbya=1.625 -vgamrep=f"
  if ($?gamma) then
    set gfargs1 = "-vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vldmat=0 -vcbya=1.625 -vgamrep=$gamma"
  endif
  set lmargs3 = "-vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vcbya=1.625"
  set lmargs4 = "-vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vcbya=1.625"
  set lmargs5 = "-vnk=12 -vbzj=0 -vccor=f -vcbya=1.625"
  set ldauargs5 = "-vldau=1 -vudiag=1"
  set lmargs7 = "-vnc=1 -vtwoc=0 -ef=-0.025725 -vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vldmat=0 -vcbya=1.625"
  set gfargs7 = "-vnc=1 -vtwoc=0 -ef=-0.025725 -vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vldmat=0 -vcbya=1.625"
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (1 2 3 4 5 7)
  set detol2  = 1e-3
  set detol7  = 1e-3
  set rmlst1 = ({a,ctrl,log,mixm,moms,occnum,sdot,str,syml,vshft,wkp}.co q.bulk)
  set rmlst2 = ($rmlst1)
  set rmlst3 = ({a,ctrl,dmat,log,mixm,occnum,sdot,str,syml,vshft}.co q.bulk)
  set rmlst4 = ({a,ctrl,log,mixm,occnum,save,sdot,str,sv,syml,vshft}.co q.bulk)
  set rmlst5 = ({dmats,a,bnds,ctrl,jr,log,moms,occnum,qpp,sdot,sstiff,str,syml,vshft,wkp}.co q.bulk)
  set rmlst7 = ({a,ctrl,eula,log,mixm,moms,occnum,save,sdot,str,sv,syml,vshft,wkp}.co)
else if ($ext == "eras") then
  echo "$space Case ErAs: a system with partially filled f band."
  set cplst = ($testdir/{ctrl,dmats,rsta}.eras)
  set lmargs1 = "--rs=1,0 -vnk=12 -vudiag=0 -vsharm=1"
  set lmargs4 = "--rs=1,0 -vnk=6 -vudiag=0"
  set lmargs7 = "--rs=1,0 -vnk=12 -vudiag=1 -vsharm=1 -vnc=1"
  set gfargs7 = "--rs=1,0 -vnk=12 -vudiag=1 -vsharm=1 -vnc=1 --spinoronly"
  set dqtol1 = 2.6e-4
#    set detol1  = 0.006
  set dqtol2  = 8e-4
  set detol2  = .02
  set dqtol4  = 7e-3
  set detol4  = 2.7e-4
  set dqtol7 = .002 stdotol7 = 5e-5
  set detol7  = .04
  set cpagain
  set ldau
  set rmlst1 = ({ctrl,as,dmat,dmats,ea1,ec1,er,log,mixm,qpp,rsta,save,sdot,str,sv,vshft,dos,moms,wkp}.eras)
  set rmlst2 = ( $rmlst1 )
  set rmlst4 = ({as,ctrl,dmat,dmats,ea1,ec1,er,log,mixm,qpp,rsta,save,sdot,str,sv,vshft}.eras)
  set rmlst7 = ({ctrl,as,dmat,dmats,ea1,ec1,er,log,mixm,qpp,rsta,save,sdot,str,sv,vshft,dos,eula.moms,wkp,eula,moms}.eras)
else if ($ext == "fccfe") then
  echo '         Case fcc Fe: a noncollinear system.'
  set lmargs1 = "-vnk=16 -vccor=f -vsdyn=t"
  set lmargs7 = "-vnk=10 -vccor=f -vsdyn=t"
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (7)
  set rmlst7 = ({ctrl,a2,a3,a4,a,dmat,eula,log,moms,mq,sdot,str,vshft,wkp,bnds}.fccfe)
else if ($ext == "fe") then
  set nkj=16
  set lmargs4 = (-vbxcs=1 -vnk=$nkj -vccor=f -vtwoc=f -vfile=f)
  set lmargs5 = "-vnk=$nkj -vccor=f -vtwoc=f -vfile=f"
  if ($?fe2) set lmargs5 = "-vnk=$nkj -vccor=f -vtwoc=f -vfile=t"
  set ldauargs5 = "-vldau=1 -vudiag=1"
  set vextargs5 = "-vrdvext=1"
# set filel = ({a,vshft,log,mixm,sv,save,bnds,dmats,occnum,vext}.$ext)
  set cplst = ($testdir/{ctrl.$ext,a.$ext,syml.$ext,occnum.$ext,vext.$ext,shfac.$ext})
  set rmlst4 = (mixm.$ext vshft.$ext dmats.$ext tcqres.$ext bxc.$ext)
  set rmlst5 = (mixm.$ext vshft.$ext dmats.$ext tcqres.$ext bxc.$ext)
else if ($ext == "ni") then
  echo '         Case fcc Ni: a simple one-atom/cell.'
  set nkj=16
  set lmargs5 = "-vnk=$nkj -vbzj=0 -vccor=f -vtwoc=f"
# set filel = ({ctrl,sstiff,syml,wkp,moms,str,sdot,a,vshft,log,mixm,sv,save,bnds,jr}.$ext)
  set cplst = ($testdir/{ctrl.$ext,a.$ext,syml.$ext})
  set joblist = ($argv);
  if ( $#joblist == 0 ) set joblist = (5)
  set rmlst5 = ({ctrl,sstiff,syml,wkp,moms,str,sdot,a,vshft,log,mixm,sv,save,bnds,jr,eula,tcqres}.$ext)
#   set rmlst5 = ({a,bnds,ctrl,jr,log,moms,out.job10  out.job11  sdot,sstiff,str,syml,vshft,wkp}.$ext)
else if ($ext == "cdte") then
  echo '         Case CdTe: an insulator'
  set lmargs8 = "-vso=f -vbzj=0 -vccor=t -vtwoc=f -vnk=6 -vnk2=nk -vgamma=t -vadnf=f -vnl=4 -vlmxb=3 -vef=-.19"
  set lmargs8 = "-vso=f -vbzj=0 -vccor=t -vtwoc=f -vnk=6 -vnk2=nk -vsig=1 -vgamma=t -vadnf=f -vnl=4 -vlmxb=3 -vef=-.19"
  if ($?gamma) then
    set gfargs8 = "-vso=f -vbzj=0 -vccor=t -vtwoc=f -vnk=6 -vnk2=nk -vsig=1 -vgamma=$gamma -vadnf=f -vnl=4 -vlmxb=3 -vef=-.19"
  endif
# set filel = ({ctrl,c1,a1,ec1,ea1,syml,vshft,log,sigm,mixm,sv,save,wkp,sgm,gfdm,sigm-old,p0,gfdm1,str,sdot,bnds,qasa}.$ext qasa semi.mater)
  set cplst = ( $testdir/{{ctrl,syml,qasa}.$ext,semi.mater} )
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (8)
  set rmlst8 = ({a1,bnds,c1,ctrl,ea1,ec1,gfdm1,log,p0,qasa,save,sdot,sgm,sigm,str,syml,vshft,wkp}.cdte qasa semi.mater specialspecc)
else if ($ext == "fept2") then
  echo '         Case FePt2: A noncollinear GF calculation with SO coupling'
  set rmlst7 = ({ctrl,fe,log,mixm,pt,save,sdot,str,sv,vshft,vshft-bk}.$ext)
  set lmargs7 = "-vnk=12 -vso=1 -veula=1 --quit=rho"
  set gfargs7 = "-vnit=-1 -vnk=12 -vso=1 -veula=1 --quit=rho"
  if ($?gamma) then
    set gfargs7 = "$gfargs7 -vgamma=$gamma"
  endif
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (7)
  set cplst = $testdir/{ctrl,vshft}.$ext
  set twoc = 1
  set dttol7 = 5e-7  # tolerance in torque
  set dqtol7 = 2e-6  # tolerance in RMS DQ
  set dmtol7 = 1e-5  # Tolerance in the moment
  set detol7 = 1e-5  # Tolerance in ehf
else if ($ext == "fept") then
  echo '         Case FePt: A noncollinear GF calculation with downfolding'
  set nkj=8
  set rmlst1 = ({ctrl,eula,fe,jr,log,mixm,moms,pt,rsta,save,sdot,str,sv,vshft,wkp}.fept)
  set rmlst2 = ({ctrl,eula,fe,jr,log,mixm,moms,pt,rsta,save,sdot,str,sv,vshft,wkp}.fept)
  set rmlst5 = ({ctrl,sstiff,syml,wkp,moms,str,sdot,vshft,log,mixm,sv,save,bnds,jr,eula,tcqres,gfqp,fe,pt,rsta}.$ext)
  set lmargs1 = "-vnk=12 -vnsp=2 -vso=0 -vidxp=1 -vidxf=1 -vccor=f --rs=1,0"
  set lmargs2 = "-vnk=12 -vnsp=2 -vso=0 -vidxp=1 -vidxf=1 --rs=1,0"
  set lmargs7 = "-vnk=8 -vnit=0 -vnc=1 -vntht=6 --rs=1,0"
  if ($?gamma) then
    set lmargs1 = "-vnk=12 -vnsp=2 -vso=0 -vidxp=1 -vidxf=1 -vccor=f --rs=1,0 -vgamma=$gamma"
    set lmargs2 = "-vnk=12 -vnsp=2 -vso=0 -vidxp=1 -vidxf=1 --rs=1,0 -vgamma=$gamma"
    set lmargs7 = "-vnk=8 -vnit=0 -vnc=1 -vntht=6 --rs=1,0 -vgamma=$gamma"
  endif
  set lmargs5 = "-vnk=8 -vnsp=2 --rs=1,0 -vso=1 -vidxp=1 -vidxf=1 -vsharm=1"
  set joblist = ($argv)
  if ( $#joblist == 0 ) set joblist = (5 7)
  set cplst = ( $testdir/{ctrl.$ext,rsta.fept.lm.ntht=0,rsta.fept.lm.ntht=6,rsta.fept.lmgf.ntht=0,rsta.fept.lmgf.ntht=6,rsta.fept} )


  touch xxxx1111.$ext
  set rmlst6 = (`echo *.$ext | sed s/xxxx1111.$ext//`)
  set rmlst7 = (`echo *.$ext | sed s/xxxx1111.$ext//` {rsta.fept.lm.ntht=0,rsta.fept.lm.ntht=6,rsta.fept.lmgf.ntht=0,rsta.fept.lmgf.ntht=6} out.idxp2.lmgf)
  rm -f xxxx1111.$ext
  set dqtol2 = .012
  set dqtol7 = 1e-6  stdotol7 = 1e-3
  set dmtol7 = 1e-5
  set detol2 = 0.07
  set detol7 = 0.08
else
  echo test.gf: No test case for $ext
  exit -1
endif

if ( $?joblist == 0 ) then
set joblist = ($argv)
if ( $#joblist == 0 ) set joblist = (1 2 3 4 5 6 7 8 9 10 11)
endif

# echo $joblist | grep 1 >/dev/null
echo $joblist | egrep '\b1\b' >/dev/null
if ($status) goto chk1e
cat <<EOF

         --- Test 1.  Compare band sum, Fermi level to program lm, 2rd order H ---
         Tests output integrated properties of lmgf.
         Compare some results of lmgf (sumev, Fermi level, output q) to program lm.
         lmgf and lm 2C approximation produce identical results for a sufficiently fine k-mesh.

EOF
if ("$ext" == co) then
cat <<EOF
         Co, hcp:
         The axial ratio is deliberately set small to highlight
         off-diagonal elements in the density matrix.

EOF
else if ("$ext" == fept) then
cat <<EOF
         lm and lmgf Fermi levels should differ by Delta Ef = 0.002 Ry, while Delta sumev = 9.3e-05
         If 32 k divisions are used, Delta Ef = 0.000066 Ry and Delta sumev = 3.8e-5

EOF
else if ("$ext" == mnn) then
cat <<EOF
         For the MnN case, there is a shift to the gamma-representation.
         A shift to the gamma-representation is necessary owing to a pole
         in the d channel of large empty sphere.

EOF
else if ("$ext" == eras) then
cat <<EOF
         The ErAs case tests the LDA+U Green's function.
         Spherical harmonics are used, to make the density-matrix nearly diagonal.

         Independent checks you can make:
          - The dependence on the off-diagonal part of U (repeat with UDIAG=1)
            Verify energy, magnetic moment, density-matrix and bands are very similar.
            Note also that the lmgf and lm RMS change in (output-input) density-matrix track each other.
            (You can find reference files in $testdir/*udiag)

          - How closely the calculation is rotationally invariant (repeat without spherical harmonics, SHARM=0)
            Verify Fermi energy, band sum, magnetic moment, are very similar, but not the output moments.
            The charge is well given, but the first and second power moments are not.
            (You can find reference files in $testdir/*rharm)

EOF
endif

set refoutlm=$testdir/out.lm.$ext.2c refoutg=$testdir/out.lmgf.$ext.2c testoutlm=out.lm testoutg=out.lmgf
if ($?lmargs1 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk1e
endif
if ($?gfargs1 == 0) set gfargs1 = ($lmargs1)
if (! -e $refoutlm) then
  echo "$space ... skipping test : missing reference file $refoutlm"
  goto chk1e
endif
if (! -e $refoutg) then
  echo "$space ... skipping test : missing reference file $refoutg"
  goto chk1e
endif
set pass; set ldaupass ; unset ldaupass; set vextpass ; unset vextpass

query chk11 chk1e 'run this test'
chk11:
# ... Look for executables
findcmd chk11a rdcmd "$path" "$topdir"
chk11a:
findcmd chk11b lm "$path" "$topdir"
chk11b:
findcmd chk11c lmstr "$path" "$topdir"
chk11c:
findcmd chk11d lmgf "$path" "$topdir"
chk11d:
set lmgf = lmgf
if ($?MPIK) then
  set lmgf = "mpirun -n $nmpi lmgf"
endif

set twoc = 1
set gmaxdif = 0

if ($?haveout) goto chk13r

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f $rmlst1"
             rm -f $rmlst1
if ($?clean) then
  echo "$space rm -f $testoutlm $testoutg"
               rm -f $testoutlm $testoutg
  goto chk1e
endif
echo "$space cp $cplst ."
             cp $cplst .
runjob chk12 /dev/null "$lmstr --no-iactiv $ext -vnit=0 $lmargs1"
chk12:
runjob chk13 $testoutg "$lm $ext -vnit=0 -vtwoc=$twoc $lmargs1 --no-iactiv"
chk13:
echo ' '
echo "$space ... Run lm band pass to find Fermi level and band sum :"
echo "$space $lm $ext -vtwoc=$twoc -vnit=1 $lmargs1 --iactiv=no --quit=rho >$testoutlm"
             $lm $ext -vtwoc=$twoc -vnit=1 $lmargs1 --iactiv=no --quit=rho >$testoutlm
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lm returned successfully; output is in file $testoutlm"
#  set fermi = `grep '  ef ' log.$ext | tail -1 | awk '{print $7}'`
#  set sumev = `cat log.$ext | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`

chk13r:

set testout=$testoutlm refout=$refoutlm
get_resf chk1ca fermi fermiref "Fermi energy" 4 0 ";"
chk1ca:
get_resf chk1cb sumev sumevref "sumev=" 4 0 sumev=
chk1cb:
if ($?ldau > 0) then
get_resf chk1cc rmsdmlm rmsdmlmref "RMS diff in dens mat" 5 0 "mat("
chk1cc:
set rmsdmlm = `echo "$rmsdmlm" | sed 's/)//'`
set rmsdmlmref = `echo "$rmsdmlmref" | sed 's/)//'`
endif
set lmdq  = `grep 'RMS DQ=' log.$ext | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`

if (! $?quiet) then
  echo ' '
  echo "$space Fermi energy from lm               = $fermi"
  echo "$space Fermi energy of reference          = $fermiref"
  set diff = `echo $fermi $fermiref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space band sum from lm                   = $sumev"
  echo "$space band sum of reference              = $sumevref"
  set diff = `echo $sumev $sumevref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space RMS DQ (output-input) from lm      = $lmdq"

  if ($?ldau > 0) then
  echo ' '
  echo "$space RMS change in dmat from lm         = $rmsdmlm"
  echo "$space RMS change of reference            = $rmsdmlmref"
  endif

  echo ' '
else
  echo "$space Found Fermi level = $fermi, sumev = $sumev,  rms DQ between input and output moments = $lmdq"
endif

if ($?haveout) goto chk14r

if ($?add0) then
  echo -n "         ..." ; $add0 $testoutlm
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutlm
endif
call zdiffiles chk14 "CPU 1 $testoutlm $refoutlm"
chk14:

echo ' '
echo "$space ... invoking lmgf using Fermi level output from lm:"
echo "$space rm -f vshft.$ext mixm.$ext"
             rm -f vshft.$ext mixm.$ext
if ($?cpagain) then
echo "$space cp $cplst ."
             cp $cplst .
endif
echo "$space $lmgf $ext -vtwoc=$twoc -vnit=1 $gfargs1 -ef=$fermi --iactiv=no --quit=rho >$testoutg"
             $lmgf $ext -vtwoc=$twoc -vnit=1 $gfargs1 -ef=$fermi --iactiv=no --quit=rho >$testoutg
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully.  The following data were output:"
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutg
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutg
endif

chk14r:

set zval = `grep 'zval' log.$ext | tail -1 | awk '{print $2}'`

set vshftp = `grep 'gfasa  vshft:' log.$ext | head -1 | awk '{print $5}'`
set vshft  = `grep 'gfasa  vshft:' log.$ext | tail -1 | awk '{print $5}'`
#set sumevg = `grep 'sumev=' log.$ext | tail -1 | awk '{print $3}' | sed 's/sumev=//'`
set sumevg = `grep 'sumev=' $testoutg | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set sumevgr = `grep 'sumev=' $refoutg | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`

# set dc = `echo null | awk -v zval=$zval -v vshft=$vshft '{print zval*vshft}'`
# echo DC= $dc
set dc = `grep 'd\.c\. term vbar' $testoutg | awk '{print $NF}'`
set dcr = `grep 'd\.c\. term vbar' $refoutg | awk '{print $NF}'`
# echo DC= $dc

set sumevgdc = `echo $sumevg $dc | awk '{print ($1-$2)}'`
set ediff = `echo $sumevgdc $sumev  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set gfdq  = `grep 'RMS DQ=' log.$ext | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`

if (! $?quiet) then
  echo ' '
  echo "$space Pade estimate for Fermi level shift = $vshftp"
  echo "$space             final Fermi level shift = $vshft"

  echo ' '
  echo "$space band sum from lmgf                  = $sumevg"
  echo "$space band sum of reference               = $sumevgr"
  set diff = `echo $sumevg $sumevgr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference       = $diff"

  echo ' '
  echo "$space d.c. term zval*vshft from lmgf      = $dc"
  echo "$space d.c. term reference                 = $dcr"

  echo ' '
  echo "$space lm band structure energy            = $sumev"
  echo "$space lmgf band sum - d.c.                = $sumevgdc"
  set diff = `echo $sumevgdc $sumev  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference       = $diff"

  echo ' '
  echo "$space RMS DQ (output-input) from lm       = $lmdq"
  echo "$space RMS DQ (output-input) from lmgf     = $gfdq"

endif

set testout=$testoutg refout=$refoutg
if ($?ldau > 0) then
get_resf chk1cd rmsdmlmgf rmsdmlmgfref "RMS diff in dens mat" 5 0 "mat("
chk1cd:
set rmsdmlmgf = `echo "$rmsdmlmgf" | sed 's/)//'`
set rmsdmlmgfref = `echo "$rmsdmlmgfref" | sed 's/)//'`
echo "$space RMS change in dmat from lmgf        = $rmsdmlmgf"
echo "$space RMS change in dmat from lm          = $rmsdmlm"
#  echo "$space RMS change of reference            = $rmsdmlmgfref"
endif

#  if ($?ldau > 0) goto chk15
call qprint chk15 "$space The differences should go away as the qpt and energy mesh are made finer."
chk15:

if ($?add0) then
  echo -n "         ..." ; $add0 $testoutg
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutg
endif

call zdiffiles chk16 "CPU 1 $testoutg $refoutg - $mcexcl~excl=gippad~excl=gfzerq~excl=spin"
chk16:

chk1c:
echo ' '
call qprint chk1ce "$space ... automatic pass checks :"
chk1ce:
if ($?dqtol1 == 0) set dqtol1 = 2.5e-4
echo -n "$space RMS dq ($gfdq) generated by lmgf within tol $dqtol1 of lm ($lmdq)? ... "
if (`echo $lmdq $gfdq $dqtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

# set ediff = .1

if ($?detol1 == 0) set detol1 = .0003
echo -n "$space sumev ($sumevgdc) generated by lmgf within tol $detol1 of lm? ... "
if (`echo $ediff 0 $detol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

if ($?ldau > 0) then
set dmatstol = 1e-4
echo -n "$space RMS change in dmat ($rmsdmlmgf) generated by lmgf within tol $dmatstol of lm ($rmsdmlm)? ... "
if (`echo $rmsdmlmgf $rmsdmlm $dmatstol | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
#zcmpmfiles_res_0 chk1dm "... Max deviation in dmats.$ext from reference" $dmatstol pass 4 dmats.$ext $testdir/dmats.$ext
chk1dm:
endif

if (! $?gmaxdif) then
else if (! $?mcx) then
  echo "$space ... mcx not installed ... no check on stdout"
else
  if (! $?stdotol1) set stdotol1 = 1e-6
  echo -n "$space maximum numerical difference in stdout ($gmaxdif) <= tol ($stdotol1) ? ... "
  if (`echo ' ' | awk -v maxdif=$gmaxdif -v tol=$stdotol1 '{print (maxdif <= tol)}'` == 1) then
    echo yes
  else
    echo no "(max diff = $gmaxdif)"
    unset pass
  endif
endif

if ($?clean) then
else if ($?pass) then
    echo "$space test 1 PASSED ($ext)"
    set jobspassed = ($jobspassed 1)
else
    echo "$space test 1 FAILED ($ext)"
    set failed = ($failed 1)
endif

chk1e:

echo $joblist | grep 2 >/dev/null
if ($status) goto chk2e
cat <<EOF

         --- Test 2.  Compare band sum, Fermi level to program lm, 3rd order H ---
         Tests output integrated properties of lmgf.
         Compare some results of lmgf (sumev, Fermi level, output q) to program lm.
         Similar to Test 1, but uses 3rd order potential functions.
         lmgf and lm the 3C approximation produce similar, but not identical results.

EOF
if ("$ext" == mnpt) then
cat <<EOF
         In the mnpt case, the GF is transformed to the gamma-representation,
         demonstrating that the transformation does not alter the GF.

EOF
endif
if ("$ext" == mnn) then
cat <<EOF
         For the mnn case, it checks downfolding (idx=2).

         A shift to the gamma-representation is necessary owing to a pole
         in the d channel of large empty sphere.

         Also, the large-memory (--mem) switch is checked.

EOF
endif
if ("$ext" == eras) then
cat <<EOF
         The ErAs case tests the LDA+U Green's function.
         Spherical harmonics are used, to make the density-matrix nearly diagonal.

EOF
endif

set refoutlm=$testdir/out.lm.$ext.3c refoutg=$testdir/out.lmgf.$ext.3c testoutlm=out.lm testoutg=out.lmgf
if ($?lmargs1 == 0 && $?lmargs2 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk2e
endif
if ($?lmargs2 == 0) then
  set lmargs2 = ($lmargs1)
endif
if ($?gfargs2 == 0) set gfargs2 = ($lmargs2)
if (! -e $refoutlm) then
  echo "$space ... skipping test : missing reference file $refoutlm"
  goto chk2e
endif
if (! -e $refoutg) then
  echo "$space ... skipping test : missing reference file $refoutg"
  goto chk2e
endif
set pass; set ldaupass ; unset ldaupass; set vextpass ; unset vextpass

query chk21 chk2e 'run this test'
chk21:
# ... Look for executables
findcmd chk21a rdcmd "$path" "$topdir"
chk21a:
findcmd chk21b lm "$path" "$topdir"
chk21b:
findcmd chk21c lmstr "$path" "$topdir"
chk21c:
findcmd chk21d lmgf "$path" "$topdir"
chk21d:
if ($?MPIK) then
  set lmgf = "mpirun -n $nmpi lmgf"
endif

set twoc = 0
set gmaxdif = 0
echo "$space ... set up ASA strux and starting potential"
# echo "$space rm -f $filel"
#              rm -f $filel
echo "$space rm -f $rmlst2"
             rm -f $rmlst2
if ($?clean) then
  echo "$space rm -f $testoutlm $testoutg"
               rm -f $testoutlm $testoutg
  goto chk2e
endif
echo "$space cp $cplst ."
             cp $cplst .
runjob chk22 /dev/null "$lmstr --no-iactiv $ext -vnit=0 $lmargs2"
chk22:
runjob chk23 $testoutg "$lm --no-iactiv $ext -vnit=0 -vtwoc=$twoc $lmargs2"
chk23:
echo ' '
echo "$space ... Run lm band pass to find Fermi level and band sum :"
echo "$space $lm $ext -vtwoc=$twoc -vnit=1 $lmargs2 --iactiv=no --quit=rho >$testoutlm"
             $lm $ext -vtwoc=$twoc -vnit=1 $lmargs2 --iactiv=no --quit=rho >$testoutlm
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lm returned successfully; output is in file $testoutlm"
#  set fermi = `grep '  ef ' log.$ext | tail -1 | awk '{print $7}'`
#  set sumev = `grep '  sev ' log.$ext | tail -1 | awk '{print $9}'`
#  set sumev = `cat log.$ext | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set testout=$testoutlm refout=$refoutlm
get_resf chk2ca fermi fermiref "Fermi energy" 4 0 ";"
chk2ca:
get_resf chk2cb sumev sumevref "sumev=" 4 0 sumev=
chk2cb:
if ($?ldau > 0) then
get_resf chk2cc rmsdmlm rmsdmlmref "RMS diff in dens mat" 5 0 "mat("
chk2cc:
set rmsdmlm = `echo "$rmsdmlm" | sed 's/)//'`
set rmsdmlmref = `echo "$rmsdmlmref" | sed 's/)//'`
endif
grep 'RMS DQ=' log.$ext | tail -1 | awk '{print $9}'  | sed 's/DQ=//' > /dev/null
if ($status == 1) then
  echo abort ... no dq in log file
endif
set lmdq  = `grep 'RMS DQ=' log.$ext | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`
if (! $?quiet) then
  echo ' '
  echo "$space Fermi energy from lm               = $fermi"
  echo "$space Fermi energy of reference          = $fermiref"
  set diff = `echo $fermi $fermiref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space band sum from lm                   = $sumev"
  echo "$space band sum of reference              = $sumevref"
  set diff = `echo $sumev $sumevref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space RMS DQ (output-input) from lm      = $lmdq"

  if ($?ldau > 0) then
  echo ' '
  echo "$space RMS change in dmat from lm         = $rmsdmlm"
  echo "$space RMS change of reference            = $rmsdmlmref"
  endif

  echo ' '
else
  echo "$space Found Fermi level = $fermi, sumev = $sumev,  rms DQ between input and output moments = $lmdq"
endif

if ($?add0) then
  echo -n "         ..." ; $add0 $testoutlm
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutlm
endif
call zdiffiles chk24 "CPU 1 $testoutlm $refoutlm"
chk24:

echo ' '
echo "$space ... invoking lmgf using Fermi level output from lm:"
echo "$space rm -f vshft.$ext mixm.$ext"
             rm -f vshft.$ext mixm.$ext
if ($?cpagain) then
echo "$space cp $cplst ."
             cp $cplst .
endif
echo "$space $lmgf $ext -vtwoc=$twoc -ef=$fermi -vnit=1 $gfargs2 --iactiv=no --quit=rho >$testoutg"
             $lmgf $ext -vtwoc=$twoc -ef=$fermi -vnit=1 $gfargs2 --iactiv=no --quit=rho >$testoutg
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully.  The following data were output:"
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutg
endif

set zval = `grep 'zval' log.$ext | tail -1 | awk '{print $2}'`
set vshftp = `grep 'gfasa  vshft:' log.$ext | head -1 | awk '{print $5}'`
set vshft  = `grep 'gfasa  vshft:' log.$ext | tail -1 | awk '{print $5}'`
#set sumevg = `grep 'sumev=' log.$ext | tail -1 | awk '{print $3}' | sed 's/sumev=//'`
set sumevg = `grep 'sumev=' $testoutg | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set sumevgr = `grep 'sumev=' $refoutg | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`

# set dc = `echo null | awk -v zval=$zval -v vshft=$vshft '{print zval*vshft}'`
# echo DC= $dc
set dc = `grep 'd\.c\. term vbar' $testoutg | awk '{print $NF}'`
set dcr = `grep 'd\.c\. term vbar' $refoutg | awk '{print $NF}'`
# echo DC= $dc

set sumevgdc = `echo $sumevg $dc | awk '{print ($1-$2)}'`
set ediff = `echo $sumevgdc $sumev  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set gfdq  = `grep 'RMS DQ=' log.$ext | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`

if (! $?quiet) then
  echo ' '
  echo "$space Pade estimate for Fermi level shift = $vshftp"
  echo "$space             final Fermi level shift = $vshft"

  echo ' '
  echo "$space band sum from lmgf                  = $sumevg"
  echo "$space band sum of reference               = $sumevgr"
  set diff = `echo $sumevg $sumevgr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference       = $diff"

  echo ' '
  echo "$space d.c. term zval*vshft from lmgf      = $dc"
  echo "$space d.c. term reference                 = $dcr"

  echo ' '
  echo "$space lm band structure energy            = $sumev"
  echo "$space lmgf band sum - d.c.                = $sumevgdc"
  set diff = `echo $sumevgdc $sumev  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference       = $diff"

  echo ' '
  echo "$space RMS DQ (output-input) from lm       = $lmdq"
  echo "$space RMS DQ (output-input) from lmgf     = $gfdq"

endif

set testout=$testoutg refout=$refoutg
if ($?ldau > 0) then
get_resf chk2cd rmsdmlmgf rmsdmlmgfref "RMS diff in dens mat" 5 0 "mat("
chk2cd:
set rmsdmlmgf = `echo "$rmsdmlmgf" | sed 's/)//'`
set rmsdmlmgfref = `echo "$rmsdmlmgfref" | sed 's/)//'`
echo "$space RMS change in dmat from lmgf       = $rmsdmlmgf"
echo "$space RMS change in dmat from lm         = $rmsdmlm"
#  echo "$space RMS change of reference            = $rmsdmlmgfref"
endif

if ($?add0) then
  echo -n "         ..." ; $add0 $testoutg
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutg
endif
# Too difficult to maintain error tolerance across different OS ... just use zdiffiles
call zdiffiles chk25 "CPU 1 $testoutg $refoutg - $mcexcl~excl=gippad~excl=gfzerq~excl=spin"
chk25:


chk2c:
echo ' '
call qprint chk2ce "$space ... automatic pass checks :"
chk2ce:

#  compare_res chk2cb "RMS dq generated by lmgf compare to lm" $gfdq $lmdq .003 pass
#  chk2cb:

if ($?dqtol2 == 0) set dqtol2 = .003
echo -n "$space RMS dq ($gfdq) generated by lmgf within tol $dqtol2 of lm ($lmdq)? ... "
if (`echo $lmdq $gfdq $dqtol2 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

set nbas = `grep 'nbas =' $testoutg | awk '{k=match($0,"nbas = [0-9]*") ; print substr($0,k+6,RLENGTH-6)}'`
if ($?detol2 == 0) set detol2 = `echo .0075 $nbas | awk '{print $1*$2}'`
echo -n "$space sumev generated by lmgf within tol $detol2 of lm? ... "
if (`echo $ediff 0 $detol2 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

if ($?ldau > 0) then
set dmatstol = .004
echo -n "$space RMS change in dmat ($rmsdmlmgf) generated by lmgf within tol $dmatstol of lm ($rmsdmlm)? ... "
if (`echo $rmsdmlmgf $rmsdmlm $dmatstol | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
#zcmpmfiles_res_0 chk2dm "... Max deviation in dmats.$ext from reference" $dmatstol pass 4 dmats.$ext $testdir/dmats.$ext
chk2dm:
endif

if (! $?gmaxdif) then
else if (! $?mcx) then
  echo "$space ... mcx not installed ... no check on stdout"
else
  if (! $?stdotol2) set stdotol2 = 1e-6
  echo -n "$space maximum numerical difference in stdout ($gmaxdif) <= tol ($stdotol2) ? ... "
  if (`echo ' ' | awk -v maxdif=$gmaxdif -v tol=$stdotol2 '{print (maxdif <= tol)}'` == 1) then
    echo yes
  else
    echo no "(max diff = $gmaxdif)"
    unset pass
  endif
endif

if ($?clean) then
else if ($?pass) then
    echo "$space test 2 PASSED ($ext)"
    set jobspassed = ($jobspassed 2)
else
    echo "$space test 2 FAILED ($ext)"
    set failed = ($failed 2)
endif

chk2e:

echo $joblist | grep 3 >/dev/null
if ($status) goto chk3e
# echo "$space ... density matrix checks have been disabled"
goto chk3e
cat <<EOF

         --- Test 3.  Generate the density matrix, 3rd order ASA H ---
         The site-diagonal and off-diagonal density matrix can be calculated.
         This test will compare against density-matrix calculated without the use
         of symmetry operations.
EOF

set pass
if ($?lmargs3 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk3e
endif
if ( $?have_mc == 0 && $?quiet) then
  echo "$space ... skipping test : no automatic checks unless mcx program installed"
  goto chk3e
endif
if ($?gfargs3 == 0) set gfargs3 = ($lmargs3)
query chk31 chk3e 'run this test'
chk31:
# ... Look for executables
findcmd chk31a rdcmd "$path" "$topdir"
chk31a:
findcmd chk31b lm "$path" "$topdir"
chk31b:
findcmd chk31c lmstr "$path" "$topdir"
chk31c:
findcmd chk31d lmgf "$path" "$topdir"
chk31d:

set twoc = 0
set refout=$testdir/out.lm.$ext.3c refoutg=$testdir/out.lmgf.$ext.3c testout=out.lmgf

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f $rmlst3"
             rm -f $rmlst3
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk3e
endif
echo "$space cp $cplst ."
             cp $cplst .
runjob chk32 /dev/null "lmstr $ext -vnit=0 $lmargs3"
chk32:
runjob chk33 $testout "lmgf $ext -vnit=0 -vtwoc=$twoc $gfargs3 --no-iactiv"
chk33:
echo ' '
if (! -e $testdir/sdmat.$ext.fbz) then
  echo "$space ... skipping site-diagonal test : missing file $testdir/sdmat.$ext.fbz"
  goto chk34a
endif
query chk34 chk34a "generate the site-diagonal density-matrix"
chk34:
if ($?quiet) echo "$space ... generate the site-diagonal density-matrix"
echo "$space rm -f vshft.$ext"
             rm -f vshft.$ext
echo "$space lmgf $ext -vnit=1 -vtwoc=$twoc -vldmat=1 $gfargs3 --iactiv <<EOF >$testout"
             lmgf $ext -vnit=1 -vtwoc=$twoc -vldmat=1 $gfargs3 --iactiv <<EOF >$testout




q
EOF
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully.  Site density-matrix is in file dmat.$ext"
if ($?add0) then
 echo -n "         ..." ; $add0 $testout
 echo -n "         ..." ; $add0 dmat.$ext
endif

cmp dmat.$ext $testdir/sdmat.$ext.fbz >/dev/null
call zcmpnfiles chk3c1 "8 dmat.$ext $testdir/sdmat.$ext.fbz"
chk3c1:
echo -n "$space files dmat.$ext and $testdir/sdmat.$ext.fbz equivalent to 8 digits? ..."
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 6 digits? ... "
  call zcmpnfiles chk3c2 "6 dmat.$ext $testdir/sdmat.$ext.fbz"
chk3c2:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif

chk34a:
echo ' '
if (! -e $testdir/dmat.$ext.fbz) then
  echo "$space ... skipping site-diagonal test : missing file $testdir/dmat.$ext.fbz"
  goto chk36
endif
query chk35 chk36 "generate the off-diagonal density-matrix"
chk35:
if ($?quiet) echo "$space ... generate the off-diagonal density-matrix"
echo "$space rm -f vshft.$ext"
             rm -f vshft.$ext
echo "$space lmgf $ext -vnit=1 -vtwoc=0 -vldmat=2 $gfargs3 --iactiv <<EOF >$testout"
             lmgf $ext -vnit=1 -vtwoc=0 -vldmat=2 $gfargs3 --iactiv <<EOF >$testout




q
EOF
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully.  Off-diagonal density-matrix is in file dmat.$ext"
if ($?add0) then
 echo -n "         ..." ; $add0 $testout
 echo -n "         ..." ; $add0 dmat.$ext
endif

chk3c:
call zcmpnfiles chk3c3 "8 dmat.$ext $testdir/dmat.$ext.fbz"
chk3c3:
echo -n "$space files dmat.$ext and $testdir/dmat.$ext.fbz equivalent to 8 digits? ..."
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 6 digits? ... "
  call zcmpnfiles chk3c4 "6 dmat.$ext $testdir/dmat.$ext.fbz"
chk3c4:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif

chk36:
echo ' '
if (! -e $testdir/sdmat.$ext.fbz-sph) then
  echo "$space ... skipping site-diagonal test : missing file $testdir/sdmat.$ext.fbz-sph"
  goto chk38
endif
query chk37 chk38 "generate site-diagonal density-matrix using spherical harmonics "
chk37:
if ($?quiet) echo "$space ... generate the site-diagonal density-matrix, spherical harmonics"
echo "$space rm -f vshft.$ext"
             rm -f vshft.$ext
echo "$space lmgf $ext -vnit=1 -vtwoc=0 -vldmat=1 -vsharm=1 $gfargs3 <<EOF >$testout"
             lmgf $ext -vnit=1 -vtwoc=0 -vldmat=1 -vsharm=1 $gfargs3 <<EOF >$testout




q
EOF
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully.  Off-diagonal density-matrix is in file dmat.$ext"
if ($?add0) then
 echo -n "         ..." ; $add0 $testout
 echo -n "         ..." ; $add0 dmat.$ext
endif

chk3cc:
call zcmpnfiles chk3c8 "8 dmat.$ext $testdir/sdmat.$ext.fbz-sph"
chk3c8:
echo -n "$space files dmat.$ext and $testdir/sdmat.$ext.fbz-sph equivalent to 8 digits? ..."
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 6 digits? ... "
  call zcmpnfiles chk3c9 "6 dmat.$ext $testdir/sdmat.$ext.fbz-sph"
chk3c9:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif
chk38:

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 3 PASSED ($ext)"
    set jobspassed = ($jobspassed 3)
else
    echo "$space test 3 FAILED ($ext)"
    set failed = ($failed 3)
endif
chk3e:

echo $joblist | grep 4 >/dev/null
if ($status) goto chk4e

cat <<EOF

         --- Test 4.  Iterate to self-consistency, 3rd order ASA H ---
         Tests iterations towards self-consistency.  See also test 6.

EOF
if ("$ext" == gas) then
cat <<EOF
         Checks special to the GaAs case :

         * screened output density, model dielectric function

         * freezing of potential shifts each iteration (useful for insulators)

         * ellipse+line integration scheme (using for insulators)

         * scaling of sphere radii to fill space

         * gamma-representation

EOF
endif
if ("$ext" == mnn) then
cat <<EOF
         For the mnn case, it checks downfolding (idx=2).  Also,
         A shift to the gamma-representation is necessary owing to a pole
         in the d channel of large empty sphere.

EOF
endif
if ("$ext" == mnpt) then
cat <<EOF
         For the mnpt case, test downfolds the Mn f-orbitals.  See
         $testdir/out.lmgf.mnpt.sc.3c.idx=1 for the test without downfolding.

EOF
endif
if ("$ext" == fe) then
cat <<EOF
         The Fe test checks option HAM_BXSCAL=1.
         The self-consistent magnetic moment is reduced to 2.20 from 2.26 (HAM_BXSCAL=0).

EOF
endif
if ("$ext" == eras) then
cat <<EOF
         The ErAs case tests the LDA+U Green's function.
         Spherical harmonics are used, to make the density-matrix nearly diagonal.

         In this case, convergence is difficult to obtain with
         third-order potential functions (no problems with 2nd order).

EOF
endif
set pass

set gmaxdif = 0
set twoc = 0
if ("$ext" == eras) then
set twoc = 0
endif
set refout=$testdir/out.lmgf.$ext.sc.3c testout=out.lmgf
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk4e
endif
if ($?lmargs4 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk4e
endif
if ($?gfargs4 == 0) set gfargs4 = ($lmargs4)
query chk41 chk4e 'run this test'
chk41:
# ... Look for executables
findcmd chk41a rdcmd "$path" "$topdir"
chk41a:
findcmd chk41b lm "$path" "$topdir"
chk41b:
findcmd chk41c lmstr "$path" "$topdir"
chk41c:
findcmd chk41d lmgf "$path" "$topdir"
chk41d:
if ($?MPIK) then
  set lmgf = "mpirun -n $nmpi lmgf"
endif

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f $rmlst4"
             rm -f $rmlst4
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk4e
endif
echo "$space cp $cplst ."
             cp $cplst .
runjob chk42 /dev/null "lmstr $ext -vnit=0 $lmargs4"
chk42:
runjob chk43 $testout "lmgf $ext -vnit=0 $gfargs4 --no-iactiv"
chk43:
echo "$space ... invoke lmgf for 10 iterations"

runjob chk44 $testout "$lmgf $ext -vnit=10 -vtwoc=$twoc $gfargs4 --no-iactiv"
chk44:
# if ($?quiet) goto chk4c
echo "$space Program lmgf returned successfully."
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif
echo " "
echo "$space ... Compare integrated DOS to file $refout":
cat $testout | awk '{if ($2 == "integrated") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p
echo ---
zcat  $refout | awk '{if ($2 == "integrated") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p

call showout chk46 SV
chk46:
call showout chk47 CPU
chk47:
call zdiffiles chk49 "CPU 1 $testout $refout"
chk49:

chk4c:
echo ' '
call qprint chk4ca "$space ... automatic pass checks :"
chk4ca:
set dqend = `cat $testout | grep SV: | tail -1 | awk '{print $3}' | sed 's/D/e/'`
if ($?dqtol4 == 0) set dqtol4 = 1e-4

set dqend = `cat $testout | grep SV: | tail -1 | awk '{print $3}' | sed 's/D/e/'`
compare_res_0 chk4cd "RMS DQ in moments of last iteration" $dqend $dqtol4 pass
chk4cd:

if ($?detol4 == 0) set detol4 = 1e-5
set etest = `cat $testout      | grep "  it" | awk '{print $8}' | tail -1`
set eref  = `zcat $refout | grep "  it" | awk '{print $8}' | tail -1`
compare_res chk4cf "ehf of last iteration" $etest $eref $detol4 pass
chk4cf:

if (! $?gmaxdif) then
else if (! $?mcx) then
  echo "$space ... mcx not installed ... no check on stdout"
else
  if (! $?stdotol2) set stdotol2 = 1e-6
  echo -n "$space maximum numerical difference in stdout ($gmaxdif) <= tol ($stdotol2) ? ... "
  if (`echo ' ' | awk -v maxdif=$gmaxdif -v tol=$stdotol2 '{print (maxdif <= tol)}'` == 1) then
    echo yes
  else
    echo no "(max diff = $gmaxdif)"
    unset pass
  endif
endif

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 4 PASSED ($ext)"
    set jobspassed = ($jobspassed 4)
else
    echo "$space test 4 FAILED ($ext)"
    set failed = ($failed 4)
endif
chk4e:

echo $joblist | grep 5 >/dev/null
if ($status) goto chk5e

# NB: for this test to work, you must EITHER define lmargs5 OR set variable fermi in the setup.
# In the former case, the Fermi level is found by running lm; in the latter it uses the preset value.

cat <<EOF

         --- Test 5.  Magnetic exchange coupling, 3rd order ASA H ---
         This test tests linear response calculations of magnetic exchange.
EOF
if ("$ext" == mnpt) then
cat <<EOF

         The mnpt test downfolds the f orbitals, demonstrating calculation of J in a downfolded repsn.

         This test also uses the spin-averaged gamma-representation (gamrep=2).
         Exchange parameters should not be affected except for J_00
         (gamrep=1 and gamrep=2 should produce the same J).

EOF
endif
if ("$ext" == fe) then
cat <<EOF

         Fe tests a case with one atom in the unit cell.

         It generates an RPA estimate for Tc, and resolves
         it by the radial component of q (file tcqres.fe).
         To see a picture of how contributions to MF and RPA Tc get resolved by q,
         do the following (you must have plot installed)
            cp gf/test/tcqres.fe .
            fplot -f gf/test/plot.tcqres

       * The test uses a 16x16x16 k-mesh.  It is not sufficiently fine to converge
         the exchange parameters.
         The MF and RPA Tc are (see out.lm.job11 after job executes) are:

             MF vs RPA estimate for Tc.  Tc k-resolved (file tcqres)
             mesh  2/3 1/N sum_q [J(q)-J(0)]     2/3 [1/N sum_q [J(q)-J(0)]^-1]^-1
                          meV     K                     meV     K
                        98.151  1139.0                72.941   846.5
              2x        98.151  1139.0                72.399   840.2
         With nk=24 divisions, the result is
                          meV     K                     meV     K
                        90.289  1047.8                64.788   751.9
              2x        90.289  1047.8                64.055   743.4
         With nk=32 divisions, the result is
                          meV     K                     meV     K
                        89.785  1042.0                61.357   712.1
              2x        89.785  1042.0                62.143   721.2
         With nk=48 divisions, the result is
                          meV     K                     meV     K
                        92.197  1069.9                62.152   721.3
              2x        92.197  1069.9                61.857   717.8


       * You can repeat this calculation with an (artificial) 2-atom supercell.  Do
           cp gf/test/test.gf/site.fe .
           gf/test/test.gf --fe2 fe 5
         In the two-atom supercell the results are:
           Mean-field estimate:  Tc = 10.025 mRy =  1054.8 K
           RPA estimate: Tc = 65.440 meV =  759.4 K

       * The exchange calculation is repeated with an LDA+U potential on the Fe d channels.
         A Hubbard U=0.1 Ry is used.  Occupation numbers are set to:
            0.49 for minority t2g
            0.31 for minority eg
            0.84 for majority t2g
            0.92 for majority eg
         These are the occupation numbers resulting from an LDA calculation.

         The effect of U is to push the majority d bands down and the minority d
         bands up slightly.  But the magnetic moment and fermi level change, so the
         net effect is to leave the majority d levels roughly fixed, while the minority
         levels get pushed up.  Self-consistency with the addition of U (not done here)
         modifies this effect.

       * The calculation may be repeated once more with shifts in C and Delta added
         so as to (approximately) reproduce QSGW energy bands.  Shifts were generated
         using the Levenberg-Marquardt algorithm (see input file testing/ctrl.fe)
         and are contained in file vext.fe.

EOF
endif
if ("$ext" == ni) then
cat <<EOF

         To draw a picture of the SW's in meV (you must have the FPLOT package installed), do:
            echo 0 450 5 10 | plbnds -fplot -ef=0 -scl=13.6 bnds.ni
            fplot -disp -f plot.plbnds

EOF
else if ("$ext" == fevg) then
cat <<EOF

         Tests exchange interactions for an alloy in the CPA approximation.

EOF

else if ("$ext" == fept) then
cat <<EOF

         Tests exchange interactions for FePt, with spin orbit coupling included

EOF


else if ("$ext" == femnpt) then
cat <<EOF

         Tests exchange interactions for a complex alloy in the CPA approximation.

EOF
else if ("$ext" == co) then
cat <<EOF

         Co has 2 atoms/cell; spin waves are generated for this two-atom case.

       * A similar calculation in the fcc lattice yields the following estimates for Tc (K):
               nk   Mean-Field    RPA
               12     1547.0     1200.9
               16     1435.1     1207.1
               24     1523.6     1209.3
               32     1496.7     1219.5
         Compare to the hcp case:
               12     1562.5     1291.4
               16     1506.1     1258.9
               24     1559.6     1285.7

       * This test is repeated with an LDA+U potential on the Co d channels.
         A Hubbard U=0.1 Ry is used.  Occupation numbers (0.6 for minority,
         0.9 for majority) are close to the self-consistent result.

         The effect of U is to leave the minority d bands almost unchanged, but to
         push the majority bands down -0.04Ry.  The primary result is to increase
         the magnetic moment to 1.7 and the exchange coupling by about 20%.
         Essentially the same result without the addition of U can be obtained by
         manually shifting the minority enu and C parameters by -0.04 Ry.

EOF
endif
if ("$ext" == mnn) then
cat <<EOF

         The mnn test downfolds the d orbitals on anions and empty spheres,
         and the p orbital on one empty sphere.

         A shift to the gamma-representation is necessary owing to
         a pole in the d channel of large empty sphere.

EOF
endif

set pass; set ldaupass ; unset ldaupass; set vextpass ; unset vextpass
set twoc = 0
set gmaxdif = 0
set refout=$testdir/out.$ext.3c.nk=$nkj testout=out.lm
if ($?lmargs5 == 0) then
  if ($?gfargs5 == 0) then
    echo "$space ... skipping test : case $ext has not been set up"
    goto chk5e
  endif
endif
if ($?gfargs5 == 0) then
  set gfargs5 = ($lmargs5)
endif
if (! -e $refout.job10) then
  echo "$space ... skipping test : missing file $refout.job10"
  goto chk5e
endif

if ($?lmargs5 > 0) then
  set lmargs50 = ($lmargs5)
endif
set gfargs50 = ($gfargs5)

query chk51 chk5e 'run this test'
chk51:
# ... Look for executables
findcmd chk51a rdcmd "$path" "$topdir"
chk51a:
findcmd chk51b lm "$path" "$topdir"
chk51b:
findcmd chk51c lmstr "$path" "$topdir"
chk51c:
findcmd chk51d lmgf "$path" "$topdir"
chk51d:
if ($?MPIK) then
   echo " "
   echo "$space sorry, MPIK doesn't work yet for this branch ..."
   goto chk5e
endif

# exchange interactions not parallelized ... iogfrs isn't ready
# if ($?MPIK) then
#   set lmgf = "mpirun -n $nmpi lmgf"
# endif

echo "$space ... set up ASA strux and starting potential"
# echo "$space rm -f $filel"
#              rm -f $filel
echo "$space rm -f $rmlst5" out.lm.job1  out.lm.job10  out.lm.job11
             rm -f $rmlst5  out.lm.job1  out.lm.job10  out.lm.job11
if ($?clean) then
  echo "$space rm -f $testout out.job10 out.job11 out.job1 qj0z j0z dj0dz qz"
               rm -f $testout out.job10 out.job11 out.job1 qj0z j0z dj0dz qz
  goto chk5e
endif
echo "$space cp $cplst ."
             cp $cplst .
# if ($ext == "fept") then
#   echo "$space cp rsta.fept.lmgf.ntht=0 rsta.fept"
#                cp rsta.fept.lmgf.ntht=0 rsta.fept
# endif

if ($?lmargs5 == 0) then
runjob chk52x /dev/null "lmstr $ext -vnit=0 $gfargs5"
endif
runjob chk52 /dev/null "lmstr $ext -vnit=0 $lmargs5"
chk52:
echo "$space ... run band program to find the Fermi level"
echo "$space lm $ext -vnit=0 $lmargs5 -vtwoc=$twoc --iactiv >out.lm"
             lm $ext -vnit=0 $lmargs5 -vtwoc=$twoc --iactiv <<END >out.lm
1

q
END
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lm returned successfully; output is in file out.lm"
set fermi = `grep '  ef ' log.$ext | tail -1 | awk '{print $7}'`
echo "$space fermi level generated by lm program = $fermi"
chk52x:
if ("$ext" == co || "$ext" == fe) then
  echo "$space rm -f dmats.$ext"
               rm -f dmats.$ext
endif
# If set up via lmgf
if ($?gfargs5x) then
echo "$space ... setup with lmgf"
if ("$ext" == femnpt) then
echo "$space echo ef=0.0000000  vconst=0.0742239 > vshft.$ext"
echo         echo ef=0.0000000  vconst=0.0742239 > vshft.$ext
echo "$space lmgf $ext $gfargs5x > out.lmgf"
             lmgf $ext $gfargs5x > out.lmgf
endif
echo "$space ... done setup"
endif

runjob chk53 $testout.job10 "$lmgf $ext $gfargs5 -vgfmod=10 --no-iactiv -ef=$fermi+.000"
chk53:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout.job10
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout.job10
endif

set gfargx = "$gfargs5"
if (! $?ldaupass && ! $?vextpass && -r $testdir/tcqres.$ext) then
  set gfargx = "$gfargs5 --tcqres"
endif

runjob chk54 $testout.job11 "lmgf $ext $gfargx -vgfmod=11 --no-iactiv -ef=$fermi"
chk54:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout.job11
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout.job11
endif

if (-e $refout.job24) then
  runjob chk55 $testout.job24 "lmgf $ext $gfargs5 -vgfmod=24 --no-iactiv -ef=$fermi"
  chk55:
else
  echo " "
  echo "$space ... skipping longitudinal susceptibility  : missing file $refout.job24"
endif

echo ' '
echo "$space ... Compare sum rule to lines in file $refout.job10":
cat $testout.job10 | grep 'Sum rule'
echo '---'
cat $refout.job10  | grep 'Sum rule'

echo ' '
echo "$space ... Compare charge, mag. moment to lines in file $refout.job10":
cat $testout.job10 | awk '{if ($3 == "amom") {print; getline; print}}'
echo '---'
cat $refout.job10  | awk '{if ($3 == "amom") {print; getline; print}}'

echo ' '
echo "$space ... J_ij to lines in file $refout.job11":
cat $testout.job11 | awk '{if ($3 == "J_ij") {print; getline; print}}'
echo '---'
cat $refout.job11  | awk '{if ($3 == "J_ij") {print; getline; print}}'

echo ' '
call zdiffilesx chk56 "CPU 1 $testout.job10 $refout.job10"
chk56:

call zdiffiles chk57 "CPU 1 $testout.job11 $refout.job11"
chk57:

if (-e $refout.job24) then
call zdiffiles chk58 "CPU 1 $testout.job24 $refout.job24"
chk58:
endif

chk5c:
echo ' '
call qprint chk5ca "$space ... automatic pass checks :"
chk5ca:

# check file tcqres.$ext
switch ("$gfargx")
case "*--tcqres*":
call zcmpnfiles chk57a "6 tcqres.$ext $testdir/tcqres.$ext"
chk57a:
echo -n "$space ... files tcqres.$ext $testdir/tcqres.$ext equivalent to 6 digits? ... "
if ($retval == 0) then
  echo yes
else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
  echo ok "($retval difference(s) of $ncharfile)"
else
  echo no "($retval difference(s) of $ncharfile)"
  unset pass
endif
endsw

set refjr = jr.$ext.$nkj
if ($?ldaupass) set refjr = jr.$ext.$nkj.ldau
if ($?vextpass) set refjr = jr.$ext.$nkj.vext
if ($?add0) then
  echo -n "         ..." ; $add0 jr.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer jr.$ext
endif
call zcmpnfiles chk5cb "8 jr.$ext $testdir/$refjr"
chk5cb:
echo -n "$space ... files jr.$ext $testdir/$refjr equivalent to 8 digits? ... "
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 6 digits? ... "
  call zcmpnfiles chk5cc "6 jr.$ext $testdir/$refjr"
chk5cc:
 if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif

if (-e $refout.job24) then
set refchiq = chiq.$ext.$nkj
call zcmpnfiles chk5cd "8 chiq.$ext $testdir/$refchiq"
chk5cd:
echo -n "$space ... files chiq.$ext $testdir/$refchiq equivalent to 8 digits? ... "
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 6 digits? ... "
  call zcmpnfiles chk5ce "6 chiq.$ext $testdir/$refchiq"
chk5ce:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif

set refchiq = chiq2.$ext.$nkj
call zcmpnfiles chk5cf "8 chiq2.$ext $testdir/$refchiq"
chk5cf:
echo -n "$space ... files chiq2.$ext $testdir/$refchiq equivalent to 8 digits? ... "
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 6 digits? ... "
  call zcmpnfiles chk5cg "6 chiq2.$ext $testdir/$refchiq"
chk5cg:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif
endif

set bndref = bnds.$ext
if ($?ldaupass) set bndref = bnds.$ext.ldau
if ($?vextpass) set bndref = bnds.$ext.vext

if (-e syml.$ext && -e bnds.$ext && -e $testdir/$bndref) then

if ($?add0) then
  echo -n "         ..." ; $add0 bnds.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer bnds.$ext
endif
call zcmpnfiles chk5ch "4 bnds.$ext $testdir/$bndref"
chk5ch:
echo -n "$space ... spin wave files bnds.$ext and $testdir/$bndref equivalent to 4 digits? ..."
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 3 digits? ... "
  call zcmpnfiles chk5ch2 "3 bnds.$ext $testdir/$bndref"
chk5ch2:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (1000*ndiff/ntot<1.)}'` == 1) then
      echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif
endif  # If there are bands

if (! -e $refout.job10.zdependence) then
  echo "$space ... skipping test : missing file $refout.job10.zdependence"
  goto chk5e0
endif
if (! -e $refout.job11.zdependence) then
  echo "$space ... skipping test : missing file $refout.job11.zdependence"
  goto chk5e0
endif
if ($?ldaupass | $?vextpass) goto chk5e0

echo "$space"
echo "$space"
echo "$space The following check illustrates the generation of the density of states and"
echo "$space energy derivative of exchange J0, dJ0/dE, in a uniform mesh along the real axis."
echo "$space"
echo "$space The DOS is generated invoking lmgf in the standard mode (GF MODE=1)",
echo "$space but using a special energy contour (EMESH contour type 2)."
echo "$space This contour is documented in file doc/gf.html"
echo "$space"
echo "$space dJ0/dz is generated invoking lmgf in the exchange mode (GF MODE=10)",
echo "$space and also using the special energy contour (EMESH contour type 2)."
echo "$space"
echo "$space After completion, shell script $testdir/getJq0z is invoked to extract"
echo "$space q(z) and dJ0/dz from the output files.  It can also numerically integrate"
echo "$space dJ0/dz to obtain J0(z) as a function of z."

if (! $?have_mc) then
echo "$space (The integration is omitted now because a recent version of the"
echo "$space matrix calculator program 'mcx' needs be installed in your path)"
endif
echo "$space"

query chk5z1 chk5e0 'run test illustrating energy dependence'
chk5z1:

echo "$space ... use EMESH for energy integration along real axis"
echo "$space     (the next line copies input file, stripping out EMESH mode 10)"
echo "$space cat $testdir/ctrl.$ext | awk '"'{if (match($1,"EMESH=.*") && $2 == "10") next; else print}'"' >ctrl.$ext"
             cat $testdir/ctrl.$ext | awk '{if (match($1,"EMESH=.*") && $2 == "10") next; else print}' >ctrl.$ext

echo "$space ... invoke lmgf, mode 1, to illustrate trapezoidal integration of N(E)"
runjob chk5z2 $testout.job1 "lmgf $ext $gfargs5 -vgfmod=1 --no-iactiv -vef=$fermi"
chk5z2:
echo ' '
if ($?add0) then
  echo -n "         ..." ; $add0 $testout.job1
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout.job1
endif
call zdiffilesx chk5z3 "CPU 1 $testout.job1 $refout.job1.zdependence 1e-4"
chk5z3:

echo "$space ... invoke lmgf, mode 10 and 11 to generate dJ(z)/dz"
runjob chk5z4a $testout.job10 "lmgf $ext $gfargs5 -vgfmod=10 --no-iactiv -vef=$fermi"
chk5z4a:
runjob chk5z4b $testout.job11 "lmgf $ext $gfargs5 -vgfmod=11 --no-iactiv -vef=$fermi"
chk5z4b:

echo ' '
if ($?add0) then
  echo -n "         ..." ; $add0 $testout.job10
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout.job10
endif
call zdiffiles chk5z5 "CPU 1 $testout.job10 $refout.job10.zdependence"
chk5z5:

if ($?add0) then
  echo -n "         ..." ; $add0 $testout.job11
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout.job11
endif
call zdiffiles chk5z6 "CPU 1 $testout.job11 $refout.job11.zdependence"
chk5z6:

echo "$space ... invoke script $testdir/getJq0z to extract the energy-dependence"
echo "$space     of the integrated charge qz from $testout.job1 (saved in file qz)"
echo "$space     and dJ0/dz from $testout.job11 (saved in file dj0dz)."
echo "$space"
echo "$space     File $testdir/j0z.$nkj.$ext-elliptical-summary contains a summary"
echo "$space     of the generation of qz and j0 by integration on elliptical"
echo "$space     contour at successive endpoints.  It provides a useful comparison"
echo "$space     to the files qz and djdz (or j0z) generated by $testdir/getJq0z."
echo "$space     $testdir/getJq0z -facJ0=1000 -ef0=$fermi $testout.job11 $testout.job1"
echo -n "$space "
                 $testdir/getJq0z -facJ0=1000 -ef0=$fermi $testout.job11 $testout.job1

if ($?have_fplot && $?have_ghostscript && $?slow && -e $testdir/j0z.$nkj.$ext-elliptical-summary && ! $?quiet) then
#  query chk5z7 chk5zc "plot generated qz with elliptical integration in file $testdir/j0z.$nkj.$ext-elliptical-summary"
#  chk5z7:
#    echo "$space     $plot -pr10 -disp qz -lt 2 $testdir/j0z.$nkj.$ext-elliptical-summary"
#                     $plot -pr10 -disp qz -lt 2 $testdir/j0z.$nkj.$ext-elliptical-summary
query chk5z8 chk5zc "plot j0z of first atom with elliptical integration in file $testdir/j0z.$nkj.$ext-elliptical-summary"
chk5z8:
  echo "$space     $plot -pr10 -disp j0z -lt 2 -ord x3/.001 -s circ $testdir/j0z.$nkj.$ext-elliptical-summary"
                   $plot -pr10 -disp j0z -lt 2 -ord x3/.001 -s circ $testdir/j0z.$nkj.$ext-elliptical-summary
endif

chk5zc:
echo ' '
call qprint chk5zca "$space ... automatic pass checks :"
chk5zca:
call zcmpnfiles chk5zcb "8 jr.$ext $testdir/jr.$ext.$nkj.zdependence"
chk5zcb:
echo -n "$space ... files jr.$ext $testdir/jr.$ext.$nkj.zdependence equivalent to 8 digits? ... "
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 6 digits? ... "
  call zcmpnfiles chk5zcc "6 jr.$ext $testdir/jr.$ext.$nkj.zdependence"
chk5zcc:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
    echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif

chk5e0:
if ($?ldauargs5 && ! $?ldaupass && $?pass) then
  set ldaupass
  set lmargs5 = ($ldauargs5 $lmargs50)
  set gfargs5 = ($ldauargs5 $gfargs50)
  echo " "
  echo "$space Repeat test $ext with LDA+U parameters set."
  set refout=$testdir/out.$ext.3c.nk=$nkj.ldau
  query chk51 chk5e1 'run this test'
endif
chk5e1:

if ($?vextargs5 && ! $?vextpass && $?pass) then
  set vextpass
  set lmargs5 = ($vextargs5 $lmargs50)
  set gfargs5 = ($vextargs5 $gfargs50)
  echo " "
  echo "$space Repeat test $ext with potential shifts read from file vext."
  set refout=$testdir/out.$ext.3c.nk=$nkj.vext
  query chk51 chk5e 'run this test'
endif

if (! $?gmaxdif) then
else if (! $?mcx) then
  echo "$space ... mcx not installed ... no check on stdout"
else
  if (! $?stdotol5) set stdotol5 = 1e-6
  echo -n "$space maximum numerical difference in stdout ($gmaxdif) <= tol ($stdotol5) ? ... "
  if (`echo ' ' | awk -v maxdif=$gmaxdif -v tol=$stdotol5 '{print (maxdif <= tol)}'` == 1) then
    echo yes
  else
    echo no "(max diff = $gmaxdif)"
    unset pass
  endif
endif

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 5 PASSED ($ext)"
    set jobspassed = ($jobspassed 5)
else
    echo "$space test 5 FAILED ($ext)"
    set failed = ($failed 5)
endif

chk5e:


echo $joblist | grep 6 >/dev/null
if ($status) goto chk6e

cat <<EOF

         --- Test 6.  Linear response ---
         Test linear response branch of code.
EOF
if ("$ext" == nife) then
cat <<EOF
         The nife case tests the noncollinear linear response case;
         the ground state is noncollinear.  (Warning: this test takes a while!)

         It is interesting to observe the evolution of energy with iteration.
         You can extract the magnetization and energy from the save file using
           cat save.nife | awk '{print NR,\$8,\$9}' | sed s/mmom=// | sed s/ehf=//

         After the initial spike in the first iterations (from large deviations
         in self-consistency in the charge), the energy drops smoothly as
         the magnetic angles relax.  The magnetization also falls smoothly.
         This test only finds a metastable state.
EOF
endif
if ("$ext" == mnn) then
cat <<EOF
         The Mn magnetic moment is very sensitive to environment;
         you can see this by the large eigenvalue in the response matrix:
            zcat $testdir/psta.mnn >out ; mcx out -herm -evl
         This system is difficult to converge to self-consistency owing to
         the 'finicky' magnetic degrees of freedom.  Because of this
         sensitivity, this test updates the response matrix frequently.
         (Test 4 succeeds in this case by using Anderson mixing.)

         Also a shift to the gamma-representation is necessary owing to
         a pole in the d channel of large empty sphere.

EOF
endif

set pass
set refout=$testdir/out.lmgf.$ext.sc-scr.3c testout=out.lmgf
if ($?lmargs6 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk6e
 else if (! -e $refout) then
   echo "$space ... skipping test : missing reference file $refout"
   goto chk6e
endif
if ($?gfargs6 == 0) set gfargs6 = ($lmargs6)

# lm -vnit=0 -vscr=11 nife --no-iactiv
# find initial potential shift for Fermi level
# lmgf -vnit=1 -vscr=0 nife --iactiv
# blank
#

set pass
query chk61 chk6e 'run this test'
chk61:

# ... Look for executables
findcmd chk61a rdcmd "$path" "$topdir"
chk61a:
findcmd chk61b lm "$path" "$topdir"
chk61b:
findcmd chk61c lmstr "$path" "$topdir"
chk61c:
findcmd chk61d lmgf "$path" "$topdir"
chk61d:
if ($?MPIK) then
   echo " "
   echo "$space sorry, MPIK doesn't work yet for this branch ..."
   goto chk6e
endif

set twoc = 0

echo "$space ... set up ASA strux and starting potential"
if ($?touchl) then
 touch $touchl
endif

# echo "$space rm -f $filel"
#              rm -f $filel

echo "$space rm -f $rmlst6"
             rm -f $rmlst6
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk6e
endif
echo "$space cp $cplst ."
             cp $cplst .

runjob chk62 /dev/null "lmstr $ext -vnit=0 $lmargs6"
chk62:
echo "$space ... set up initial potential and get approximate Fermi level"
echo "$space $lmgf $ext -vnit=0 -vscr=0 $gfargs6 --iactiv <<EOF >out.lmgf"
             $lmgf $ext -vnit=0 -vscr=0 $gfargs6 --iactiv <<EOF >out.lmgf
1

sf


q
EOF
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully; output is in file out.lmgf"


echo "$space ... set up initial potential, intra-atomic U, and initial response function"
echo "$space lmgf $ext -vnit=0 -vscr=11 $gfargs6 --iactiv <<EOF >>out.lmgf"
             lmgf $ext -vnit=0 -vscr=11 $gfargs6 --iactiv <<EOF >>out.lmgf
1

q
EOF
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully; output is in file out.lmgf"

echo "$space ... saving file psta.$ext to psta.$ext.init"
cp psta.$ext psta.$ext.init

echo "$space ... invoke lmgf for self-consistency cycle"
runjob chk65 '>>'$testout "$lmgf $ext $gfargs6 --no-iactiv"
chk65:

call showout chk66 SV
chk66:
call showout chk67 CPU
chk67:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif
call zdiffiles chk69 "CPU 3 $testout $refout"
chk69:


chk6c:
if ($?dqtol6 == 0) set dqtol6 = 5e-5
if ($?dptol6 == 0) set dptol6 = 1e-8
if ($?dmtol6 == 0) set dmtol6 = 1e-4
if ($?detol6 == 0) set detol6 = 1e-5

echo ' '
call qprint chk6ca "$space ... automatic pass checks :"
chk6ca:

call zcmpnfiles chk6cb "5 psta.$ext.init $testdir/psta.$ext"
chk6cb:
echo -n "$space ... files psta.$ext.init and $testdir/psta.$ext compare to 6 digits? ... "
if ($retval == 0) then
  echo yes
else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
  echo ok "($retval difference(s) of $ncharfile)"
else
  echo no "($retval difference(s) of $ncharfile)"
  unset pass
endif

set dqend = `cat $testout | grep SV: | tail -1 | awk '{print $3}' | sed 's/D/e/'`
compare_res_0 chk6cd "RMS DQ in moments of last iteration" $dqend $dqtol6 pass
chk6cd:

set mmend = `cat $testout | grep mmom=  | tail -1 | awk '{match($0,"mmom=[^ ]*"); print substr($0,RSTART+5,RLENGTH-5)}'`
set mmref  = `zcat $refout | grep mmom=  | tail -1 | awk '{match($0,"mmom=[^ ]*"); print substr($0,RSTART+5,RLENGTH-5)}'`
compare_res chk6ce "magnetic moment of last iteration" $mmend $mmref $dmtol6 pass
chk6ce:

set etest = `cat $testout      | grep "  it" | awk '{print $6}' | tail -1`
set eref  = `zcat $refout | grep "  it" | awk '{print $6}' | tail -1`
compare_res chk6cf "ehf of last iteration" $etest $eref $detol6 pass
chk6cf:

if (! $?gmaxdif) then
else if (! $?mcx) then
  echo "$space ... mcx not installed ... no check on stdout"
else
  if (! $?stdotol6) set stdotol6 = 1e-6
  echo -n "$space maximum numerical difference in stdout ($gmaxdif) <= tol ($stdotol6) ? ... "
  if (`echo ' ' | awk -v maxdif=$gmaxdif -v tol=$stdotol6 '{print (maxdif <= tol)}'` == 1) then
    echo yes
  else
    echo no "(max diff = $gmaxdif)"
    unset pass
  endif
endif

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 6 PASSED ($ext)"
    set jobspassed = ($jobspassed 6)
else
    echo "$space test 6 FAILED ($ext)"
    set failed = ($failed 6)
endif

chk6e:

echo $joblist | grep 7 >/dev/null
if ($status) goto chk7e

cat <<EOF

         --- Test 7.  Noncollinear special checks ---
EOF
if ("$ext" == fccfe) then
cat <<EOF

         The fccfe test case checks some noncollinear branches.

         If 2nd order potential functions are used, lmgf should
         produce the same results as lm for converged k and energy meshes.
         (To verify this, set twoc=1 in test 7 for fccfe, script $testfile)

         The symmetry operations

EOF
set twoc=0
endif

if ("$ext" == bi2te3) then
cat <<EOF

         The bi2te3 test contains four stacked layers of Bi2Te3 separated by a vacuum layer.
         It demonstrates a Bi2Te3 surface, and compares SO coupling in lm and lmgf.

EOF
if (! $?quiet) then
cat <<EOF
         To draw a picture of the lmgf-generated spectral function,
         invoke the following after this test completes:
           SpectralFunction.sh -e -f 10 -w -2:2

EOF
endif # quiet
set bndref = bnds.$ext
endif  # bi2te3

if ("$ext" == fept2) then
cat <<EOF

         The fept2 test computes G with SO coupling included for:
         (1) a reference frame rotated about y by 45 degrees, with symetry operation
         (2) a reference frame rotated about y by 45 degrees, without symop
         (3) Euler angles applied to both Fe and Pt atoms, same rotation, without symop

EOF
if (! $?quiet) then
cat <<EOF
         This check establishes that (1),(2), and (3) yield the same the physical quantities
         (total energy density, spin- and orbital- moments, torque)

         There is only one symop comparing tests (1) and (2).
         As a further check on correct implemention of symops, do the following and
         compare files out and out2 generated for the unrotated case.
         This establishes that symops correctly generate the diagonal part of G in an SO coupled case.

           lmstr        ctrl.fept2 -vnit=-1 -vnk=6 -vso=1 -vrot=0 --quit=rho
           cp vshft-bk.fept2 vshft.fept2
           rm -f mixm.fept eula.fept
           lmgf --nosym ctrl.fept2 -vnit=-1 -vnk=6 -vso=1 -vrot=0 --quit=rho > out
           cp vshft-bk.fept2 vshft.fept2
           rm -f mixm.fept eula.fept
           lmgf         ctrl.fept2 -vnit=-1 -vnk=6 -vso=1 -vrot=0 --quit=rho > out2
           diff out out2

EOF
endif # quiet
endif  # fept2

if ("$ext" == fept) then
cat <<EOF

         Case FePt: Show noncollinear code reproduces collinear case by comparing:
         1. FePt (Euler angles all zero)
         2. Rotate quantization axis for Pt to -z and simultaneously flip Pt magnetic moment.

         This test also demonstrates a noncollinear calculation with downfolding :
         In the main test, the Fe f and Pt f are downfolded.
         In a more stringent test, the Fe p,f and Pt f are downfolded.

         It also checks a fully relativistic, noncollinear calculation.

EOF
if (! $?quiet) then
cat <<EOF
         Note: the lm and lmgf densities were made self-consistent separately.
         The RMS DQ should be very small for both  of lm and both runs of lmgf.
         This is the key test that checks whether the code is working properly.
         Note that lm and lmgf Fermi levels are a little different.
         Also downfolding is handled differently in the lm and lmgf cases;
         the band sum energies (and total energies) are not very close.

         For more extensive tests, try
           gf/test/test.ncgfso
           gf/test/test.ncgfmst

EOF
endif # quiet
set twoc=0
endif # fept

if ("$ext" == co) then
cat <<EOF

         This check of hcp Co tests symmetrization of both orbital and spinor parts.
         The Euler angles have been set to be consistent with a 2-fold rotation around z.

EOF
if (! $?quiet) then
cat <<EOF
         The result should be the same whether or not the symmetry operation is used.
         To confirm this, make the density self-consistent after this test completes:
           lmgf co -vconvc=1e-6 -vldmat=2 -vtwoc=0 -ef=-0.021546 -vnit=20 -vtwoc=0 -vnc=1 -vtwoc=0 -ef=-0.025725 -vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vcbya=1.625 --iactiv=no
         Then make one another pass, this time suppressing symmetry operations:
           lmgf --nosym co -vconvc=1e-6 -vldmat=2 -vtwoc=0 -ef=-0.021546 -vnit=1 -vtwoc=0 -vnc=1 -vtwoc=0 -ef=-0.025725 -vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vcbya=1.625 --iactiv=no
         You should see that the RMS DQ is very small.

         You can also confirm that the symmetry is kept because the spinors are also rotated.
         To see this, edit ctrl.co, remove the comment ('#') in this line
           SYMGRP  r2z(0,0,-0.8125) # NOSPIN=t
         This tells lmgf to rotate the orbitals but not the spinors when symmetrizing.
         Do one further iteration
           lmgf co -vconvc=1e-6 -vldmat=2 -vtwoc=0 -ef=-0.021546 -vnit=1 -vtwoc=0 -vnc=1 -vtwoc=0 -ef=-0.025725 -vnk=12 -vbzj=0 -vccor=f -vgfmod=1 -vcbya=1.625 --iactiv=no
         The output density is now far from self-consistent.


EOF
endif
set twoc=0
endif

if ("$ext" == eras) then
cat <<EOF

         The eras test checks the noncollinear LDA+U branch.
         Spherical harmonics are used, to make the density-matrix nearly diagonal.
         Also, the diagonal-only approximation is used.

         Compare some results of lmgf (sumev, Fermi level, output q) to program lm.
         lmgf and lm 2C approximation produce nearly identical results for a fine k-mesh.
         (To verify this, set twoc=1 in test 7 for eras, script $testfile)

EOF
set twoc=0
endif
set pass
if (! $?stdotol7) set stdotol7 = 3e-6
#goto xxx
if ($?lmargs7 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk7e
endif
if ($?gfargs7 == 0) set gfargs7 = ($lmargs7)
query chk71 chk7e 'run this test'
chk71:
# ... Look for executables
findcmd chk71a rdcmd "$path" "$topdir"
chk71a:
findcmd chk71b lm "$path" "$topdir"
chk71b:
findcmd chk71c lmstr "$path" "$topdir"
chk71c:
findcmd chk71d lmgf "$path" "$topdir"
chk71d:
if ($?MPIK) then
  set lmgf = "mpirun -n $nmpi lmgf"
endif

if ($twoc == 0) then
set refoutlm=$testdir/out.lm.$ext.nc1.3c refoutgf=$testdir/out.lmgf.$ext.nc1.3c
else
set refoutlm=$testdir/out.lm.$ext.nc1.2c refoutgf=$testdir/out.lmgf.$ext.nc1.2c
endif
set testoutlm=out.lm testoutgf=out.lmgf
if ($ext == "bi2te3") then
  set refoutlm = "$refoutgf"
  set testoutlm = "$testoutgf"
endif

if (! $?quiet) then
echo ' '
echo "$space output file for lm = $testoutlm"
echo "$space output file for GF = $testoutgf"
echo "$space reference          = $refoutgf"
if ("$ext" == fept) then
echo "$space downfolding test   = out.idxp2.lmgf"
echo "$space fully relativistic = out.rel2.$ext"
endif
echo
endif

if ($?haveout && "$ext" == fept2) goto chk7haveout
if ($?haveout) goto chk7ca

if (! -e $refoutlm) then
  echo "$space ... skipping test : missing reference file $refoutlm"
  goto chk7e
endif
if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chk7e
endif
echo "$space rm -f $rmlst7"
             rm -f $rmlst7
if ($?clean) then
  echo "$space rm -f $testoutlm $testoutgf out.rel2.lmgf"
               rm -f $testoutlm $testoutgf out.rel2.lmgf
  goto chk7e
endif

if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chk7e
endif
if ("$ext" != bi2te3 && ! $?quiet) then
echo "$space NB: density-matrix is hermitian, but does not come out exactly so"
echo "$space owing to integration errors in the complex plane."
echo " "
endif

if ($twoc == 0) then
echo "$space ... set up ASA strux and starting potential for 3c hamiltonian"
else
echo "$space ... set up ASA strux and starting potential for 2c hamiltonian"
endif

touch ctrl.$ext
echo "$space rm -f *.$ext"
             rm -f *.$ext
if ("$ext" == fccfe) then
echo "$space cp $testdir/ctrl.$ext ."
             cp $testdir/ctrl.$ext .
echo "$space cp $testdir/eula.ran eula.$ext"
             cp $testdir/eula.ran eula.$ext
else if ("$ext" == co) then
echo "$space cp $testdir/{{ctrl,eula}.co,q.bulk} ."
             cp $testdir/{{ctrl,eula}.co,q.bulk} .
else if ("$ext" == fept) then
echo "$space cp $testdir/{ctrl}.fept ."
             cp $testdir/{ctrl}.fept .
echo "$space cp $cplst ."
             cp $cplst .
else if ("$ext" == fept2) then
echo "$space cp $cplst ."
             cp $cplst .
else if ("$ext" == bi2te3) then
echo "$space cp $cplst ."
             cp $cplst .
else
echo "$space cp $cplst ."
             cp $cplst .
echo "$space cp $testdir/eula.ran eula.$ext"
             cp $testdir/eula.ran eula.$ext
endif

if ("$ext" == fept) then
# goto XXXX
runrdcmd chk7ca %11f $testoutlm "-cat:BNDTEST --noerr ctrl.$ext"
endif # fept

if ("$ext" == fept2) then
if ($?MPIK) then
runrdcmd chk7b1 %11f $testoutgf "-cat:GFMTEST --noerr ctrl.$ext"
else
runrdcmd chk7b1 %11f $testoutgf "-cat:GFTEST --noerr ctrl.$ext"
endif
chk7b1:
call zdiffiles chk7haveout "CPU 1 $testoutgf $refoutgf"
endif # fept2

if ($ext == "bi2te3") then
if ($?MPIK) then
runrdcmd chk7ca %11f $testoutlm "-cat:GFMTEST --noerr ctrl.$ext"
endif
runrdcmd chk7ca %11f $testoutlm "-cat:GFTEST --noerr ctrl.$ext"
endif #bi2te3

runjob chk73c1 /dev/null "lmstr $lmargs7 $ext"
chk73c1:

# Set up potential parameters
runjob chk73c2 $testoutgf "lmgf $ext -vnit=0 -vtwoc=$twoc $gfargs7 --no-iactiv"
chk73c2:

echo "$space ... Run lm band pass to find Fermi level and band sum :"
echo "$space lm $ext -vnit=1 -vtwoc=$twoc $lmargs7 >$testoutlm"
             lm $ext -vnit=1 -vtwoc=$twoc $lmargs7 <<END >$testoutlm





q
END
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lm returned successfully; output is in file $testoutlm"
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutlm
endif

#    set fermi = `grep Fermi $testoutlm | tail -1 | awk '{print $4}' | sed 's/;//'`
#    set sumev = `grep 'Sum occ' $testoutlm | tail -1 | awk '{print $4}' | sed 's/,//'`

chk7ca:

set testout=$testoutlm refout=$refoutlm
# Fermi level is 4th argument in lines containing "Fermi energy"; return last occurence
get_resf chk7cb fermi fermiref "Fermi energy" 4 0 ";"
chk7cb:
get_resf chk7cc sumev sumevref "bands:" 4 0  ";"
chk7cc:

if ($?ldau > 0) then
get_resf chk7cd rmsdmlm rmsdmlmref "RMS diff in dens mat" 5 0 "mat("
chk7cd:
set rmsdmlm = `echo "$rmsdmlm" | sed 's/)//'`
set rmsdmlmref = `echo "$rmsdmlmref" | sed 's/)//'`
endif

set lmdq  = `grep 'RMS DQ=' log.$ext | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`
if ("$ext" == fept || "$ext" == bi2te3) then
set lmdq  = `grep 'RMS DQ=' log.$ext | awk '{print $9}'  | sed 's/DQ=//'`
endif

if (! $?quiet) then
  echo ' '
  echo "$space Fermi energy from lm               = $fermi"
  echo "$space Fermi energy of reference          = $fermiref"
  set diff = `echo $fermi $fermiref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space band sum from lm                   = $sumev"
  echo "$space band sum of reference              = $sumevref"
  set diff = `echo $sumev $sumevref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space RMS DQ (output-input) from lm      = $lmdq[1]"
  if ($#lmdq > 1) then
  echo "$space ditto, from second invocation      = $lmdq[2]"
  endif

  if ($?ldau > 0) then
  echo ' '
  echo "$space RMS change in dmat from lm         = $rmsdmlm"
  echo "$space RMS change of reference            = $rmsdmlmref"
  endif

  echo ' '
else
  echo "$space Found Fermi level = $fermi, sumev = $sumev,  rms DQ between input and output moments = $lmdq"
endif

if ($?haveout) goto chk7haveout

if ($?add0) then
  echo -n "         ..." ; $add0 $testoutlm
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutlm
endif

call zdiffiles chk73c4 "CPU 1 $testoutlm $refoutlm $stdotol7"
chk73c4:

if ("$ext" == bi2te3) then
  goto chk7ga
endif # bi2te3

if ("$ext" == fept) then
# XXXX:
rm -f log.$ext
runrdcmd chk73c5 %11f $testoutgf "-cat:GFTEST --noerr ctrl.$ext"
chk73c5:
echo
echo "         ... redo last calculation folding down the Fe p orbitals"
runrdcmd chk73c6 %11f out.idxp2.lmgf "-cat:DNFTEST --noerr ctrl.$ext"
chk73c6:
call zdiffiles chk73c7 "CPU 1 out.idxp2.lmgf $refoutgf:h/out.idxp2.$ext - $mcexcl~excl=gippad~excl=gfzerq~excl=spin"
chk73c7:

echo
echo "         ... redo calculation fully relativistically (rotate spinors only!)"
runrdcmd chk73c8 %11f out.rel2.lmgf "-cat:FRTEST --noerr ctrl.$ext"
chk73c8:
call zdiffiles chk73c9 "CPU 1 out.rel2.lmgf $refoutgf:h/out.rel2.$ext"
chk73c9:
goto chk7ga

endif # fept

echo ' '
echo "$space ... invoke lmgf using Fermi level output from lm:"
echo "$space rm -f vshft.$ext mq.$ext"
             rm -f vshft.$ext mq.$ext
if ($?cpagain) then
echo "$space cp $cplst ."
             cp $cplst .
endif
echo "$space $lmgf $ext -vldmat=2 -vtwoc=$twoc -ef=$fermi -vnit=1 -vtwoc=$twoc $gfargs7 <<EOF >out.lmgf"
             $lmgf $ext -vldmat=2 -vtwoc=$twoc -ef=$fermi -vnit=1 -vtwoc=$twoc $gfargs7 <<EOF >out.lmgf






q
EOF
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully.  The following data were output:"
chk7ga:
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutgf
 echo -n "         ..." ; $add0 dmat.$ext
endif

chk7haveout:

set testout=$testoutgf refout=$refoutgf
set zval = `grep 'zval' log.$ext | tail -1 | awk '{print $2}'`
set vshft = `grep 'gfasa  vshft:' log.$ext | tail -1 | awk '{print $3}'`
if ($ext == "bi2te3") then
get_resf chk7gc vshft vshftref "IOVSHF" 6 0 vconst=
chk7gc:
endif
set dc = `echo null | awk -v zval=$zval -v vshft=$vshft '{print zval*vshft}'`  # zval*vshft

if ($ext == "bi2te3") then
get_resf chk7gb sumevg sumevgref "total" 4 5 ";"
endif
get_resf chk7gb sumevg sumevgref "sumev=" 4 0 sumev=
chk7gb:

if ("$ext" == fept2) then

  set amoms   = (`extract-lines 'START LMGF' Exit 1 $testout  | vextract . mmom | tail -1 `)
  set amomc   = (`extract-lines 'START LMGF' Exit 2 $testout  | vextract . mmom | tail -1 `)
  set amoma   = (`extract-lines 'START LMGF' Exit 3 $testout  | vextract . mmom | tail -1 `)
  set amomcr  = (`extract-lines 'START LMGF' Exit 2 $refoutgf | vextract . mmom | tail -1 `)
  set amomar  = (`extract-lines 'START LMGF' Exit 3 $refoutgf | vextract . mmom | tail -1 `)

  set ehfs   = (`extract-lines 'START LMGF' Exit 1 $testout  | vextract h ehf `)
  set ehfc   = (`extract-lines 'START LMGF' Exit 2 $testout  | vextract h ehf `)
  set ehfa   = (`extract-lines 'START LMGF' Exit 3 $testout  | vextract h ehf `)
  set ehfcr  = (`extract-lines 'START LMGF' Exit 2 $refoutgf | vextract h ehf `)
  set ehfar  = (`extract-lines 'START LMGF' Exit 3 $refoutgf | vextract h ehf `)

  set dqs    = (`extract-lines 'START LMGF' Exit 1 $testoutgf | grep 'RMS DQ' | vextract . DQ`)
  set dqc    = (`extract-lines 'START LMGF' Exit 2 $testoutgf | grep 'RMS DQ' | vextract . DQ`)
  set dqa    = (`extract-lines 'START LMGF' Exit 3 $testoutgf | grep 'RMS DQ' | vextract . DQ`)

  if ($?have_mc) then
    set torqs = (`extract-lines --quiet 'MAGCPA:' IORBTM: 1 $testoutgf | grep Torq | tail -1 | awk '{print $2, $3, $4}' | mcx -f3f12.8 . -w:nohead .`)
    set torqc = (`extract-lines --quiet 'MAGCPA:' IORBTM: 2 $testoutgf | grep Torq | tail -1 | awk '{print $2, $3, $4}' | mcx -f3f12.8 . -w:nohead .`)
    set torqa = (`extract-lines --quiet 'MAGCPA:' IORBTM: 3 $testoutgf | grep Torq | tail -1 | awk '{print $2, $3, $4}' | mcx -f3f12.8 . -w:nohead .`)
  else
    set torqs = (`extract-lines --quiet 'MAGCPA:' IORBTM: 1 $testoutgf | grep Torq | tail -1 | awk '{print $2, $3, $4}'`)
    set torqc = (`extract-lines --quiet 'MAGCPA:' IORBTM: 2 $testoutgf | grep Torq | tail -1 | awk '{print $2, $3, $4}'`)
    set torqa = (`extract-lines --quiet 'MAGCPA:' IORBTM: 3 $testoutgf | grep Torq | tail -1 | awk '{print $2, $3, $4}'`)
  endif

  set doms = (`extract-lines --quiet 'IORBTM:' LMGF: 1 $testoutgf | grep L+ | tail -1 | awk '{print $NF}'`)
  set domc = (`extract-lines --quiet 'IORBTM:' LMGF: 2 $testoutgf | grep L+ | tail -1 | awk '{print $NF}'`)
  set doma = (`extract-lines --quiet 'IORBTM:' LMGF: 3 $testoutgf | grep L+ | tail -1 | awk '{print $NF}'`)

  echo " "
  echo    "$space Compare rotated coordinate system against rotated spins"
  echo -n "$space magnetic moment,  rotated coordinates   = $amomc"
  if (`echo $amomc $amomcr | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  unset pass
  echo " ... test failed (reference value = $amomcr)"
  endif

  echo -n "$space magnetic moment,  rotated coord+symop   = $amoms"
  if (`echo $amoms $amomc | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with nosymops"
  else
  unset pass
  echo " ... test failed (nosymop value = $amomac)"
  endif

  echo -n "$space magnetic moment,  rotated spins         = $amoma"
  if (`echo $amoma $amomar | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  unset pass
  echo " ... test failed (reference value = $amomar)"
  endif

  echo " "
  echo -n "$space HF energy,  rotated coordinates         = $ehfc"
  if (`echo $ehfc $ehfcr | awk -v tol=$detol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  unset pass
  echo " ... test failed (reference value = $ehfcr)"
  endif
  echo -n "$space HF energy,  rotated coord+symops        = $ehfs"
  if (`echo $ehfa $ehfs | awk -v tol=$detol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with nosymops"
  else
  unset pass
  echo " ... test failed (nosymop value = $ehfaa)"
  endif
  echo -n "$space HF energy,  rotated spins               = $ehfa"
  if (`echo $ehfa $ehfar | awk -v tol=$detol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  unset pass
  echo " ... test failed (reference value = $ehfar)"
  endif

  echo " "
  echo "$space y component of Pt torque, rotated coord = $torqc[2]"
  echo "$space y component of Pt torque, with symops   = $torqs[2]"
  echo "$space y component of Pt torque, rotated spins = $torqa[2]"
  echo -n "$space torques agree within tol ($dttol7)? ... "
  if (`echo $torqc[2] $torqa[2] | awk -v tol=$dttol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
    echo " yes"
  else
    unset pass
    echo " no"
  endif

  echo " "
  echo "$space Pt orbital moment, rotated coord        = $domc"
  echo "$space Pt orbital moment, rotated coord+symop  = $doms"
  echo "$space Pt orbital moment, rotated spins        = $doma"
  echo -n "$space agree within tol ($dqtol7)? ... "
  if (`echo $domc $doma | awk -v tol=$dqtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
    echo " ... yes"
  else
    unset pass
    echo " ... no"
  endif

  echo
  echo -n "$space RMS DQ ($dqc), rotated coordinates less than tol ($dqtol7)? ... "
  if (`echo $dqc 0 $dqtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
    echo yes
  else
    echo no
    unset pass
  endif

  echo -n "$space RMS DQ ($dqs), rotated coordinates+symops less than tol ($dqtol7)? ... "
  if (`echo $dqs 0 $dqtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
    echo yes
  else
    echo no
    unset pass
  endif

  echo -n "$space RMS DQ ($dqa), rotated spins less than tol ($dqtol7)? ... "
  if (`echo $dqa 0 $dqtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
    echo yes
  else
    echo no
    unset pass
  endif

  goto chk7e0
endif # fept2

set ediff  = `echo $sumevg $sumev  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`                    # lmgf-lm difference in sumev
set ediff  = `echo $sumevg $sumev $dc | awk '{{k=($1-($2+$3))>0?($1-($2+$3)):(($2+$3)-$1);} print k}'`  # ditto, with d.c. included
set gfdq   = `grep 'RMS DQ=' log.$ext | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`
set sumevgdc = `echo $sumevg $dc | awk '{print ($1-$2)}'`

if ("$ext" == fept) then
  set efgf = `grep efermi $testout | awk '{print $7}' | tail -1`                               # FePt lm and lmgf use different Fermi level
  set dc = `echo null | awk -v zval=$zval -v ef=$fermi -v efgf=$efgf '{print zval*(efgf-ef)}'` # change in sumev from delta Ef: lm vs lmgf
  set ediff  = `echo $sumevg $sumev $dc | awk '{{k=($1-($2+$3))>0?($1-($2+$3)):(($2+$3)-$1);} print k}'`
  set sumevgdc = `echo $sumevg $dc | awk '{print ($1-$2)}'`
  set vshft = `grep 'gfasa  vshft:' log.$ext | awk '{print $3}'`                               # 6 entries: sr(FM,AFM) sr+dnf(3x) fr(FM)
  set gfdq   = `grep 'RMS DQ=' $testout | awk '{print $9}'  | sed 's/DQ=//'`                   # pick up 2 numbers from sr, 1st and 2nd invocation

  set amomc   = (`extract-lines 'START LMGF' Exit 2 $testout | vextract . mmom `)
  set amoma   = (`extract-lines 'START LMGF' Exit 4 $testout | vextract . mmom `)
  set amomf   = (`extract-lines 'START LMGF' Exit 2 out.rel2.lmgf | vextract h mmom `)
  set amomc   = (`extract-lines 'START LMGF' Exit 2 $testout | vextract . mmom `)
  set amomcr  = (`extract-lines 'START LMGF' Exit 2 $refoutgf  | vextract . mmom `)
  set amomar  = (`extract-lines  'START LMGF' Exit 4 $refoutgf  | vextract . mmom `)
  set amomfr  = (`extract-lines 'START LMGF' Exit 2 $refoutgf:h/out.rel2.$ext | vextract h mmom `)
  set amomp   = (`grep mmom= out.idxp2.lmgf | tail -1 | vextract . mmom `)
  set amompr  = (`grep mmom= $topdir/gf/test/out.idxp2.fept | tail -1 | vextract . mmom `)

  set ehfc   = (`extract-lines 'START LMGF' Exit 2 $testout | vextract x ehf `)
  set ehfa   = (`extract-lines 'START LMGF' Exit 4 $testout | vextract x ehf `)
  set ehff   = (`extract-lines 'START LMGF' Exit 2 out.rel2.lmgf | vextract h ehf `)
  set ehfc   = (`extract-lines 'START LMGF' Exit 2 $testout | vextract x ehf `)
  set ehfcr  = (`extract-lines 'START LMGF' Exit 2 $refoutgf  | vextract x ehf `)
  set ehfar  = (`extract-lines  'START LMGF' Exit 4 $refoutgf  | vextract x ehf `)
  set ehffr  = (`extract-lines 'START LMGF' Exit 2 $refoutgf:h/out.rel2.$ext | vextract h ehf `)
  set ehfp   = (`grep ehf= out.idxp2.lmgf | tail -1 | vextract x ehf `)
  set ehfpr  = (`grep ehf= $topdir/gf/test/out.idxp2.fept | tail -1 | vextract x ehf `)

  echo " "
  echo    "$space Compare lmgf with Pt collinear, anticollinear"
  echo -n "$space SR magnetic moment,  collinear     = $amomc"
  if (`echo $amomc $amomcr | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $amomcr"
  endif

  echo -n "$space SR magnetic moment,  anticollinear = $amoma"
  if (`echo $amoma $amomar | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $amomar"
  endif

  echo -n "$space FR magnetic moment,  anticollinear = $amomf"
  if (`echo $amomf $amomfr | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $amomfr"
  endif


  echo -n "$space SR anticollinear, p downfolded     = $amomp"
  if (`echo $amomp $amompr | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $amompr"
  endif

# set diff = `echo $amomc $amoma | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  set diff = `echo $amomc $amoma | awk '{{k=($1-$2)} print k}'`
  echo    "$space SR(coll) SR(anticoll) difference   = $diff"

  set diff = `echo $amoma $amomf | awk '{{k=($1-$2)} print k}'`
  echo    "$space SR-FR difference (both anticoll)   = $diff"


  echo " "
  echo -n "$space SR HF energy,  collinear           = $ehfc"
  if (`echo $ehfc $ehfcr | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $ehfcr"
  endif

  echo -n "$space SR HF energy,  anticollinear       = $ehfa"
  if (`echo $ehfa $ehfar | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $ehfar"
  endif

  echo -n "$space FR HF energy,  anticollinear       = $ehff"
  if (`echo $ehff $ehffr | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $ehffr"
  endif


  echo -n "$space SR anticollinear, p downfolded     = $ehfp"
  if (`echo $ehfp $ehfpr | awk -v tol=$dmtol7 '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
  echo " ... agrees with reference"
  else
  echo " ... reference value = $ehfpr"
  endif

# set diff = `echo $ehfc $ehfa | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  set diff = `echo $ehfc $ehfa | awk '{{k=($1-$2)} print k}'`
  echo    "$space SR(coll) SR(anticoll) difference   = $diff"

  set diff = `echo $ehfa $ehff | awk '{{k=($1-$2)} print k}'`
  echo    "$space SR-FR difference (both anticoll)   = $diff"

  echo ' '
endif  # FePt

if (! $?quiet) then
  echo "$space Pade estimate for Fermi level shift = $vshft[1]"
  if ($#vshft > 1) then
  echo "$space ditto, from second invocation       = $vshft[2]"
  endif
  echo "$space lm band structure energy            = $sumev"
  echo "$space gf band structure energy            = $sumevg"
  echo "$space d.c. term zval*vshft                = $dc"
  echo "$space difference - d.c.                   = $ediff"
endif # quiet

if ($ext == "bi2te3") then
  goto chk73c8e
endif

set testout=$testoutgf refout=$refoutgf
if ($?ldau > 0) then
get_resf chk7ce rmsdmlmgf rmsdmlmgfref "RMS diff in dens mat" 5 0 "mat("
chk7ce:
set rmsdmlmgf = `echo "$rmsdmlmgf" | sed 's/)//'`
set rmsdmlmgfref = `echo "$rmsdmlmgfref" | sed 's/)//'`
echo "$space RMS change in dmat from lmgf       = $rmsdmlmgf"
echo "$space RMS change in dmat from lm         = $rmsdmlm"
#  echo "$space RMS change of reference            = $rmsdmlmgfref"
endif  # ldau

echo "$space lm deviation in output-input charge = $lmdq[1]"
if ($#lmdq > 1) then
echo "$space ditto, from second invocation       = $lmdq[2]"
endif
echo "$space gf deviation in output-input charge = $gfdq[1]"
if ($#lmdq > 1) then
echo "$space ditto, from second invocation       = $gfdq[2]"
endif

#  xxx:
#  set testoutlm=out.lm testoutgf=out.lmgf
if ($?quiet) then
else
  echo ' '
  echo "$space ... Magnetic components output by lmgf:"
  cat $testoutgf | awk '{if ($1 == "AMAGNC:" && $5 == "density") {print;getline;print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,7p
  echo ' '
  echo "$space ... compared to output of lm:"
  cat $testoutlm | awk '{if ($1 == "AMAGNC:" && $5 == "density") {print;getline;print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,7p

  echo ' '
  echo "$space ... Magnetic forces output by lmgf:"
  cat $testoutgf | awk '{if ($1 == "local") {print;getline;print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p
  echo "$space ... compared to output of lm:"
  cat $testoutlm | awk '{if ($1 == "local") {print;getline;print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p

  call showout chk73c8c CPU
chk73c8c:
  echo ' '
  call zdiffiles chk73c8d "CPU 1 $testoutgf $refoutgf"
chk73c8d:

endif  # quiet

chk73c8e:

echo ' '
call qprint chk7aca "$space ... automatic pass checks :"
chk7aca:
if ($?detol7 == 0) set detol7 = .0003
if ($?dqtol7 == 0) set dqtol7 = .001
echo -n "$space RMS dq ($gfdq[1]) generated by lmgf within tol $dqtol7 of lm ($lmdq[1])? ... "
if (`echo $lmdq[1] $gfdq[1] $dqtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

if ($#lmdq > 1 && $#gfdq > 1) then
echo -n "$space RMS dq ($gfdq[2]) generated by lmgf within tol $dqtol7 of lm ($lmdq[2])? ... "
if (`echo $lmdq[2] $gfdq[2] $dqtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
endif

echo -n "$space sumev ($sumevgdc) generated by lmgf within tol ($detol7) of lm? ... "
if (`echo $ediff 0 $detol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

if ("$ext" == fccfe) then
set dphigf = `cat $testoutgf | awk '{if ($1 == "local") {getline;print $2}}'`
set dthegf = `cat $testoutgf | awk '{if ($1 == "local") {getline;print $3}}'`
set dphilm = `cat $testoutlm | awk '{if ($1 == "local") {getline;print $2}}'`
set dthelm = `cat $testoutlm | awk '{if ($1 == "local") {getline;print $3}}'`
set ddphi = `echo $dphigf $dphilm | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ddthe = `echo $dthegf $dthelm | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

if (`echo $ddphi $ddthe | awk '{print ($1<.03 && $2<.0005)}'`) then
  echo "$space magnetic force 1st atom generated by lmgf within .0005 of  lm? ..." yes
else
  echo "$space magnetic force 1st atom generated by lmgf within .0005 of  lm? ..." no
  unset pass
endif
endif # fccfe

# if (-e syml.$ext && -e bnds.$ext && -e $testdir/$bndref) then
if ("$ext" == bi2te3) then
if ($?add0) then
  echo -n "         ..." ; $add0 bnds.$ext
else if ($?poszer) then
  echo -n "         ..." ; $poszer bnds.$ext
endif
call zcmpnfiles chk7ch "4 bnds.$ext $testdir/$bndref"
chk7ch:
echo -n "$space lm-generated band file bnds.$ext and $testdir/$bndref equivalent to 4 digits? ..."
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 3 digits? ... "
  call zcmpnfiles chk7ch2 "3 bnds.$ext $testdir/$bndref"
chk7ch2:
  if ($retval == 0) then
    echo yes
  else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (1000*ndiff/ntot<1.)}'` == 1) then
      echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif
endif  # If there are bands

if ("$ext" == bi2te3) then

#   echo
#   echo "$space ... checking spectral function files against reference"
  set fn = (spf.$ext spf1.$ext)
  set nl = 500
  set ndig = 5
  if ! ($?dostol) set dostol = 1e-3
  chk7ch3:
    if ($?add0) then
      echo -n "         ..." ; $add0 $fn[1]
    else if ($?poszer) then
      echo -n "         ..." ; $poszer $fn[1]
    endif
    zcmpmfiles_res_tol chk7ch4 "Max deviation in ($nl lines) $fn[1] from reference $testdir/$fn[1]" $dostol pass $ndig $fn[1] $testdir/$fn[1] $nl
    chk7ch4:
    shift fn
  if ($#fn > 0) goto chk7ch3

endif

if ("$ext" == fept) then
  echo -n "$space magnetic moment ($amomp) with Fe p downfolded within tol ($dqtol7) of reference? ... "
  if (`echo $amomp $amompr $dqtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
    echo yes
  else
    echo no
    unset pass
  endif

  set gfdq[1] = `grep 'RMS DQ=' out.idxp2.lmgf | awk '{print $9}'  | sed 's/DQ=//'`
  set gfdqr   = `grep 'RMS DQ=' $topdir/gf/test/out.idxp2.fept | awk '{print $9}'  | sed 's/DQ=//'`
  echo -n "$space RMS DQ ($gfdqr) with Fe p downfolded within tol ($dqtol7) of reference? ... "
  if (`echo $gfdq[1] $gfdqr $dqtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
    echo yes
  else
    echo no
    unset pass
  endif

  echo -n "$space magnetic moment ($amomf) for fully relativistic calculation within tol ($dmtol7) of reference? ... "
  if (`echo $amomf $amomfr $dmtol7 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
    echo yes
  else
    echo no
    unset pass
  endif

  set gfdqr = `grep 'deviation from charge neutrality:' out.rel2.lmgf | awk '{print $NF}' | head -1`
  echo -n "$space fully relativistic calculation deviation from charge neutrality ($gfdqr) within tol (5e-5) of reference? ... "
  if (`echo 0 $gfdqr 5e-5 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
    echo yes
  else
    echo no
    unset pass
  endif

endif  # FePt

# call zcmpnfiles chk7c6 "8 dmat.$ext $testdir/dmat.$ext.fbz"
# chk7c6:
# echo -n "$space files dmat.$ext and $testdir/dmat.$ext.fbz equivalent to 8 digits? ..."
# if ($retval == 0) then
#   echo yes
# else
#   echo -n "no ... to 6 digits? ... "
#   call zcmpnfiles chk7c7 "6 dmat.$ext $testdir/dmat.$ext.fbz"
# chk7c7:
#   if ($retval == 0) then
#     echo yes
#   else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (100*ndiff/ntot<1.)}'` == 1) then
#     echo ok "($retval difference(s) of $ncharfile)"
#   else
#     echo no "($retval difference(s) of $ncharfile)"
#     unset pass
#   endif
# endif

endif

chk7e0:

if (! $?gmaxdif) then
else if (! $?mcx) then
  echo "$space ... mcx not installed ... no check on stdout"
else
  echo -n "$space maximum numerical difference in stdout ($gmaxdif) <= tol ($stdotol7) ? ... "
  if (`echo ' ' | awk -v maxdif=$gmaxdif -v tol=$stdotol7 '{print (maxdif <= tol)}'` == 1) then
    echo yes
  else
    echo no "(max diff = $gmaxdif)"
    unset pass
  endif
endif

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 7 PASSED ($ext)"
    set jobspassed = ($jobspassed 7)
else
    echo "$space test 7 FAILED ($ext)"
    set failed = ($failed 7)
endif

chk7e:

echo $joblist | grep 8 >/dev/null
if ($status) goto chk8e
if (! -e $topdir/sx) goto chk8e

cat <<EOF

         --- Test 8.  Test generation of screened-exchange self-energy ---

EOF
if ("$ext" == cdte) then
cat <<EOF
         Compute screened exchange self-energy sigma for CdTe, 3rd order potential functions

EOF
endif
set pass

set twoc = 0
set refout=$testdir/out.lmgf.$ext.sx.3c refoutlm=$testdir/out.lm.$ext.sx.3c
set testout=out.lmgf testoutlm=out.lm
if ($?lmargs8 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk8e
endif
if ($?gfargs8 == 0) set gfargs8 = ($lmargs8)
query chk81 chk8e 'run this test'
chk81:
# ... Look for executables
findcmd chk81a rdcmd "$path" "$topdir"
chk81a:
findcmd chk81b lm "$path" "$topdir"
chk81b:
findcmd chk81c lmstr "$path" "$topdir"
chk81c:
findcmd chk81d lmgf "$path" "$topdir"
chk81d:

echo "$space ... set up ASA strux and starting potential"
# echo "$space rm -f $filel"
#              rm -f $filel
echo "$space rm -f $rmlst8"
             rm -f $rmlst8
if ($?clean) then
  echo "$space rm -f $testout $testoutlm"
               rm -f $testout $testoutlm
  if (-e qasa.$ext) then
  echo "$space rm qasa.$ext qasa"
               rm qasa.$ext qasa
  endif
  goto chk8e
endif
echo "$space cp $cplst ."
             cp $cplst .
if (-e qasa.$ext) then
echo "$space cp qasa.$ext qasa"
             cp qasa.$ext qasa
endif

runjob chk82 /dev/null "lmstr $ext -vnit=0 $lmargs8"
chk82:
runjob chk83 $testout "lmgf $ext -vnit=0 -vsx=11 $gfargs8 --no-iactiv --quit=band"
chk83:
echo "$space ... invoke lmgf to make screened-exchange self-energy"
runjob chk84 $testout "lmgf $ext -vnit=1 -vsx=11 $gfargs8 --no-iactiv --quit=band"
chk84:
echo "$space ... screened-exchange self-energy second pass"
runjob chk84a '>>'$testout "lmgf $ext -vnit=1 -vsx=11 $gfargs8 --no-iactiv --quit=band"
chk84a:
echo "$space Program lmgf returned successfully."
if ($?add0) then
 echo -n "         ..." ; $add0 $testout
endif
echo "$space ... invoke lm to generate screened-exchange energy bands"
runjob chk85 $testoutlm "lm $ext -vnit=1 -vsig=1 $gfargs8 --no-iactiv --quit=band --band:fn=syml"
chk85:
echo "$space Program lm returned successfully."
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutlm
endif
if ($?quiet) goto chk8c

echo " "
echo "$space ... Compare integrated DOS to file $refout":
cat $testout | awk '{if ($2 == "integrated") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p
echo ---
zcat  $refout | awk '{if ($2 == "integrated") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p

echo ' '
echo "$space ... Compare evals at Gamma to file $refoutlm":
cat      $testoutlm | awk -v ncnt=0 '{if ($1 == "SECMAT:" && $7 == 0 && $8 == 0 && $9 == 0) {ncnt += 1; if (ncnt==1) {print;getline;print;getline;print;getline;print}}}'
echo ---
zcat $refoutlm | awk -v ncnt=0 '{if ($1 == "SECMAT:" && $7 == 0 && $8 == 0 && $9 == 0) {ncnt += 1; if (ncnt==1) {print;getline;print;getline;print;getline;print}}}'

call zdiffiles chk89 "CPU 1 $testout $refout"
chk89:
call zdiffiles chk8a "CPU 1 $testoutlm $refoutlm"
chk8a:

chk8c:
echo ' '
set bndref = bnds.$ext
call qprint chk8ca "$space ... automatic pass checks :"
chk8ca:
call zcmpnfiles chk8c1 "4 bnds.$ext $testdir/$bndref"
chk8c1:
echo -n "$space files bnds.$ext and $testdir/$bndref equivalent to 4 digits? ..."
if ($retval == 0) then
  echo yes
else
  echo -n "no ... to 3 digits? ... "
  call zcmpnfiles chk8c2 "3 bnds.$ext $testdir/$bndref"
chk8c2:
  if ($retval == 0) then
    echo yes
    else if (`echo ' ' | awk -v ndiff=$retval -v ntot=$ncharfile '{print (1000*ndiff/ntot<1.)}'` == 1) then
      echo ok "($retval difference(s) of $ncharfile)"
  else
    echo no "($retval difference(s) of $ncharfile)"
    unset pass
  endif
endif


#  set etest = `cat       $testout|grep sumev=  | grep ehk= | tail -1 | awk '{print $2}' | sed s/ehf=//`
#  set eref  = `zcat $refout |grep sumev=  | grep ehk= | tail -1 | awk '{print $2}' | sed s/ehf=//`
#  compare_res chk8cf "ehf of last iteration" $etest $eref 1e-5 pass
#  chk8cf:

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 8 PASSED ($ext)"
    set jobspassed = ($jobspassed 8)
else
    echo "$space test 8 FAILED ($ext)"
    set failed = ($failed 8)
endif
chk8e:

echo $joblist | grep 9 >/dev/null
if ($status) goto chk9e

cat <<EOF

         --- Test 9.  Disordered Local Moments  ---
EOF
if ("$ext" == fepd) then
cat <<EOF

         The fepd test checks the disordered local moments form of CPA.

EOF
set twoc=0
endif
set pass
# set testoutgf=out.lmgf ; goto xxx
if ($?lmargs9 == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chk9e
endif
if ($?gfargs9 == 0) set gfargs9 = ($lmargs9)
query chk91 chk9e 'run this test'
chk91:
# ... Look for executables
findcmd chk91a rdcmd "$path" "$topdir"
chk91a:
findcmd chk91b lm "$path" "$topdir"
chk91b:
findcmd chk91c lmstr "$path" "$topdir"
chk91c:
findcmd chk91d lmgf "$path" "$topdir"
chk91d:


set refoutgf=$testdir/out.lmgf.$ext
set testoutgf=out.lmgf
if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chk9e
endif
echo "$space rm -f $rmlst9"
             rm -f $rmlst9
if ($?clean) then
  echo "$space rm -f $testoutgf"
               rm -f $testoutgf
  goto chk9e
endif

if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chk9e
endif

echo "$space cp $cplst ."
             cp $cplst .

runjob chk93c1 /dev/null "lmstr $ext"
chk93c1:

set lmgf = lmgf
if ($?MPIK) then
  set lmgf = "mpirun -n $nmpi lmgf"
endif

echo "$space $lmgf $ext $gfargs9 >$testoutgf"
             $lmgf $ext $gfargs9 >$testoutgf
set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully."
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutgf
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutgf
endif

# xxx:


set dlmq = (`grep Qv= out.lmgf | tail -9 | vextract . Qv `)
set etotA  = ` grep  'ehf= ' $testoutgf  | awk '{print $6}' `
set entr  = ` grep  'Entropy=' $testoutgf | sed s/Entropy=//`
# set e0    = ` grep -A 3 'Legendre coefficients' $testoutgf | sed -n 1,4p | tail -1 | awk '{print $3}'`
# set e1    = ` grep -A 3 'Legendre coefficients' $testoutgf | sed -n 1,4p | tail -1 | awk '{print $4}'`
# set e2    = ` grep -A 3 'Legendre coefficients' $testoutgf | sed -n 1,4p | tail -1 | awk '{print $5}'`
# set tK    = ` grep -A 3 'Legendre coefficients' $testoutgf | sed -n 1,4p | tail -1 | awk '{print $6}'`
set sumevg = `grep 'sumev=' $testoutgf | grep -v Omega | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set mmom   = `grep 'mmom=' $testoutgf | grep -v Omega | tail -1 | awk '{match($0,"mmom=[^ ]*"); print substr($0,RSTART+5,RLENGTH-5)}'`
set ehf    = `grep 'ehf=' $testoutgf | grep -v Omega | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set gfdq   = `grep 'RMS DQ=' $testoutgf | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`

set dlmqr  = (`zcat $testdir/out.lmgf.fepd | grep Qv= | tail -9 | vextract . Qv `)
set etotAr = `zcat $testdir/out.lmgf.fepd | grep -A 3 'Legendre coefficients' | sed -n 1,4p | tail -1 | awk '{print $2}'`
set etotAr  = `zcat $testdir/out.lmgf.fepd | grep 'ehf= ' | awk '{print $6}' `
set entrr  = `zcat $testdir/out.lmgf.fepd | grep  'Entropy=' | sed s/Entropy=//`
set e0r    = `zcat $testdir/out.lmgf.fepd | grep -A 3 'Legendre coefficients' | sed -n 1,4p | tail -1 | awk '{print $3}'`
set e1r    = `zcat $testdir/out.lmgf.fepd | grep -A 3 'Legendre coefficients' | sed -n 1,4p | tail -1 | awk '{print $4}'`
set e2r    = `zcat $testdir/out.lmgf.fepd | grep -A 3 'Legendre coefficients' | sed -n 1,4p | tail -1 | awk '{print $5}'`
set tKr    = `zcat $testdir/out.lmgf.fepd | grep -A 3 'Legendre coefficients' | sed -n 1,4p | tail -1 | awk '{print $6}'`
set sumevgr = `zcat $testdir/out.lmgf.fepd | grep 'sumev=' | grep -v Omega | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set mmomr   = `zcat $testdir/out.lmgf.fepd | grep 'mmom=' | grep -v Omega | tail -1 | awk '{match($0,"mmom=[^ ]*"); print substr($0,RSTART+5,RLENGTH-5)}'`
set ehfr    = `zcat $testdir/out.lmgf.fepd | grep 'ehf=' | grep -v Omega | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set gfdqr   = `zcat $testdir/out.lmgf.fepd | grep 'RMS DQ=' | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`


echo "$space gf deviation in output-input charge = $gfdq"

if (! $?quiet) then
  echo ' '
  echo "$space DLM average charge                 = $dlmq"
  echo "$space                    reference       = $dlmqr"
  set diff = `echo $dlmq[1] $dlmqr[1] | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space DLM entropy                        = $entr"
  echo "$space                    reference       = $entrr"
  set diff = `echo $entr $entrr | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space Magnetic moment for system         = $mmom"
  echo "$space                    reference       = $mmomr"
  set diff = `echo $mmom $mmomr | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space Total energy for system            = $etotA"
  echo "$space                    reference       = $etotAr"
  set diff = `echo $etotA $etotAr | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

echo ' '
call zdiffiles chk96 "CPU 1 $testoutgf $refoutgf"
chk96:

endif

echo ' '
call qprint chk9aca "$space ... automatic pass checks :"
chk9aca:
if ($?detol9 == 0) set detol9 = 5e-6
if ($?dqtol9 == 0) set dqtol9 = 1e-6

echo -n "$space RMS dq ($gfdq) generated by lmgf within tol $dqtol9 of reference ($gfdqr)? ... "
if (`echo $gfdq $gfdqr $dqtol9 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space output moment ($mmom) generated by lmgf within tol $dqtol9 of reference ($mmomr)? ... "
if (`echo $mmom $mmomr $dqtol9 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space sumev ($sumevg) generated by lmgf within tol $detol9 of reference ($sumevgr)? ... "
if (`echo $sumevg $sumevgr $detol9 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space entropy ($entr) generated by lmgf within tol $detol9 of reference ($entrr)? ... "
if (`echo $entr $entrr $detol9 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space etot ($etotA) generated by lmgf within tol $detol9 of reference ($etotAr)? ... "
if (`echo $etotA $etotAr $detol9 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 9 PASSED ($ext)"
    set jobspassed = ($jobspassed 9)
else
    echo "$space test 9 FAILED ($ext)"
    set failed = ($failed 9)
endif

chk9e:

echo $joblist | grep 10 >/dev/null
if ($status) goto chkae

cat <<EOF

         --- Test 10.  CPA: chemical CPA ---
EOF
if ("$ext" == fev) then
cat <<EOF

         The fev test checks the chemical CPA form of CPA.

EOF
else if ("$ext" == fevg) then
cat <<EOF
         The fevg test does this test in the spin-averaged gamma representation.
         Compare outputs of the fev and fevg tests to see representation-dependence.

EOF
else if ("$ext" == fe2b) then
cat <<EOF
         The fe2b test demonstrates a CPA calculation with perturbative SO coupling.

EOF
endif
set twoc=0
set pass

set testoutgf=out.lmgf
if ($?lmargsa == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chkae
endif
if ($?haveout) goto chka4
if ($?gfargsa == 0) set gfargsa = ($lmargsa)
query chka1 chkae 'run this test'
chka1:
# ... Look for executables
findcmd chka1a rdcmd "$path" "$topdir"
chka1a:
findcmd chka1b lm "$path" "$topdir"
chka1b:
findcmd chka1c lmstr "$path" "$topdir"
chka1c:
findcmd chka1d lmgf "$path" "$topdir"
chka1d:

set lmgf = lmgf
if ($?MPIK) then
  set lmgf = "mpirun -n $nmpi lmgf"
endif

set refoutgf=$testdir/out.lmgf.$ext
if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chkae
endif

if ("$rmlsta" == "all") then
  touch ctrl.$ext
  set rmlsta = `ls *.$ext`
endif
echo "$space rm -f $rmlsta"
             rm -f $rmlsta
if ($?clean) then
  echo "$space rm -f $testoutgf"
               rm -f $testoutgf
  goto chkae
endif

if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chkae
endif

if ($cplst[1]:e == "tar") then
  echo "$space cat $cplst | tar xf -"
         echo `cat $cplst | tar tf -`
               cat $cplst | tar xf -
else
  echo "$space cp $cplst ."
               cp $cplst .
endif

runjob chka3c1 /dev/null "lmstr $gfargsa $ext"
chka3c1:

runjob chka3c2 $testoutgf "$lmgf $ext $gfargsa"
chka3c2:
if ($retval == -2) goto chkae

if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully."
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutgf
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutgf
endif

chka4:

set dlmq = (`grep Qv= out.lmgf | tail -2 | vextract . Qv `)
set etotA  = ` grep  'ehf= ' $testoutgf  | awk '{print $6}' `
set entr  = ` grep  'Entropy=' $testoutgf | sed s/Entropy=//`
set sumevg = `grep 'sumev=' $testoutgf | grep LMGF | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set mmom   = `grep 'mmom=' $testoutgf | grep -v Omega | tail -1 | awk '{match($0,"mmom=[^ ]*"); print substr($0,RSTART+5,RLENGTH-5)}'`
set ehf    = `grep 'ehf=' $testoutgf | grep -v Omega | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set gfdq   = `grep 'RMS DQ=' $testoutgf | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`

set dlmqr  = (`zcat $testdir/out.lmgf.$ext | grep Qv= | tail -2 | vextract . Qv `)
set etotAr = `zcat $testdir/out.lmgf.$ext | grep -A 3 'Legendre coefficients' | sed -n 1,4p | tail -1 | awk '{print $2}'`
set etotAr  = `zcat $testdir/out.lmgf.$ext | grep 'ehf= ' | awk '{print $6}' `
set entrr  = `zcat $testdir/out.lmgf.$ext | grep  'Entropy=' | sed s/Entropy=//`
set sumevgr = `zcat $testdir/out.lmgf.$ext | grep 'sumev=' | grep LMGF | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set mmomr   = `zcat $testdir/out.lmgf.$ext | grep 'mmom=' | grep -v Omega | tail -1 | awk '{match($0,"mmom=[^ ]*"); print substr($0,RSTART+5,RLENGTH-5)}'`
set ehfr    = `zcat $testdir/out.lmgf.$ext | grep 'ehf=' | grep -v Omega | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set gfdqr   = `zcat $testdir/out.lmgf.$ext | grep 'RMS DQ=' | tail -1 | awk '{print $9}'  | sed 's/DQ=//'`


echo "$space gf deviation in output-input charge = $gfdq"

if (! $?quiet) then
  echo ' '
  echo "$space DLM average charge                 = $dlmq"
  echo "$space                    reference       = $dlmqr"
  set diff = `echo $dlmq[1] $dlmqr[1] | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space DLM entropy                        = $entr"
  echo "$space                    reference       = $entrr"
  set diff = `echo $entr $entrr | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space Magnetic moment for system         = $mmom"
  echo "$space                    reference       = $mmomr"
  set diff = `echo $mmom $mmomr | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

  echo ' '
  echo "$space Total energy for system            = $etotA"
  echo "$space                    reference       = $etotAr"
  set diff = `echo $etotA $etotAr | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference      = $diff"

echo ' '
call zdiffiles chka6 "CPU 1 $testoutgf $refoutgf"
chka6:

endif

echo ' '
call qprint chkaaca "$space ... automatic pass checks :"
chkaaca:
if ($?detola == 0) set detola = 5e-6
if ($?dqtola == 0) set dqtola = 1e-6

echo -n "$space RMS dq ($gfdq) generated by lmgf within tol $dqtola of reference ($gfdqr)? ... "
if (`echo $gfdq $gfdqr $dqtola | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space output moment ($mmom) generated by lmgf within tol $dqtola of reference ($mmomr)? ... "
if (`echo $mmom $mmomr $dqtola | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space sumev ($sumevg) generated by lmgf within tol $detola of reference ($sumevgr)? ... "
if (`echo $sumevg $sumevgr $detola | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space etot ($etotA) generated by lmgf within tol $detola of reference ($etotAr)? ... "
if (`echo $etotA $etotAr $detola | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 10 PASSED ($ext)"
    set jobspassed = ($jobspassed 10)
else
    echo "$space test 10 FAILED ($ext)"
    set failed = ($failed 10)
endif

chkae:

echo $joblist | grep 11 >/dev/null
if ($status) goto chkbe

cat <<EOF

         --- Test 11.  Spectral functions ---
EOF
if ("$ext" == fe2b) then
cat <<EOF
         The fe2b test demonstrates calculation of the k-resolved spectral function in the CPA.

EOF
if (! $?quiet) then
cat <<EOF
         The Green's function is calculated on a uniform mesh of points on the real axis.
         The input files are constructed so that the Omega potential is self-consistent
         on this mesh.  No intermediate iterations should be required.

         To draw a picture of the spectral function, let this calculation complete,
         then do:
           SpectralFunction.sh -e -f 10 -w -5:3
         This should create {$ext}_UP.eps and {$ext}_DN.eps .

EOF
endif
set twoc=0
endif
set pass

set testoutgf=out.lmgf refoutgf=$testdir/out.specfun.lmgf.$ext
if ($?codecheck) then
  set refoutgf = $testdir/out.specfun.lmgf.codecheck.$ext
endif
if ($?lmargsb == 0) then
  echo "$space ... skipping test : case $ext has not been set up"
  goto chkbe
endif
if ($?gfargsb == 0) set gfargsb = ($lmargsb)
query chkb1 chkbe 'run this test'
chkb1:
if ($?haveout) goto chkb4
# ... Look for executables
findcmd chkb1a rdcmd "$path" "$topdir"
chkb1a:
findcmd chkb1b lm "$path" "$topdir"
chkb1b:
findcmd chkb1c lmstr "$path" "$topdir"
chkb1c:
findcmd chkb1d lmgf "$path" "$topdir"
chkb1d:

set lmgf = lmgf
if ($?MPIK) then
  set lmgf = "mpirun -n $nmpi lmgf"
endif

if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chkbe
endif
if ($rmlstb[1] == "all") then
  touch ctrl.$ext
  set rmlstb = `ls *.$ext`
endif
echo "$space rm -f $rmlstb"
             rm -f $rmlstb
if ($?clean) then
  echo "$space rm -f $testoutgf"
               rm -f $testoutgf
  goto chkbe
endif

if (! -e $refoutgf) then
  echo "$space ... skipping test : missing reference file $refoutgf"
  goto chkbe
endif

if ($cplstsf[1]:e == "tar") then
  echo "$space cat $cplstsf | tar xf -"
         echo `cat $cplstsf | tar tf -`
               cat $cplstsf | tar xf -
else
  echo "$space cp $cplstsf ."
               cp $cplstsf .
endif

runjob chkb3c1 /dev/null "lmstr $gfargsb $ext"
chkb3c1:

echo ""
echo "$space ... Make the Omega potential self-consistent"
runjob chkb3c2 $testoutgf "$lmgf -vsc=0 $ext $gfargsb"
chkb3c2:
if ($retval == -2) goto chkb3c3
# echo "$space $lmgf -vsc=0 $ext $gfargsb >$testoutgf"
#              $lmgf -vsc=0 $ext $gfargsb >$testoutgf
# set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully."
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutgf
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutgf
endif

chkb3c3:
echo ""
echo "$space ... Make the Spectral function"
runjob chkb3c4 '>>'$testoutgf "$lmgf -vsc=-1 $ext $gfargsb"
chkb3c4:
if ($retval == -2) goto chkbe
# echo "$space $lmgf -vsc=-1 $ext $gfargsb >>$testoutgf"
#              $lmgf -vsc=-1 $ext $gfargsb >>$testoutgf
# set retval = $status
if ($retval != 0) goto cleanup
echo "$space Program lmgf returned successfully."
if ($?add0) then
 echo -n "         ..." ; $add0 $testoutgf
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testoutgf
endif

chkb4:

set Nomega = (` cat $testoutgf | grep -c repeat=F`)
set Nomegar = (`zcat $refoutgf | grep -c repeat=F`)

set Nomegai = (` cat $testoutgf | grep -c repeat=T`)
set Nomegair = (`zcat $refoutgf | grep -c repeat=T`)

if ($?have_mc) then
set ckspf = `mcx spf.$ext -s1/nr -rsum | grep -v rows`
set diffs = `echo $ckspf $ckspfr | mcx -nc=4 . -t -e1 x1-x2 -t -abs -csum | grep -v rows`
endif

if (! $?quiet) then
   echo ' '
   echo "$space Number of converged CPA Omaga iterations    = $Nomega"
   echo "$space                                reference    = $Nomegar"
   set diff = `echo $Nomega[1] $Nomegar[1] | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
   echo "$space                               difference    = $diff"

   echo ' '
   echo "$space Number of intermediate CPA Omaga iterations = $Nomegai"
   echo "$space                            reference        = $Nomegair"
   set diff = `echo $Nomegai[1] $Nomegair[1] | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
   echo "$space                           difference        = $diff"

   if ($?have_mc) then
   echo ' '
   echo "$space Checksum spf.$ext : " $ckspf
   echo "$space         reference : " $ckspfr
   echo "$space  sum of abs(diff) : " $diffs
   endif

echo ' '
call zdiffiles chkb6 "CPU 1 $testoutgf $refoutgf"
chkb6:

endif

call qprint chkbaca "$space ... automatic pass checks :"
chkbaca:

echo -n "$space Number of converged CPA Omaga iterations ($Nomega) compare to reference ($Nomegar)? ... "
if (`echo $Nomega $Nomegar 0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space Number of intermediate CPA Omaga iterations ($Nomegai) compare to reference ($Nomegair)? ... "
if (`echo $Nomegai $Nomegair 1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

set dstolb = 1e-5
if ($?have_mc) then
  set cksum = `echo $diffs | mcx . -max | grep -v rows`
echo -n "$space Max change in checksum of spectral function relative to reference ($cksum) less than tol ($dstolb)? ... "
if (`echo $cksum 0 $dstolb | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

endif


echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 11 PASSED ($ext)"
    set jobspassed = ($jobspassed 11)
else
    echo "$space test 11 FAILED ($ext)"
    set failed = ($failed 11)
endif

chkbe:

# --- Summary ---
echo ' '
if ($?clean) then
    exit 0
else if ($#failed <= 1) then
    echo "$space all tests ($jobspassed) PASSED ($ext)"
    echo " "
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
    set retval = -2; if ($?noexec) goto $quitjob
    echo " "
    $callarg
    set retval = $status
    if ($retval != 0) goto cleanup
    goto $quitjob
  endif

  if (`echo $outfile | awk '{print substr($1,1,2)}'` == '>>') then
    set appfile = `echo $outfile | awk '{print substr($1,3)}'`
    echo "$space $callarg  >> $appfile"
    set retval = -2; if ($?noexec) goto $quitjob
    echo "$space $callarg" >> $appfile
    $callarg >> $appfile
    set retval = $status
  else
    echo "$space $callarg  > $outfile"
    set retval = -2; if ($?noexec) goto $quitjob
    echo "$space $callarg" > $outfile
    $callarg >> $outfile
    set retval = $status
  endif
  if ($retval != 0) goto cleanup
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
if ($found == 'no' && $make_path != "no") then
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

# ---------------- runrdcmd --------------
exit
runrdcmd:
  set quitjob=$retcall
  if ($outfile == ".") then
    echo "$space Invoking rdcmd will execute the following job(s):"
    $rdcmd -f:$rdcmdfmt --n $callarg
    echo "$space $rdcmd '-f:rdcmd:%2f' $callarg"
#                  $rdcmd '-f:rdcmd:%2f' $callarg
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
      echo "$space Invoking rdcmd will execute the following job(s):"
      $rdcmd -f:$rdcmdfmt --n $callarg
      echo "$space $rdcmd '-f:#rdcmd:%2f' $callarg  >& $outfile"
      $rdcmd '-f:rdcmd:%2f' $callarg >& $outfile
      set retval = $status
      if ($retval == 0) then
        echo "$space Job(s) completed successfully; output in $outfile"
      endif
    endif
  endif

  if ($retval == 0) then
    if ($?add0) then
      echo -n "         ..." ; $add0 $outfile
    endif
    goto $quitjob
  else
    echo "$space ...oops... the following command returned with nonzero exit status:"
    echo -n "$space   "
    grep rdcmd: $outfile | tail -1 | sed 's/rdcmd:  //'
    goto cleanup
  endif

# ---------------- cleanup --------------
exit
cleanup:
  if ($retval != 0) echo "$space job returned with error status $retval"
  if ($retval != 0) echo "$space ... $testfile aborting"
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
    set nend = `grep $endstr $files[1] | wc | awk '{print $1}'`
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

# ---------------- zcmpmfiles_res_0 --------------
# Compares two files, stripping all but numerical fields.
# Checks for max absolute difference and unsets $passvar if difference<$tol
# Files with .gz or .Z extensions are assumed to be gzipped.
# usage: zcmpnfiles_res_0 retcall keyword testvar tol passvar ndig srcfile reffile
# See also zcmpmfiles_res_tol, which accomplished the same thing but extra argument nlines
#   retcall      : return to this point in script on exit
#   keyword      : label (for printout)
#   tol          : tolerance in maximum allowed deviation
#   passvar      : (output) $passvar is unset if |testvar|<tol
#   ndig         : number of digits numbers in file are stripped to
#   srcfile      : first file to compare
#   reffile      : second file to compare
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

# ---------------- zcmpmfiles_res_mc --------------
# Compares numerical arguments in two files, using mcx -cmpf
# If mcx is not installed, calls zcmpnfiles_res_tol with the same arguments but the last
# usage: zcmpnfiles_res_mc retcall keyword tol passvar ndig srcfile reffile nlines count
#   retcall      : return to this point in script on exit
#   keyword      : label (for printout)
#   tol          : tolerance in maximum allowed deviation
#   passvar      : $passvar is unset if numerical value of word > tol, count occurences
#   ndig         : number of digits numbers in file are stripped to
#   srcfile      : first file to compare
#   reffile      : second file to compare
#   nlines       : number of lines to compare. Use 0 to include all lines.
#   count        : maximum number of deviations to permit before unsetting passvar.
#
# Checks for max absolute difference and unsets $passvar if difference<$tol fewer than count times
# If variable exclude is set, it is used as a regular expression to exclude lines containing it
# Example:
# zcmpmfiles_res_mc chk1ck "Max deviation in bnds.$ext from reference" $bndstol pass 4 bnds.$ext $testdir/bnds.$ext.gz 0 1
exit
zcmpmfiles_res_mc:
  set quitjobl=$retcall
# echo $retcall $keyword $tol $?passvar $ndig $srcfile $reffile $nlines $count

  if (! $?mcx) then
    echo "$space OOPS! missing mcx calculator ... skipping test $keyword"
    goto $quitjobl
  endif

  set fn1 = $tmpdir/tmp_compnfile_1
  set fn2 = $tmpdir/tmp_compnfile_2

  if ($?exclude) then
    cat   $srcfile | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" > $fn1
    cat   $reffile | sed 's:\([1-9]\)-:\1 -:g' | egrep -v "$exclude" > $fn2
  else
    cat   $srcfile | sed 's:\([1-9]\)-:\1 -:g' > $fn1
    cat   $reffile | sed 's:\([1-9]\)-:\1 -:g' > $fn2
  endif

  set nl
  if ($nlines != 0) set nl = "~ln=$nlines"
# echo $mcx -cmpf$nl~fn1=$fn1~fn2=$fn2~tol=1d-$ndig~max~verb
  set retval = `$mcx -cmpf$nl~fn1=$fn1~fn2=$fn2~tol=1d-15~max`
# echo retval $retval
  set ndiff = `$mcx -cmpf$nl~fn1=$fn1~fn2=$fn2~tol=1d-$ndig~ndiff`
# echo ndiff $ndiff
# echo $count

  echo -n "$space $keyword ($retval) within tol ($tol)? ... "
  if (`echo $retval 0 | awk -v tol=$tol '{print ($1 <= tol)}'`) then
    echo yes
  else if ($count > 0) then
    echo -n "no ... fewer than $count occurences? ..."
    if (`echo $ndiff 0 | awk -v count=$count '{print ($1 <= count)}'`) then
      echo yes
    else
      echo no
      unset $passvar
    endif
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
#   nlines=#     : (optional)
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
    call zcmpmfiles zcmpmfilesx "$ndig $srcfile $reffile"
  else
    call zcmpmfiles zcmpmfilesx "nlines=$nlines $ndig $srcfile $reffile"
  endif
zcmpmfilesx:
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

# wc $fn1 $fn2

  set retval = `diff -y --width=300 $fn1 $fn2 | grep '|' | awk -v top=0 '{n=split($0,a,"|"); n1=split(a[1],b1); n2=split(a[2],b2); { j=0; while (j++ < n1) if (j <= n1 && j<=n2) {x = (b1[j]-b2[j])>0?(b1[j]-b2[j]):(b2[j]-b1[j]); top = (top-x)>0?top:x; }}} END {printf "%12.4e\n", top}'`
  rm -f $fn1 $fn2
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

# ---------------- get_resf --------------
# Extracts one element of a line in files $testout and $refout containing a keyword.
# Variables testout and refout point to file names and must be set beforehand ($refout is gzipped file)
# usage: get_resf retcall testvar refvar keyword arg_number occur_number sed_strn
#   Variables testout and refout referring to file names must be set
#   testvar      : put result from file $testout into this variable
#   refvar       : put result from file $refout (compressed) into this variable
#   keyword      : string line must contain
#   arg_number   : extracts $arg_number'th entry in line, as defined by awk
#   occur_number : argument from $occur_number'th line; if zero, use last line
#   sed_strn     : purge this string from result before assigning
exit
get_resf:
  set quitjob=$retcall
# echo $retcall $testvar $refvar $keyword $arg_number $occur_number $sed_strn
#   echo \
#   "set $testvar =" '`grep' '"'"$keyword"'"' $testout '| awk -v ncnt=0 -v num='$arg_number '-v count='$occur_number "'{ncnt+=1; if (ncnt==count || count == 0) {print "'$'"num}}'" '| sed' '"'"s/$sed_strn//"'"' ' | tail -1`'
  set $testvar = `grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  set $refvar = `zcat $refout | grep "$keyword" | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//" | tail -1`
  goto $quitjob

# ---------------- cnvt_d_fmt --------------
# converts exponential format #.##D## or #.##d## to #.##E##
# usage: cnvt_d_fmt retcall testvar testval
exit
cnvt_d_fmt:
  set quitjob = $retcall
  set $testvar = `echo $testval | sed s/D/E/ | sed s/d/E/`
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

all:
  set allargs = (`echo $argv |  sed s/--all//g`)
  set mater_lst = (co gas fccfe mnpt mnn ni cdte eras fepd fev femnpt nife)
  echo "$space $testfile : checks for materials :"  $mater_lst
  set jobargs joblist
  while ($#allargs)
    set jobargs = ($jobargs `echo $allargs[1] | awk '{if (match($1,"-") && (RSTART == 1)) print $1}'`)
    set joblist = ($joblist `echo $allargs[1] | awk '{if (match($1,"-") && (RSTART == 1)) { } else {print $1}}'`)
    shift allargs
  end
  set failed; unset failed
  foreach i ($mater_lst)
     echo  ' '
     echo  ' ' $testfile $jobargs $i $joblist
               $testfile $jobargs $i $joblist
     set retval = $status
     if ($retval != 0) then
       if (! $?failed) set failed
       set failed = ($failed $i)
#        echo " $testfile : failed test $i ... aborting"
#        exit -1
     endif
  end

  if ($?failed) then
    echo "$space $testfile : these scalar relativistic tests FAILED:  $failed"
    exit -1
  endif
  echo "$space $testfile : all scalar relativistic tests PASSED ($mater_lst)"

  echo " "
  echo " $testfile : fully relativistic checks: invoke "  $testdir/test.frgf $jobargs --all
  echo  " " $testdir/test.frgf $jobargs --all
            $testdir/test.frgf $jobargs --all
  set retval = $status
  if ($retval != 0) then
    echo " $testdir/test.frgf: fully relativistic tests FAILED (ni)"
    exit -1
  endif
  exit

# ---------------- List tests --------------
showtests:
cat <<EOF
  Usage: invoke with:
    test.gf material-name job-list

   Material:
         co: elemental hexagonal metal
        gas: insulator
      fccfe: noncollinear checks. The spinor symmetrization of G is suppressed
       mnpt: various kinds of tests, incl division of H into lower,intermediate, high blocks
        mnn: various kinds of tests, incl exchange interactions
         ni: tests magnetic exchange coupling
       cdte: tests screened exchange self-energy
       eras: tests LDA+U GF.  In the noncollinear test, spinor part is not symmetrized.
       fepd: tests disordered local moments
        fev: tests chemical CPA
     femnpt: tests chemical CPA
       nife: tests noncollinear spin statics: relaxation of Euler angles
         ... the following are not calculated in the checks gf/test/test.gf --all
       fept: tests noncollinear code to compare to equivalent collinear structure
       fe2b: tests CPA with spin orbit coupling
       fevg: same as fev, but in gamma representation (CPA Omega files are supplied in gamma)
     bi2te3: four monolayers of Bi2Te3 separated by a vacuum.  Compare SO-coupled spectral function to lm bands

  jobs:   1: standard pass, 2nd order GF
          2: standard pass, 3rd order GF
          3: generate site density matrix
          4: iterates GF to self-consistency
          5: magnetic exchange interactions
          6: tests linear response
          7: tests of noncollinear branch
          8: tests linear response to get self-energy
          9: tests CPA (disordered Local Moments)
         10: tests CPA (chemical CPA)
         11: tests spectral function in CPA (chemical CPA)
EOF
exit

# ---------------- usage: --------------
usage:
cat <<EOF
 usage: $0 [switches] [material] [testcase-list | --all]
        e.g., "$0 co 1 2"
        'material' can be one of : co gas fccfe mnpt mnn ni cdte eras fepd fev nife fe
        If 'material' is missing, test.gf uses co
        Switches:
        --list       lists the tests you can run
        --no-iactive runs tests without prompting user
        --quiet      runs tests with minimal output and without prompting user
        --all        run through a default list of test cases
                     Invoke test.gf as :  "$0 [switches] --all"
        --noplot     skip any steps that generate a plot'
        --clean      clean up files generated by this script
        --gamma=#    Run lmgf in some other screening transformation : # can be one of 0 1 2 4 5
        --add0       add suppressed or leading zeros in output for real nubmers \`.nnn'
        --poszer     strips (-) sign from numbers represented as 0
        --whichexec  prints out which lmgf executable it finds in path and exits
        --MPIK[=#]   run lmgf in MPIK mode
                     If # is given, it specifies the number of processors.
                     The MPI job is run as
                     mpirun -n # lmgf ...

EOF
exit -1
