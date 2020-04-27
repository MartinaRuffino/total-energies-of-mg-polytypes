#!/bin/tcsh -f

# This file is a shell script testing operation of noncollinear ASA package

alias call 'set retcall = \!\!:2 ; set callarg = \!\!:3 ; goto \!\!:1'
alias runjob 'set retcall = \!\!:1; set outfile = \!\!:2 ; set callarg = \!\!:3 ; goto runjob'
alias runrdcmd 'set retcall = \!\!:1; set rdcmdfmt = \!\!:2 ; set outfile = \!\!:3 ; set callarg = \!\!:4 ; goto runrdcmd'
alias findcmd  'set retcall = \!\!:1 ; set prog_cmd = \!\!:2 ; set path_name = \!\!:3 ; set make_path = \!\!:4 ; goto findcmd'
alias query 'set retcall = \!\!:1 ; set retcall2 = \!\!:2 ; set callarg = \!\!:3 ; goto query'
alias compare_res 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set refvar = \!\!:4 ; set tol = \!\!:5 ; set passvar = \!\!:6 ; goto compare_res'
alias compare_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set testvar = \!\!:3 ; set tol = \!\!:4 ; set passvar = \!\!:5 ; goto compare_res_0'
alias compare_resf 'set retcall = \!\!:1; set testvar = \!\!:2 ; set refvar = \!\!:3 ; set keyword = \!\!:4  ; set arg_number = \!\!:5 ; set occur_number = \!\!:6 ; set sed_strn = \!\!:7 ; goto compare_resf'
alias zcmpmfiles_res_0 'set retcall = \!\!:1; set keyword = \!\!:2 ; set tol = \!\!:3 ; set passvar = \!\!:4 ; set ndig = \!\!:5 ; set srcfile = \!\!:6 ; set reffile = \!\!:7 ; goto zcmpmfiles_res_0 '
alias cnvt_d_fmt  'set retcall = \!\!:1; set testvar = \!\!:2 ; set testval = \!\!:3 ; goto cnvt_d_fmt'

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

#alias mpix mpirun
set eojob

# Prepend current working-directory, top-level and related dir to path
set path = ($cwd $topdir $topdir/testing $path)

# --- Pick off switches ---
while (`echo $1 | sed -e 's/\(.\).*/\1/' `  ==  "-")

  set arg1 = $1; shift
  if ($?verb) echo test.pgf: parsing switch $arg1
  switch ($arg1)
    case "--quiet":
      set quiet
      unset slow
      breaksw
    case "--clean":
      set clean
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
      set MPIK
      breaksw
    case "--no-iact*":
      unset slow
      breaksw
    case "--list":
      goto showtests
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

    case "--all":
      set allargs = `echo $allargs | sed s/--all//g | sed -e 's/\([0-9][0-9]*\)//g' | sed -e 's/-add/-add0/g'`
      if ($?MPIK) then
        echo "$space test.so not set up for MPI job ... skipping"
      else
        echo "    ... invoke $testdir/test.so for spin-orbit coupling checks ..."
        echo "$space $testdir/test.so $allargs"
        $testdir/test.so $allargs
        set retval = $status
        if ($retval != 0) then
          echo "$space $testfile : SO tests FAILED ... aborting"
          exit -1
        endif
      endif

      set joblist
      while (`echo $1 | sed -e 's/\([0-9][0-9]*\)/-/'`  ==  "-")
        set joblist = ($joblist $1)
        shift
      end
      if ($joblist[1] == "") set joblist = (1 2 3 4 5 6 7 8)
      set pass
      set failed
      echo "$space invoke $testfile $allargs $joblist"
                          $testfile $allargs $joblist
      set retval = $status
      if ($retval != 0) then
        unset pass
#       set failed = ($failed $i)
#         echo "$space $testfile : some tests FAILED ... abort"
        exit -1
      endif
      exit

    default:
      echo unrecognized switch $arg1
      goto usage
  endsw

end

set joblist = (`echo $argv | sed s/fe2//`)
if ($#joblist == 0 ) set joblist = (1 2 3 4 5 6 7 8)

echo $joblist | grep 1 >/dev/null
if ($status) goto chk1e

cat <<EOF

         --- Test case 1: 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests that equilibrium 3k structure is self-consistent.
         Test DESTROYS or overwrites *.fccfe
EOF
if ($?quiet) then
else
cat <<EOF
         Inspect input file for other magnetic structures to test.
EOF
endif
set pass
set refout=$testdir/out.fccfe.3k.equil testout=out.fccfe
echo ' '
query chk11 chk1e 'run this test'
chk11:
# ... Look for executables
findcmd chk11a rdcmd "$path" "$topdir"
chk11a:
findcmd chk11b lm "$path" "$topdir"
chk11b:
findcmd chk11c lmstr "$path" "$topdir"
chk11c:

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f *.fccfe"
             touch ctrl.fccfe
             rm -f *.fccfe
if ($?clean) then
  echo "$space rm -f $testout $testout-dos"
               rm -f $testout $testout-dos
  goto chk1e
endif
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .
runjob chk12 $testout "lmstr fccfe"
chk12:
runjob chk13 $testout "lm -vnit=0 -vsdyn=t fccfe --no-iactive"
chk13:
runjob chk14 $testout "lm -vnit=1 -vsdyn=t fccfe --no-iactive"
chk14:
echo "$space Program lm returned successfully."
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set dq = `cat $testout | grep SV: | tail -1 | awk '{print $3}' | sed 's/D/e/'`
awk '{if ($1 == "MAGTRQ:") {getline;{getline;if ($3*$3 > 1e-12) {exit -1}};{getline;if ($3*$3 > 1e-12) {exit -1}};{getline;if ($3*$3 > 1e-12) {exit -1}};{getline;if ($3*$3 > 1e-12) {exit -1}}}}' $testout
set forcesmall = (! $status)

if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

if ($?quiet) then
  else
  echo ' '
  echo "$space ...show magnetic forces"
  awk '{if ($1 == "MAGTRQ:") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' $testout
  echo ' '

call showout chk15a SV
chk15a:

endif # quiet

call zdiffiles chk15 "CPU 1 $testout $refout"
chk15:

set refout=$testdir/out.fccfe.3k.ss.equil testout=out.fccfe
echo ' '
echo "$space ... repeat, for 3k+SS"
echo "$space rm -f *.fccfe"
             touch ctrl.fccfe
             rm -f *.fccfe
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .
runjob chk16a $testout "lmstr fccfe"
chk16a:
runjob chk16b $testout "lm -vnit=0 fccfe -vqss=.5 -vsdyn=f fccfe --no-iactive"
chk16b:
if (! $?MPIK) then
runjob chk16c $testout "lm -vnit=1 fccfe -vqss=.5 -vsdyn=f fccfe --no-iactive"
else
runjob chk16c $testout "mpirun -n 8 lm -vmet=2 -vnit=1 fccfe -vqss=.5 -vsdyn=f fccfe --no-iactive"
endif
chk16c:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set dqss = `cat $testout | grep SV: | tail -1 | awk '{print $3}' | sed 's/D/e/'`

if ($?quiet) then
else

  call showout chk18a SV
chk18a:

  echo ' '
  call zdiffiles chk18 "CPU 1 $testout $refout"
chk18:

endif

set refout=$testdir/out.fccfe.2ss.dnf testout=out.fccfe
echo ' '
echo "$space ... Check ++-- spin structure using SS, downfolding Fe p orbitals"
if ($?quiet) then
else
cat <<EOF
         This 2-atom SS calculation will be compared to a four-atom ++-- structure,
         showing that the total energies of the two converge to the same value.

         Also the site-summed, m-resolved partial DOS are evaluated.
         Use the following commands for each partial DOS generated:
           echo 15 5 -4 4 | pldos -ef=0 -escl=13.6 -fplot -lst='1;3;5;7;9;11;13;15;17' -lst2 dos.qss.fccfe
           echo 30 5 -4 4 | pldos -ef=0 -escl=13.6 -fplot -lst='1;3;5;7;9;11;13;15;17' -lst2 dos.fccfe
         and you can confirm that the ++-- DOS is twice the SS DOS, i.e.
         the DOS in files dosp.dat and dosp2.dat generated by the two pldos commands
         differ by a factor of almost exactly 2.

EOF
endif

echo "$space rm -f *.fccfe"
             touch ctrl.fccfe
             rm -f *.fccfe
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .

runjob chk19a $testout "lmstr fccfe -cstrx=2ss"
chk19a:
runjob chk19b $testout "lm -vnit=0 fccfe -vqss=0.5 -cstrx=2ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2"
chk19b:
#  echo "$space cp a.fccfe a1.fccfe"
#               cp a.fccfe a1.fccfe
#  echo "$space cp a.fccfe a2.fccfe"
#               cp a.fccfe a2.fccfe
if (! $?MPIK) then
runjob chk19c $testout "lm -vnit=1 fccfe -vqss=0.5 -cstrx=2ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2 --noinv -vmet=2 --pdos~group=1:2"
else
runjob chk19c $testout "mpirun -n 8 lm -vnit=1 fccfe -vqss=0.5 -cstrx=2ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2 --noinv -vmet=2 --pdos~group=1:2"
endif
chk19c:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

runjob chk19c2 $testout-dos "lmdos -vnit=1 fccfe -vqss=0.5 -cstrx=2ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2 --noinv -vmet=2 --pdos~group=1:2 --dos:npts=501:window=-1,0.5"
chk19c2:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout-dos
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout-dos
endif
echo "$space cp dos.fccfe dos.qss.fccfe"
             cp dos.fccfe dos.qss.fccfe
if ($?add0) then
  echo -n "         ..." ; $add0 dos.qss.fccfe
else if ($?poszer) then
  echo -n "         ..." ; $poszer dos.qss.fccfe
endif

set ehf2ss = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehfref = `zcat $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
#  set ehf2ss = `cat $testout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`
#  set ehfref = `zcat $refout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`

set ehfref48 = -0.8884846
set ediff = `echo $ehf2ss $ehfref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediff48 = `echo $ehf2ss $ehfref48  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`


if ($?quiet) then
else

  echo " "
  echo "$space calculated harris energy =          $ehf2ss"
  echo "$space compare to reference (32 divisions) $ehfref"
  echo "$space difference                           $ediff"
  echo "$space compare to reference (48 divisions) $ehfref48"
  echo "$space difference                           $ediff48"
  echo "$space NB: compare to total energy of ++-- spin configuration calculated from 4-atom structure (next test)"

  echo ' '
  call zdiffiles chk19d "CPU 1 $testout $refout"
chk19d:
endif

set refout=$testdir/out.fccfe.4ss.dnf testout=out.fccfe
echo ' '
echo "$space ... Check ++-- spin structure using 4-atom cell, downfolding Fe p orbitals."
echo "$space     Starting moments are identical to SS. Confirm that the 4-atom energy exactly doubles that"
echo "$space     of the 2-atom SS structure, and that the two are self-consistent in the same potential."
echo "$space rm -f {a1,a2,ctrl,dos,eula,log,moms,mq,out,save,sdot,str,sv,shfac,vshft,wkp}.fccfe"
             rm -f {a1,a2,ctrl,dos,eula,log,moms,mq,out,save,sdot,str,sv,shfac,vshft,wkp}.fccfe
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .

runjob chk19e $testout "lmstr fccfe -cstrx=4ss"
chk19e:
runjob chk19f $testout "lm -vnit=0 fccfe -vqss=0.0 -cstrx=4ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2"
chk19f:
if (! $?MPIK) then
runjob chk19g $testout "lm -vnit=1 fccfe -vqss=0.0 -cstrx=4ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2 --noinv -vmet=2 --pdos~group=1:4"
else
runjob chk19g $testout "mpirun -n 8 lm -vnit=1 fccfe -vqss=0.0 -cstrx=4ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2 --noinv -vmet=2 --pdos~group=1:4"
endif
chk19g:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif
runjob chk19g2 $testout-dos "lmdos -vnit=1 fccfe -vqss=0.0 -cstrx=4ss -vsdyn=f --iactiv=no -vnk=32 -vidxp=2 --noinv -vmet=2 --pdos~group=1:4 --dos:npts=501:window=-1,0.5"
chk19g2:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout-dos
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout-dos
endif

set ehf4ss = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehfref = `zcat $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
#  set ehf4ss = `cat $testout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`
#  set ehfref = `zcat $refout | grep LM: | tail -1 | awk '{print $2}' | sed 's/ehf=//'`

set ehfref48 = -1.7769874
set ediff4 = `echo $ehf4ss $ehfref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediff48 = `echo $ehf4ss $ehfref48  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

if ($?quiet) then
else

  echo " "
  echo "$space calculated harris energy =          $ehf4ss"
  echo "$space compare to reference (32 divisions) $ehfref"
  echo "$space difference                           $ediff4"
  echo "$space compare to reference (48 divisions) $ehfref48"
  echo "$space difference                           $ediff48"
  echo "$space NB: last energy matches ++-- spin configuration calculated from 4ss structure"

  echo ' '
  call zdiffiles chk19h "CPU 1 $testout $refout"
chk19h:

endif

# ... automatic pass checks
chk1p:
echo ' '
call qprint chk1pa "$space ... automatic pass checks :"
chk1pa:

echo -n "$space rms dq (=$dq) < 1e-6 in 3k structure ? ... "
if (`echo $dq | awk '{print ($1 < 1e-6)}'`) then
  echo  yes
else
  echo no
  unset pass
endif

echo -n "$space magnetic forces < 1e-6 in 3k structure ? ... "
if ($forcesmall) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space rms dq (=$dqss) < 1e-6 in 3k+SS structure ? ... "
if (`echo $dqss | awk '{print ($1 < 1e-6)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo -n "$space difference in Harris energy for 2SS (=$ediff) < 2e-6 ? ... "
if (`echo $ediff | awk '{print ($1 < 2e-6)}'`) then
  echo yes
else
  echo no
  unset pass
endif

set ndig = 3
call zcmpnfiles chk1pb "$ndig dos.qss.fccfe $testdir/pdos.qss0.5.fccfe"
chk1pb:
echo -n "$space ... files dos.qss.fccfe and $testdir/pdos.qss0.5.fccfe equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  echo  no
  unset pass
endif

set ndig = 3
call zcmpnfiles chk1pc "$ndig dos.fccfe $testdir/pdos.qss0.0.fccfe"
chk1pc:
echo -n "$space ... files dos.fccfe and $testdir/pdos.qss0.0.fccfe equivalent to $ndig digits? ... "
if ($retval == 0) then
  echo  yes
else
  echo  no
  unset pass
endif


echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 1  PASSED"
else
    echo "$space test 1 FAILED"
    set failed = ($failed 1)
endif

chk1e:

echo $joblist | grep 2 >/dev/null
if ($status) goto chk2e
cat <<EOF

         --- Test case 2: SD in 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests spin statics using sdmod=1 for 4 atoms/cell Fe.
         Test DESTROYS or overwrites *.fccfe.
         NB: the HF and HK energy functionals don't exactly agree because
         the atom potentials are artifically averaged (GRP2=1).

EOF
set eojob=chk2p  sdmod=1  eulafile=eula.ran
goto sdjob
chk2p:
if ($?clean) then
  goto chk2e
endif
if ($?quiet) then
else
echo "$space ... automatic pass checks :"
endif

if ($?drmsqtol2 == 0) set drmsqtol2 = 1e-5
compare_res chk2cb "last iter RMS dq" $dqn $dqnr $drmsqtol2 pass
chk2cb:

if ($?drmsetol2 == 0) set detol2 = 2e-6
compare_res chk2cc "last iter ehf" $etest $eref $detol2 pass
chk2cc:

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 2 PASSED"
else
    echo "$space test 2 FAILED"
    set failed = ($failed 2)
endif

chk2e:
if ($eojob == "chk3p") goto chk3e
if ($eojob == "chk4p") goto chk4e

echo $joblist | grep 3 >/dev/null
if ($status) goto chk3e
cat <<EOF

         --- Test case 3: SD in 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests spin statics using sdmod=11 for 4 atoms/cell Fe.
         Test DESTROYS or overwrites *.fccfe
         NB: the HF and HK energy functionals don't exactly agree because
         the atom potentials are artifically averaged (GRP2=1).

EOF
set eojob=chk3p sdmod=11  eulafile=eula.ran
goto sdjob
chk3p:
if ($?clean) then
  goto chk3e
endif
if ($?quiet) then
else
echo "$space ... automatic pass checks :"
endif

if ($?drmsqtol3 == 0) set drmsqtol3 = 1e-5
compare_res chk3cb "last iter RMS dq" $dqn $dqnr $drmsqtol3 pass
chk3cb:

if ($?drmsetol3 == 0) set detol3 = 1e-5
compare_res chk3cc "last iter ehf" $etest $eref $detol3 pass
chk3cc:


echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 3 PASSED"
else
    echo "$space test 3 FAILED"
    set failed = ($failed 3)
endif

chk3e:

echo $joblist | grep 4 >/dev/null
if ($status) goto chk4e
cat <<EOF

         --- Test case 4: SD in 3k structure in fccfe ($testdir/ctrl.fccfe)  ---
         Tests spin statics using sdmod=1011 for 4 atoms/cell Fe.
         Test also uses l-dependent Euler angles.
         Test DESTROYS or overwrites *.fccfe
EOF
set eojob=chk4p sdmod=1011 eulafile=eula-l.ran
goto sdjob
chk4p:
if ($?clean) then
  goto chk4e
endif
echo ' '
if ($?quiet) then
else
echo "$space ... automatic pass checks :"
endif

if ($?drmsmtol4 == 0) set drmsmtol4 = 1e-5
compare_res chk4cb "last iter RMS change in euler angles" $dmn $dmnr $drmsmtol4 pass
chk4cb:

if ($?drmsetol4 == 0) set detol4 = 1e-5
compare_res chk4cc "last iter ehf" $etest $eref $detol4 pass
chk4cc:

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 4 PASSED"
else
    echo "$space test 4 FAILED"
    set failed = ($failed 4)
endif

chk4e:

echo $joblist | grep 5 >/dev/null
if ($status) goto chk5e
cat <<EOF

         --- Test case 5: Applied magnetic field ($testdir/ctrl.fe)  ---
         Tests the application of an external magnetic field, bcc Fe, in a 2 atom superlattice.
         (Test DESTROYS or overwrites *.fe)

EOF
set pass
set refout=$testdir/out.fe.matchcoll testout=out.fe
echo ' '
query chk51 chk5e 'run this test'
chk51:
# ... Look for executables
findcmd chk51a rdcmd "$path" "$topdir"
chk51a:
findcmd chk51b lm "$path" "$topdir"
chk51b:
findcmd chk51c lmstr "$path" "$topdir"
chk51c:

echo " "
echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f {ctrl,mixm,moms,wkp,site2,log,fe,fe2,eula,bfield,save,sdot,str,sv,shfac}.fe"
             rm -f {ctrl,mixm,moms,wkp,site2,log,fe,fe2,eula,bfield,save,sdot,str,sv,shfac}.fe
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk5e
endif
echo "$space cp $testdir/{ctrl,site2}.fe ."
             cp $testdir/{ctrl,site2}.fe .
echo "$space touch bfield.fe"
             touch bfield.fe
runjob chk52a $testout "lmstr -vnit=0 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1"
chk52a:
runjob chk52b $testout "lm -vnit=0 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1"
chk52b:

cat <<EOF

         ... Verify that collinear and noncollinear branches produce same result.
         A collinear calculation is compared to a noncollinear one with both
         sites at the same angle.  In the latter case amplitude of the moment
         on one site (|M|=2.290999) is half the total collinear moment.

EOF
cat >eula.fe <<in
.1 .2 .3
.1 .2 .3
in
echo "$space Use the following Euler angles file:"
cat eula.fe
runjob chk53a "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=0 -vbf=0 --quit=band"
chk53a:
if (! $?MPIK) then
runjob chk53b ">>$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=0 --quit=band"
else
runjob chk53b ">>$testout" "mpirun -n 8 lm -vmet=3 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=0 --quit=band"
chk53b:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set sumevc = `cat $testout | grep 'sumev=' | sed -n 1,1p | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set sumev0 = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set sumevr = `zcat $refout | grep 'sumev=' | sed -n 1,1p | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf0   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk0   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehfr   = `zcat $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `zcat $refout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif
set mz0 = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $9}'`

echo "$space    collinear band structure energy = $sumevc"
echo "$space noncollinear band structure energy = $sumev0"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumevc $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf0"
echo "$space          reference   Harris energy = $ehfr"
set ehfdiff = `echo $ehf0 $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk0"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk0 $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf0 $ehk0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"
echo "$space       noncollinear magnetic moment = $mz0"

echo ' '
call qprint chk54a "$space ... automatic pass checks :"
chk54a:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev0 $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf0 $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

call zdiffiles chk55a "CPU 1 $testout $refout"
chk55a:

# ... Stoner susceptibility
set refout=$testdir/out.fe.stonerI testout=out.fe
echo ' '

cat <<EOF

         ... Calculate bare longitudinal (Stoner) susceptibility.
         Test applies a field along the z axis to one site; a single band pass is made.

EOF
cat >eula.fe <<in
0 0 0
0 0 0
in
echo "$space using Euler angles file:"
cat eula.fe
cat >bfield.fe <<in
0 0 -.001/2
0 0 0
in
echo "$space using B-field file:"
cat bfield.fe

zcat $refout >$testout
set mzref  = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $5}'`
set sumevr = `cat $testout | grep 'sumev=' | sed -n 1,1p | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehfr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?MPIK) then
runjob chk56a "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
else
runjob chk56a "$testout" "mpirun -n 8 lm -vmet=3 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
endif
chk56a:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set mz     = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $5}'`
set mz2    = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "2") print $5}'`
set sumev  = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif

echo "$space noncollinear band structure energy = $sumev"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumev $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf"
echo "$space         reference    Harris energy = $ehfr"
set ehfdiff = `echo $ehf $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf $ehk | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"
echo "$space output moment on atom 1            = $mz"
echo "$space output moment on atom 2            = $mz2"
set mzdiff = `echo $mz $mz0 | awk '{{k=($1-$2)} print k}'`
echo "$space induced moment site 1 M(B)-M(B=0)  =  $mzdiff"
set mzdiff = `echo $mz2 $mz0 | awk '{{k=($1-$2)} print k}'`
echo "$space induced moment site 2 M(B)-M(B=0)  = $mzdiff"
echo "$space reference output moment            = $mzref"
set mzdiff = `echo $mz $mzref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $mzdiff"

echo ' '
call qprint chk57a "$space ... automatic pass checks :"
chk57a:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
call zdiffiles chk58a "CPU 1 $testout $refout"
chk58a:

# ... transverse susceptibility, spins aligned along x
set refout=$testdir/out.fe.bfieldx testout=out.fe
cat <<EOF

         ... Calculate bare transverse susceptibility.
         Test applies a field along the x axis to one site; a single band pass is made.

EOF
cat >eula.fe <<in
0 0 0
0 0 0
in
echo "$space using Euler angles file:"
cat eula.fe
cat >bfield.fe <<in
-.001/2 0 0
0 0 0
in
echo "$space using B-field file:"
cat bfield.fe

zcat $refout >$testout
set mzref  = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $3}'`
set sumevr = `cat $testout | grep 'sumev=' | sed -n 1,1p | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehfr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?MPIK) then
runjob chk59a "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
else
runjob chk59a "$testout" "mpirun -n 8 lm -vmet=3 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
endif
chk59a:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set mz     = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $3}'`
set mz2    = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "2") print $3}'`
set sumev  = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif

echo "$space noncollinear band structure energy = $sumev"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumev $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf"
echo "$space         reference    Harris energy = $ehfr"
set ehfdiff = `echo $ehf $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf $ehk | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"
echo "$space output (induced) moment on atom 1  = $mz"
echo "$space reference output moment            = $mzref"
set mzdiff = `echo $mz $mzref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $mzdiff"
echo "$space output (induced) moment on atom 2  = $mz2"

echo ' '
call qprint chk5aa "$space ... automatic pass checks :"
chk5aa:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
call zdiffiles chk5ba "CPU 1 $testout $refout"
chk5ba:

# ... longitudinal susceptibility, spins aligned along x
set refout=$testdir/out.fe.stonerIx testout=out.fe
cat <<EOF

         ... Bare longitudinal (Stoner) susceptibility for spins aligned along the x axis.
         Test applies a field along the x axis to one site; make a single band pass.
         Try comparing the results of this calculation to the longitudinal susceptibility
         for spins parallel to the z axis:
            zdiff nc/test/out.fe.stonerI nc/test/out.fe.stonerIx

EOF

cat >eula.fe <<in
0 pi/2 0
0 pi/2 0
in
echo "$space using Euler angles file:"
cat eula.fe
cat >bfield.fe <<in
-.001/2 0 0
0 0 0
in
echo "$space using B-field file:"
cat bfield.fe

zcat $refout >$testout
set mzref  = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $3}'`
set sumevr = `cat $testout | grep 'sumev=' | sed -n 1,1p | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehfr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?MPIK) then
runjob chk5ca "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
else
runjob chk5ca "$testout" "mpirun -n 8 lm -vmet=3 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
endif
chk5ca:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set mz     = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $3}'`
set mz2    = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "2") print $3}'`
set sumev  = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif

echo "$space noncollinear band structure energy = $sumev"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumev $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf"
echo "$space         reference    Harris energy = $ehfr"
set ehfdiff = `echo $ehf $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf $ehk | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"

echo "$space output moment on atom 1            = $mz"
echo "$space output moment on atom 2            = $mz2"
set mzdiff = `echo $mz $mz0 | awk '{{k=($1-$2)} print k}'`
echo "$space induced moment site 1 M(B)-M(B=0)  =  $mzdiff"
set mzdiff = `echo $mz2 $mz0 | awk '{{k=($1-$2)} print k}'`
echo "$space induced moment site 2 M(B)-M(B=0)  = $mzdiff"
echo "$space reference output moment            = $mzref"
set mzdiff = `echo $mz $mzref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $mzdiff"



echo ' '
call qprint chk5da "$space ... automatic pass checks :"
chk5da:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
call zdiffiles chk5ea "CPU 1 $testout $refout"
chk5ea:

# ... transverse susceptibility, spins aligned along x
set refout=$testdir/out.fe.bfieldxy testout=out.fe
cat <<EOF

         ... Calculate bare transverse susceptibility, with initial spins aligned along x.
         Test applies a field along the y axis to one site; a single band pass is made.
         Try comparing the results of this calculation to the transverse susceptibility
         for spins parallel to the z axis:
            zdiff nc/test/out.fe.bfieldxy nc/test/out.fe.bfieldx

EOF
cat >eula.fe <<in
0 pi/2 0
0 pi/2 0
in
echo "$space using Euler angles file:"
cat eula.fe
cat >bfield.fe <<in
0 -.001/2 0
0 0 0
in
echo "$space using B-field file:"
cat bfield.fe

zcat $refout >$testout
set mzref  = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $4}'`
set sumevr = `cat $testout | grep 'sumev=' | sed -n 1,1p | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehfr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehkr   = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?MPIK) then
runjob chk5fa "$testout" "lm fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
else
runjob chk5fa "$testout" "mpirun -n 8 lm -vmet=3 fe --no-iactiv -vfile=2 -vnk=8 -vnl=3 -vnc=1 -vbf=1 --quit=band"
endif
chk5fa:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set mz     = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "1") print $4}'`
set mz2    = `extract-lines  --quiet 'density matrix' "ehf=" 1 $testout | awk '{if ($1 == "2") print $4}'`
set sumev  = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"sumev=[^ ]*"); print substr($0,RSTART+6,RLENGTH-6)}'`
set ehf    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk    = `cat $testout | grep 'sumev=' | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`

if (! $?quiet) then
echo ' '
echo "$space Magnetic moments from file $testout":
extract-lines  --quiet 'density matrix' "ehf=" 1 $testout
endif

echo "$space noncollinear band structure energy = $sumev"
echo "$space    reference band structure energy = $sumevr"
set sevdiff = `echo $sumev $sumevr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $sevdiff"
echo "$space         noncollinear Harris energy = $ehf"
echo "$space         reference    Harris energy = $ehfr"
set ehfdiff = `echo $ehf $ehfr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehfdiff"
echo "$space      noncollinear Kohn-Sham energy = $ehk"
echo "$space      reference    Kohn-Sham energy = $ehkr"
set ehkdiff = `echo $ehk $ehkr  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $ehkdiff"
set ediff = `echo $ehf $ehk | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space            Harris - K-S difference = $ediff"
echo "$space output (induced) moment on atom 1  = $mz"
echo "$space reference output moment            = $mzref"
set mzdiff = `echo $mz $mzref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space                         difference = $mzdiff"
echo "$space output (induced) moment on atom 2  = $mz2"

echo ' '
call qprint chk5ga "$space ... automatic pass checks :"
chk5ga:
set sevtol1 = 1e-6
if ($?sevtol1 == 0) set sevtol1 = 1e-6
echo -n "$space sumev within tol $sevtol1 of reference? ... "
if (`echo $sumev $sumevr $sevtol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif
set ehftol = 1e-6
if ($?ehftol1 == 0) set ehftol1 = 1e-6
echo -n "$space   ehf within tol $ehftol1 of reference? ... "
if (`echo $ehf $ehfr $ehftol1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=$3)}'`) then
  echo yes
else
  echo no
  unset pass
endif

echo ' '
call zdiffiles chk5ha "CPU 1 $testout $refout"
chk5ha:

echo ' '
if ($?clean) then
else if ($?pass) then
    echo "$space test 5 PASSED"
else
    echo "$space test 5 FAILED"
    set failed = ($failed 5)
endif

chk5e:

echo $joblist | grep 6 >/dev/null
if ($status) goto chk6e

set ext = fe2

cat <<EOF

         --- Test case 6: Exchange interactions  ---
         Makes a small rotation of a spin in the presence of an external
         constraining field to zero out the d component of the spin density
         matrix.  Calculation compared to rotation in the absence of
         constraining field.  NB: a careful calculation would need to
         pay attention to k-point convergence.

EOF

if ($ext == "fe2") then
cat <<EOF
         Test case Fe2 consists of 2 Fe atoms (simple cubic unit cell).
         One atom is rotated by -theta/2 about z; the other is rotated by theta/2.
         Constraining field is taken to be along x axis.

EOF
endif

set pass

set refout=$testdir/out.$ext.coll testout=out.$ext
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk6e
endif
echo ' '
query chk61 chk6e 'run this test'
chk61:
# ... Look for executables
findcmd chk61a rdcmd "$path" "$topdir"
chk61a:
findcmd chk61b lm "$path" "$topdir"
chk61b:
findcmd chk61c lmstr "$path" "$topdir"
chk61c:

findcmd chk61d rdfile "$path" "no"
chk61d:
if ("$found" == "no") then
  echo "$space ... no rdfile in path ... skipping test"
  goto chk6e
endif

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f {a,a2,eula0,b0,ctrl,site2,mixm,moms,log,wkp,bfield,eula,save,sdot,str,sv,shfac,vshft}.$ext"
             rm -f {a,a2,eula0,b0,ctrl,site2,mixm,moms,log,wkp,bfield,eula,save,sdot,str,sv,shfac,vshft}.$ext
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto chk6e
endif
echo "$space cp $testdir/{eula0,b0,ctrl,site2}.$ext ."
             cp $testdir/{eula0,b0,ctrl,site2}.$ext .

echo "$space ... generate total energy for collinear case"
runrdcmd chk62a %11f $testout "-cat:JOBCOLL --noerr ctrl.$ext"
chk62a:
echo -n "$space shortening output ... "
echo "catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
      catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext
call zdiffiles chk62b "CPU 2 $testout $refout"
chk62b:

set ehf0ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk0ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehf0    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk0    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ediffhk = `echo $ehk0 $ehk0ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediffhf = `echo $ehf0 $ehf0ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`


if (! $?quiet) then
echo ' '
echo "$space collinear ehk = $ehk0"
echo "$space reference     = $ehk0ref"
echo "$space difference    = $ediffhk"
echo ' '

echo "$space collinear ehf = $ehf0"
echo "$space     reference = $ehf0ref"
echo "$space difference    = $ediffhf"
echo ' '
endif

compare_res chk62c "ehk" $ehk0 $ehk0ref 1e-6 pass
chk62c:
compare_res chk62d "ehf" $ehf0 $ehf0ref 1e-6 pass
chk62d:

echo " "
echo "$space ... rotation in the absence of constraining field"
set refout=$testdir/out.$ext.nc1
echo "$space rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log,shfac}.fe2"
             rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log,shfac}.fe2
echo "$space cp $testdir/{eula0,b0,ctrl,site2}.fe2 ."
             cp $testdir/{eula0,b0,ctrl,site2}.fe2 .
runrdcmd chk63a %11f $testout "-cat:JOBNC1 --noerr ctrl.$ext"
chk63a:
echo -n "$space shortening output ... "
echo "catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
      catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext
call zdiffiles chk63b "CPU 1 $testout $refout"
chk63b:

set ehf1ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk1ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehf1    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk1    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ediffhk = `echo $ehk1 $ehk1ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediffhf = `echo $ehf1 $ehf1ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

set mx1in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
set mx1dout  = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $3}'`
set mx1doutl = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $6}'`
set mx1out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
zcat $refout >$testout{}~
set rmx1in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
set rmx1dout  = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $3}'`
set rmx1doutl = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $6}'`
set rmx1out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
rm $testout{}~

if (! $?quiet) then

echo "$space input  Mx     = $mx1in"
echo "$space output Mx     = $mx1out"
echo "$space output Mx(d)  = $mx1dout"
echo "$space in loc. coord = $mx1doutl"
echo "$space reference     = $rmx1doutl"

echo ' '
echo "$space ehk           = $ehk1"
echo "$space reference     = $ehk1ref"
echo "$space difference    = $ediffhk"
echo "$space ehk-ehk(coll) = "`echo $ehk1 $ehk0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

echo ' '
echo "$space ehf           = $ehf1"
echo "$space reference     = $ehf1ref"
echo "$space difference    = $ediffhf"
echo "$space ehf-ehf(coll) = "`echo $ehf1 $ehf0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo ' '
endif

compare_res chk63c "ehk" $ehk1 $ehk1ref 1e-6 pass
chk63c:
compare_res chk63d "ehf" $ehf1 $ehf1ref 1e-6 pass
chk63d:
compare_res chk63e "Mx" $mx1out $rmx1out 1e-5 pass
chk63e:
compare_res chk63f "d part of Mx, local coordinates" $mx1doutl $rmx1doutl 1e-5 pass
chk63f:

echo " "
echo "$space ... rotation in the presence of constraining field"
set refout=$testdir/out.$ext.nc3
echo "$space rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log,shfac,vshft}.fe2"
             rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log,shfac,vshft}.fe2
echo "$space cp $testdir/{eula0,b0,ctrl,site2}.fe2 ."
             cp $testdir/{eula0,b0,ctrl,site2}.fe2 .
if (! $?MPIK) then
runrdcmd chk64a %11f $testout "-cat:JOBNC3 --noerr ctrl.$ext"
else
runrdcmd chk64a %11f $testout "-cat:JMPNC3 --noerr ctrl.$ext"
endif
chk64a:
echo -n "$space shortening output ... "
echo "catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
      catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext

call zdiffiles chk64b "CPU 1 $testout $refout"
chk64b:

set ehf3ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk3ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehf3    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk3    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ediffhk = `echo $ehk3 $ehk3ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediffhf = `echo $ehf3 $ehf3ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

set mx3in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
set mx3dout  = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $3}'`
set mx3doutl = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   2' | awk '{print $6}'`
set mx3out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $3}'`
zcat $refout >$testout{}~
set rmx3in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
set rmx3dout  = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $3}'`
set rmx3doutl = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   2' | awk '{print $6}'`
set rmx3out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $3}'`
rm $testout{}~

if (! $?quiet) then

echo "$space input  Mx     = $mx3in"
echo "$space output Mx     = $mx3out"
echo "$space output Mx(d)  = $mx3dout"
echo "$space in loc. coord = $mx3doutl"
echo "$space reference     = $rmx3doutl"

echo ' '
echo "$space ehk           = $ehk3"
echo "$space reference     = $ehk3ref"
echo "$space difference    = $ediffhk"
echo "$space ehk-ehk(coll) = "`echo $ehk3 $ehk0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space ehk-ehk(noB)  = "`echo $ehk3 $ehk1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

echo ' '
echo "$space ehf           = $ehf3"
echo "$space reference     = $ehf3ref"
echo "$space difference    = $ediffhf"
echo "$space ehf-ehf(coll) = "`echo $ehf3 $ehf0 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space ehf-ehf(noB)  = "`echo $ehf3 $ehf1 | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo ' '
endif

compare_res chk64c "ehk" $ehk3 $ehk3ref 1e-6 pass
chk64c:
compare_res chk64d "ehf" $ehf3 $ehf3ref 1e-6 pass
chk64d:
compare_res chk64e "Mx" $mx3out $rmx3out 1e-5 pass
chk64e:
compare_res chk64f "d part of Mx, local coordinates" $mx3doutl $rmx3doutl 1e-5 pass
chk64f:

echo " "
echo "$space ... Longitudinal susceptibility with estat field to remove charge contribution"
set refout=$testdir/out.$ext.nc4
echo "$space rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log,shfac,vshft}.fe2"
             rm -f {a,a2,eula0,b0,ctrl,site2,mixm,log,shfac,vshft}.fe2
echo "$space cp $testdir/{eula0,b0,ctrl,site2}.fe2 ."
             cp $testdir/{eula0,b0,ctrl,site2}.fe2 .
if (! $?MPIK) then
runrdcmd chk65a %11f $testout "-cat:JOBNC4 --noerr ctrl.$ext"
else
runrdcmd chk65a %11f $testout "-cat:JMPNC4 --noerr ctrl.$ext"
endif
chk65a:
call zdiffiles chk65b "CPU 1 $testout $refout"
chk65b:


set ehf3ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk3ref = `zcat $refout | grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehf3    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehf=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ehk3    = ` cat $testout| grep ehf | tail -1 | awk '{match($0,"ehk=[^ ]*"); print substr($0,RSTART+4,RLENGTH-4)}'`
set ediffhk = `echo $ehk3 $ehk3ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
set ediffhf = `echo $ehf3 $ehf3ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`

set m1in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $5}'`
set m2in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout | grep -E '2   \*' | awk '{print $5}'`
set m1out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '1   \*' | awk '{print $5}'`
set m2out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout | grep -E '2   \*' | awk '{print $5}'`

set q1out   = `catf -start:n=3:s=ATOM= -stop:rel:l=1 $testout | vextract . Qv | sed -n 1,1p`
set q2out   = `catf -start:n=3:s=ATOM= -stop:rel:l=1 $testout | vextract . Qv | tail -1`

zcat $refout >$testout{}~
set rm1in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $5}'`
set rm2in    = `catf -start:n=1:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '2   \*' | awk '{print $5}'`
set rm1out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '1   \*' | awk '{print $5}'`
set rm2out   = `catf -start:n=2:s=AMAGNC -stop:rel:l=6 $testout{}~ | grep -E '2   \*' | awk '{print $5}'`
set rq1out   = `catf -start:n=3:s=ATOM= -stop:rel:l=1 $testout{}~ | vextract . Qv | sed -n 1,1p`
set rq2out   = `catf -start:n=3:s=ATOM= -stop:rel:l=1 $testout{}~ | vextract . Qv | tail -1`

rm -f $testout{}~

if (! $?quiet) then

  echo "$space output m1     = $m1out"
  echo "$space reference     = $rm1out"
  set ediff = `echo $m1out $rm1out | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference    = $ediff"
  echo ' '

  echo "$space output m2     = $m2out"
  echo "$space reference     = $rm2out"
  set ediff = `echo $m2out $rm2out | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference    = $ediff"
  echo ' '

  echo "$space output q1     = $q1out"
  echo "$space reference     = $rq1out"
  set ediff = `echo $q1out $rq1out | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference    = $ediff"
  echo ' '

  echo "$space output q2     = $q2out"
  echo "$space reference     = $rq2out"
  set ediff = `echo $q2out $rq2out | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference    = $ediff"
  echo ' '

endif

# compare_res chk65c "ehk" $ehk3 $ehk3ref 1e-6 pass
# chk65c:
# compare_res chk65d "ehf" $ehf3 $ehf3ref 1e-6 pass
# chk65d:
compare_res chk65e "m1" $m1out $rm1out 1e-6 pass
chk65e:
compare_res chk65f "q1" $q1out $rq1out 1e-6 pass
chk65f:

if ($?clean) then
else if ($?pass) then
    echo "$space test 6 PASSED ($ext)"
else
    echo "$space test 6 FAILED ($ext)"
    set failed = ($failed 6)
endif
chk6e:

echo $joblist | grep 7 >/dev/null
if ($status) goto chk7e

set ext = er

cat <<EOF

         --- Test case 7: noncollinear LDA+U hamiltonian ---

EOF

if ($ext == "er") then
cat <<EOF
         Test Er compares a collinear, antiferromagnetic ASA LDA+U
         calculation in hcp Er (one spin up, the other down) to a
         two kinds of equivalent noncollinear calculations:

         1.  Moments are equal and opposite in sign.
             Spins are aligned parallel, but not along z.

         2.  Moments are equal and identical in sign.
             The spinor of the second atom is rotated 180 degrees about y
             to recover the AFM result.
EOF
endif

set pass

set testout=out.lm.{$ext}.coll refout=$testdir/out.lm.$ext.coll
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk7e
endif
echo ' '
query chk71 chk7e 'run this test'
chk71:
# ... Look for executables
findcmd chk71a rdcmd "$path" "$topdir"
chk71a:
findcmd chk71b lm "$path" "$topdir"
chk71b:
findcmd chk71c lmstr "$path" "$topdir"
chk71c:

#  findcmd chk71d rdfile "$path" "no"
#  chk71d:
#  if ("$found" == "no") then
#    echo "$space ... no rdfile in path ... skipping test"
#    goto chk7e
#  endif

echo "$space ... set up ASA strux and starting potential"
touch ctrl.$ext
echo "$space rm -f *.$ext"
             touch ctrl.$ext
             rm -f *.$ext
if ($?clean) then
  echo "$space rm -f out.lm.{$ext}.coll out.lm.{$ext}.nc1 out.lm.{$ext}.nc2 specialspec1 atparms"
               rm -f out.lm.{$ext}.coll out.lm.{$ext}.nc1 out.lm.{$ext}.nc2 specialspec1 atparms
  goto chk7e
endif
echo "$space cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 $testdir/atparms ."
             cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 $testdir/atparms .
echo "$space cp $testdir/dmats.$ext.afm dmats.$ext"
             cp $testdir/dmats.$ext.afm dmats.$ext
echo "$space cp $testdir/rsta.$ext.afm rsta.$ext"
             cp $testdir/rsta.$ext.afm rsta.$ext


echo " "
echo "$space ... AFM calculation, collinear case"
runrdcmd chk72a %11f $testout "-cat:JOBCOLL --noerr ctrl.$ext"
chk72a:
#  echo -n "$space shortening output ... "
#  echo "catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext"
#        catf --x -start:n=2:s=SECMAT -stop:s=BZWTS out.$ext >out.$ext{}~; mv out.$ext{}~ out.$ext
#  call zdiffiles chk72b "CPU 2 $testout $refout"
#  chk72b:
call zdiffiles chk7c9 "CPU -1 $testout $refout"
chk7c9:

compare_resf chk7c1 efc efcref "Fermi" 4 0 ";"
chk7c1:
compare_resf chk7c2 sevc sevcref "sumev=" 4 0 "sumev="
chk7c2:
compare_resf chk7c3 ehkc ehkcref "LM:" 3 0 "ehk="
chk7c3:
compare_resf chk7c4 amomc amomcref "ATOM=" 6 0 mom=
chk7c4:

compare_res chk7c5 "Fermi level" $efc $efcref 1e-6 pass
chk7c5:
compare_res chk7c6 "Mag. moment" $amomc $amomcref 1e-5 pass
chk7c6:
compare_res chk7c7 "sum of evals" $sevc $sevcref 1e-6 pass
chk7c7:
compare_res chk7c8 "ehk" $ehkc $ehkcref 1e-6 pass
chk7c8:


echo " "
echo "$space ... Equivalent AFM calculation, both spins rotated by a common (random) angle (see site.$ext)"
set testout=out.lm.{$ext}.nc1 refout=$testdir/out.lm.$ext.nc1
echo "$space rm -f *.$ext"
             touch ctrl.$ext
             rm -f *.$ext
echo "$space cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 ."
             cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 .
echo "$space cp $testdir/dmats.$ext.afm dmats.$ext"
             cp $testdir/dmats.$ext.afm dmats.$ext
echo "$space cp $testdir/rsta.$ext.afm rsta.$ext"
             cp $testdir/rsta.$ext.afm rsta.$ext

runrdcmd chk7r0 %11f $testout "-cat:JOBNC1 --noerr ctrl.$ext"
chk7r0:
call zdiffiles chk7r9 "CPU -1 $testout $refout"
chk7r9:

compare_resf chk7r1 efr efrref "Fermi" 4 0 ";"
chk7r1:
compare_resf chk7r2 sevr sevrref "sumev=" 4 0 "sumev="
chk7r2:
compare_resf chk7r3 ehkr ehkrref "LM:" 3 0 "ehk="
chk7r3:
compare_resf chk7r4 amomr amomrref "ATOM=" 6 0 mom=
chk7r4:
compare_resf chk7r4a Myr Myrref "<My>=" 2 0 "<My>="
chk7r4a:

compare_res chk7r5 "Fermi level" $efr $efrref 1e-6 pass
chk7r5:
compare_res chk7r6 "Mag. moment" $amomr $amomrref 1e-5 pass
chk7r6:
compare_res chk7r7 "sum of evals" $sevr $sevrref 1e-6 pass
chk7r7:
compare_res chk7r8 "ehk" $ehkr $ehkrref 1e-6 pass
chk7r8:

echo " "
echo "$space ... AFM calculation by noncollinear rotation second spin through 180 degrees"
set testout=out.lm.{$ext}.nc2 refout=$testdir/out.lm.$ext.nc2
echo "$space rm -f *.$ext"
             touch ctrl.$ext
             rm -f *.$ext
echo "$space cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 ."
             cp $testdir/{ctrl,site,syml}.$ext $testdir/specialspec1 .
echo "$space cp $testdir/dmats.$ext.nc2 dmats.$ext"
             cp $testdir/dmats.$ext.nc2 dmats.$ext
echo "$space cp $testdir/rsta.$ext.nc2 rsta.$ext"
             cp $testdir/rsta.$ext.nc2 rsta.$ext
echo "$space cp $testdir/site.$ext.nc2 site.$ext"
             cp $testdir/site.$ext.nc2 site.$ext
if (! $?MPIK) then
runrdcmd chk7n0 %11f $testout "-cat:JOBNC1 --noerr ctrl.$ext"
else
runrdcmd chk7n0 %11f $testout "-cat:JMPNC1 --noerr ctrl.$ext"
endif
chk7n0:
call zdiffiles chk7n9 "CPU -1 $testout $refout"
chk7n9:

compare_resf chk7n1 efn efnref "Fermi" 4 0 ";"
chk7n1:
compare_resf chk7n2 sevn sevnref "sumev=" 4 0 "sumev="
chk7n2:
compare_resf chk7n3 ehkn ehknref "LM:" 3 0 "ehk="
chk7n3:
compare_resf chk7n4 amomn amomnref "ATOM=" 6 0 mom=
chk7n4:
compare_resf chk7n4a Myn Mynref "<My>=" 2 0 "<My>="
chk7n4a:

if (! $?quiet) then
echo ' '
echo "$space ... The following numbers are extracted from the last iteration:"
echo "$space collinear Fermi level            = $efc"
echo "$space rotated collinear Fermi level    = $efr"
echo "$space noncollinear Fermi level         = $efn"
echo "$space reference                        = $efnref"
set diff = `echo $efn $efnref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
echo "$space collinear Er mag. moment         = $amomc"
echo "$space rotated collinear Er mag. moment = $amomr"
echo "$space noncollinear Er mag. moment      = $amomn"
echo "$space reference                        = $amomnref"
set diff = `echo $amomn $amomnref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
echo "$space rotated collinear Er <My>        = $Myr"
echo "$space noncollinear Er <My>             = $Myn"

echo ' '
echo "$space collinear sum of evals           = $sevc"
echo "$space rotated collinear sum of evals   = $sevr"
echo "$space noncollinear sum of evals        = $sevn"
echo "$space reference                        = $sevnref"
set diff = `echo $sevn $sevnref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
echo "$space collinear ehk                    = $ehkc"
echo "$space rotated collinear ehk            = $ehkr"
echo "$space noncollinear ehk                 = $ehkn"
echo "$space reference                        = $ehknref"
set diff = `echo $ehkn $ehknref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
echo "$space difference                       = $diff"

echo ' '
endif

compare_res chk7n5 "Fermi level" $efn $efnref 1e-6 pass
chk7n5:
compare_res chk7n6 "Mag. moment" $amomn $amomnref 1e-5 pass
chk7n6:
compare_res chk7n7 "sum of evals" $sevn $sevnref 1e-6 pass
chk7n7:
compare_res chk7n8 "ehk" $ehkn $ehknref 1e-6 pass
chk7n8:

if ($?clean) then
else if ($?pass) then
    echo "$space test 7 PASSED ($ext)"
else
    echo "$space test 7 FAILED ($ext)"
    set failed = ($failed 7)
endif
chk7e:

echo $joblist | grep 8 >/dev/null
if ($status) goto chk8e

set ext = fe

cat <<EOF

         --- Test case 8: On-site orbital noncollinearity ($testdir/ctrl.$ext)  ---

         The self-consistent magnetic moment and total energy of
         collinear Fe is calculated, with the energy bands.  (The FM
         state is calculated in noncollinear mode so the both spins
         appear in one figure).

         Using the SAME potential and Fermi level, the bands are recalculated
         with the dxy orbital rotated by 180 degrees.
         These bands are identical at certain symmetry points in the BZ,
         e.g. on the line connecting Gamma and H points (0,0,0) and (0,0,1)

         Finally a pass is made for the self-consistent potential generated with
         with the dxy orbital spin flipped, and the bands drawn.

EOF
endif

set pass

set refout=$testdir/out.$ext.meula testout=out.$ext
if (! -e $refout) then
  echo "$space ... skipping test : missing reference file $refout"
  goto chk8e
endif
echo ' '
query chk81 chk8e 'run this test'
chk81:
# ... Look for executables
findcmd chk81a rdcmd "$path" "$topdir"
chk81a:
findcmd chk81b lm "$path" "$topdir"
chk81b:
findcmd chk81c lmstr "$path" "$topdir"
chk81c:
findcmd chk81d vextract "$path" "$topdir"
chk81d:
findcmd chk81e rdfile "$path" "no"
chk81e:
if ("$found" == "no") then
  echo "$space ... no rdfile in path ... skipping test"
  goto chk8e
endif

echo " "
echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f {ctrl,mixm,moms,wkp,site2,log,fe,fe2,eula,bfield,save,sdot,str,sv,shfac}.fe"
             rm -f {ctrl,mixm,moms,wkp,site2,log,fe,fe2,eula,bfield,save,sdot,str,sv,shfac}.fe
if ($?clean) then
  touch bnds.$ext
  echo "$space rm -f bnds*.$ext syml.$ext eulaxy.$ext $testout"
               rm -f bnds*.$ext syml.$ext eulaxy.$ext $testout
  goto chk8e
endif
echo "$space cp $testdir/{ctrl,eulaxy,syml}.$ext ."
             cp $testdir/{ctrl,eulaxy,syml}.$ext .

echo "$space ... generate band passes and three bands files"
if (! $?MPIK) then
runrdcmd chk82a %11f $testout "-cat:TSTLMNC --noerr ctrl.$ext"
else
runrdcmd chk82a %11f $testout "-cat:TSTMPNC --noerr ctrl.$ext"
endif
chk82a:
call zdiffiles chk82b "CPU -1 $testout $refout"
chk82b:

set mmom1    = `cat $testout | $vextract . mom | sed -n 1,1p`
set mmom1ref = `zcat $refout | $vextract . mom | sed -n 1,1p`
set eh1      = `cat $testout | $vextract . ehf | sed -n 1,1p`
set eh1ref   = `zcat $refout | $vextract . ehf | sed -n 1,1p`

set mmom2    = `cat $testout | $vextract . mom | tail -1`
set mmom2ref = `zcat $refout | $vextract . mom | tail -1`
set eh2      = `cat $testout | $vextract . ehf | tail -1`
set eh2ref   = `zcat $refout | $vextract . ehf | tail -1`

cnvt_d_fmt chk83e mmom1 $mmom1
chk83e:
cnvt_d_fmt chk83f mmom1ref $mmom1ref
chk83f:
cnvt_d_fmt chk83g eh1ref $eh1ref
chk83g:
cnvt_d_fmt chk83h eh1 $eh1
chk83h:

cnvt_d_fmt chk84e mmom2 $mmom2
chk84e:
cnvt_d_fmt chk84f mmom2ref $mmom2ref
chk84f:
cnvt_d_fmt chk84g eh2ref $eh2ref
chk84g:
cnvt_d_fmt chk84h eh2 $eh2
chk84h:

if ($?detol8 == 0) set detol8 = 1.0e-6
if ($?dqtol8 == 0) set dqtol8 = 5.0e-6

if (! $?quiet) then
  echo " "
  echo "$space energy of FM state               = $eh1"
  echo "$space energy of reference              = $eh1ref"
  set ediff = `echo $eh1 $eh1ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference    =  $ediff"
  echo ' '
  if ($?mmom1) then
  echo "$space magnetic moment of FM state      = $mmom1"
  echo "$space FM reference magnetic moment     = $mmom1ref"
  set mdiff = `echo $mmom1 $mmom1ref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"
  endif

  echo " "
  echo "$space energy of xy spin-flipped state  = $eh2"
  echo "$space energy of reference              = $eh2ref"
  set ediff = `echo $eh2 $eh2ref  | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space                    difference    =  $ediff"
  echo ' '
  if ($?mmom2) then
  echo "$space magnetic moment of xy S-F state  = $mmom2"
  echo "$space AFM reference magnetic moment    = $mmom2ref"
  set mdiff = `echo $mmom2 $mmom2ref | awk '{{k=($1-$2)>0?($1-$2):($2-$1);} print k}'`
  echo "$space difference                       = $mdiff"

  echo ' '
  set ediff = `echo $eh1 $eh2  | awk '{{k=($1-$2)} print k}'`
  echo "$space AFM-FM energy difference         = $ediff"
  echo ' '
  endif

endif

# pass checks
chk8c:

# ... Check that FM total energy is within tol of reference
compare_res chk8cb "Ferromagnetic ehf" $eh1 $eh1ref $detol8 pass
chk8cb:

if (! $?mmom1) goto chk8cd
compare_res chk8cd "Ferromagnetic mmom" $mmom1 $mmom1ref $dqtol8 pass
chk8cd:

# ... Check that AFM total energy is within tol of reference
compare_res chk8ce "xy-AFM ehf" $eh2 $eh2ref $detol8 pass
chk8ce:

# ... Check that AFM total energy is within tol of reference
compare_res chk8cf "xy-AFM mmom" $mmom2 $mmom2ref $dqtol8 pass
chk8cd:

compare_res chk8cf "xy-AFM ehf" $eh2 $eh2ref $detol8 pass
chk8cf:

set bndstol = 1e-4
zcmpmfiles_res_0 chk8cg "... Max deviation in bnds-coll.$ext from reference" $bndstol pass 4 bnds-coll.$ext $testdir/bnds-coll.$ext
chk8cg:
zcmpmfiles_res_0 chk8ch "... Max deviation in bndsxy0.$ext from reference" $bndstol pass 4 bndsxy0.$ext $testdir/bndsxy0.$ext
chk8ch:
zcmpmfiles_res_0 chk8ci "... Max deviation in bndsxy.$ext from reference" $bndstol pass 4 bndsxy.$ext $testdir/bndsxy.$ext
chk8ci:

if ($?clean) then
else if ($?pass) then
    echo "$space test 8 PASSED ($ext)"
else
    echo "$space test 8 FAILED ($ext)"
    set failed = ($failed 8)
endif
chk8e:

echo ' '
if ($?clean) then
    exit 0
else if ($#failed <= 1) then
    echo "$space $testfile : all noncollinear tests PASSED ($joblist)"
    echo " "
    exit 0
else
    shift failed
    echo "$space $testfile : these tests FAILED:" $failed
    echo " "
    exit -1
endif

# ---------------- sdjob --------------
# runs sd test for input $sdmod; returns to label eojob
exit
sdjob:
if ($?quiet) then
else
endif
set pass
set refout=$testdir/out.fccfe.3k.sdmod=$sdmod testout=out.fccfe
echo ' '
query chk21 chk2e 'run this test'
chk21:
# ... Look for executables
findcmd chk21a rdcmd "$path" "$topdir"
chk21a:
findcmd chk21b lm "$path" "$topdir"
chk21b:
findcmd chk21c lmstr "$path" "$topdir"
chk21c:

echo "$space ... set up ASA strux and starting potential"
echo "$space rm -f *.fccfe"
             touch ctrl.fccfe
             rm -f *.fccfe
if ($?clean) then
  echo "$space rm -f $testout"
               rm -f $testout
  goto $eojob
endif
echo "$space cp $testdir/ctrl.fccfe ."
             cp $testdir/ctrl.fccfe .
echo "$space cp $testdir/$eulafile eula.fccfe"
             cp $testdir/$eulafile eula.fccfe
runjob chk22 $testout "lmstr fccfe"
chk22:
runjob chk23 $testout "lm -vnit=0 -vsdyn=t fccfe --no-iactive"
chk23:
if (! $?MPIK) then
runjob chk24 $testout "lm -vmet=1 -vnk=5 -vnit=15 -vsdmod=$sdmod -vsdyn=t fccfe --no-iactive"
else
runjob chk24 $testout "mpirun -n 8 lm -vmet=3 -vnk=5 -vnit=15 -vsdmod=$sdmod -vsdyn=t fccfe --no-iactive"
endif
chk24:
if ($?add0) then
  echo -n "         ..." ; $add0 $testout
else if ($?poszer) then
  echo -n "         ..." ; $poszer $testout
endif

set dqn   =  `cat $testout | grep 'RMS DQ=' | grep 'file mq' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dmn   =  `cat $testout | grep 'RMS DQ=' | grep 'file ma' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dqnr  =  `zcat $refout | grep 'RMS DQ=' | grep 'file mq' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`
set dmnr  =  `zcat $refout | grep 'RMS DQ=' | grep 'file ma' | tail -1 | awk '{sub(".*RMS DQ=","");sub("last it=.*",""); print $0}'`

#set etest = `grep "LM: it" $testout | awk '{print substr($6,5)}' | tail -1`
#set eref = `zcat $refout | grep "LM: it" | awk '{print substr($6,5)}' | tail -1`
set etest = `grep ' ehf= '  $testout | grep -v last | awk '{print $6}' | tail -1`
set eref = `zcat $refout | grep ' ehf= '  | grep -v last | awk '{print $6}' | tail -1`
#  if ($sdmod >= 1000) then
#    set dq = `grep 'file ma' $testout | tail -1 | awk '{print substr($9,4)}'`
#    set etest = `grep 'ehf=' $testout | tail -1 | awk '{print substr($2,5)}'`
#    set eref = `zcat $refout | grep ehf= | tail -1 | awk '{print substr($2,5)}'`
#  endif

if ($?quiet) then
else
  echo ' '
  echo "$space rms difference output-input charge dq, last iteration = $dqn"
  echo "$space rms change in euler angles dm, last iteration =         $dmn"
  echo "$space ehf,    last iteration =                               $etest"

  echo ' '
  echo "$space Compare mag. forces 1st iteration to $refout"
  cat $testout | awk '{if ($1 == "MAGTRQ:") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p
  echo  ---
  zcat $refout | awk '{if ($1 == "MAGTRQ:") {print;getline;print;getline;print;getline;print;getline;print;getline;print}}' | sed -n 1,6p

  echo ' '
  echo "$space Compare last iterations of SV line to file $refout"
  cat $testout | grep SV: | tail -5
  echo "---"
  zcat $refout | grep SV: | tail -5
chk28a:

  echo ' '
  echo "$space Compare last iterations of av. mag. to file $refout"
  cat $testout | grep amagnc | tail -5
  echo "---"
  zcat $refout | grep amagnc | tail -5
chk28b:

  call showout chk28c CPU
chk28c:

  echo ' '
  call zdiffiles chk28d "CPU 1 $testout $refout"
chk28d:

endif
goto $eojob


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

# ---------------- runrdcmd --------------
exit
runrdcmd:
  set quitjob=$retcall
  if ($outfile == ".") then
    echo "$space Invoking rdcmd will execute the following job(s):"
    $rdcmd -f:$rdcmdfmt --n $callarg
    echo "$space $rdcmd '-f:rdcmd:%2f' $callarg"
                 $rdcmd '-f:rdcmd:%2f' $callarg
    set retval = $status
  else
    if (`echo $outfile | awk '{print substr($1,1,2)}'` == '>>') then
      set appfile = `echo $outfile | awk '{print substr($1,3)}'`
      echo "$space $callarg  >> $appfile"
      exit
#      $callarg >> $appfile
      set retval = $status
    else
      echo "$space ... the following job(s) will be executed by invoking "\""rdcmd $callarg"\"
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
  if ($retval != 0) echo "$space ... $testfile aborting"
  exit $retval

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
  echo -n "$space $keyword ($testvar) within tol ($tol) of reference ($refvar)? ... "
  if (`echo $testvar $refvar | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);tl=1.001*tol} print (k<=tl)}'`) then
    echo yes
  else
    echo no
    unset $passvar
  endif
  goto $quitjob

# ---------------- compare_res_0 --------------
# Compares a number $testvar and unsets $passvar if |testvar|<tol
# usage: compares_res_0 retcall keyword testvar tol passvar
#   keyword      : label (for printout)
#   testvar      : first number
#   tol          : tolerance
#   passvar      : $passvar is unset if |testvar|<tol
exit
compare_res_0:
  set quitjob=$retcall
#  echo $retcall $keyword $testvar $tol $passvar
 echo -n "$space $keyword ($testvar) smaller than tol ($tol)? ... "
  if (`echo $testvar 0 | awk -v tol=$tol '{{k=($1-$2)>0?($1-$2):($2-$1);} print (k<=tol)}'`) then
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
# usage: zcmpnfiles_res_0 retcall keyword testvar tol passvar ndig srcfile reffile
#   retcall      : return to this point in script on exit
#   keyword      : label (for printout)
#   tol          : tolerance in maximum allowed deviation
#   passvar      : $passvar is unset if |testvar|<tol
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

# ---------------- zcmpnfiles --------------
# Compares two files, treating each field as a number.
# call arguments should contain 3 strings: n test-file reference-file
# |n| = number of digits which numbers are truncated to.
# If n<0, sort files before comparing them
# Files with .gz or .Z extensions are assumed to be gzipped.
# Returns with retval = number of differences in reduced files
# Example :  call zcmpnfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $tmpdir/tmp_compnfile_1 $tmpdir/tmp_compnfile_2
exit
zcmpnfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)
  set digits = $zcmpnargs[1]
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
  endif
  if ("$zcmpnargs[3]:e" == 'gz' || "$zcmpnargs[3]:e" == 'Z') then
    set cat2 = 'gunzip -c'
  else
    set cat2 = cat
  endif

  if ($?lsort) then
    $cat1  $zcmpnargs[2] | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk "$a" > $fn1
    $cat2  $zcmpnargs[3] | sed s/D-/E-/g | sed s/D+/E+/g | sort | awk "$a" > $fn2
  else
    $cat1  $zcmpnargs[2] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
    $cat2  $zcmpnargs[3] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2
  endif
  set ncharfile = `wc $fn1 | awk '{print $3}'`
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
# files with .gz or .Z extensions are assumed to be gzipped.
# Returns with retval = max numerical difference
# Example :  call zcmpmfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $tmpdir/tmp_compnfile_1 $tmpdir/tmp_compnfile_2
exit
zcmpmfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)
  set digits = $zcmpnargs[1]
# set a = ' { for (i = NF; i > 0; --i) printf " %.'$digits'f", $i; printf "\n" }'
  set a = ' { for (i = 1; i <= NF; i++) { k = sprintf("%.'$digits'f",$i); if (k+k == 0) k = 0 ; printf "%s ", k}; printf "\n" }'

  set fn1 = $tmpdir/tmp_compnfile_1
  set fn2 = $tmpdir/tmp_compnfile_2
  if ("$zcmpnargs[2]:e" == 'gz' || "$zcmpnargs[2]:e" == 'Z') then
    set cat1 = 'gunzip -c'
  else
    set cat1 = cat
  endif
  if ("$zcmpnargs[3]:e" == 'gz' || "$zcmpnargs[3]:e" == 'Z') then
    set cat2 = 'gunzip -c'
  else
    set cat2 = cat
  endif

  $cat1  $zcmpnargs[2] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
  $cat2  $zcmpnargs[3] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2

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

# ---------------- compare_resf --------------
# Extracts one element of a line in files $testout and $refout containing a keyword.
# Variables testout and refout point to file names and must be set beforehand ($refout is gzipped file)
# usage: compare_resf retcall testvar refvar keyword arg_number occur_number sed_strn
#   Variables testout and refout referring to file names must be set
#   testvar      : put result from file $testout into this variable
#   refvar       : put result from file $refout (compressed) into this variable
#   keyword    	 : string line must contain
#   arg_number 	 : extracts $arg_number'th entry in line, as defined by awk
#   occur_number : argument from $occur_number'th line; if zero, use last line
#   sed_strn     : purge this string from result before assigning
exit
compare_resf:
  set quitjob=$retcall
# echo $retcall $testvar $refvar $keyword $arg_number $occur_number $sed_strn
# grep "$keyword" $testout | awk -v ncnt=0 -v num=$arg_number -v count=$occur_number '{ncnt+=1; if (ncnt==count || count == 0) {print $num}}' | sed "s/$sed_strn//"
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

# ---------------- zcmpnfiles --------------
# Compares two files, treating each field as a number.
# call arguments should contain 3 strings: no-digits test-file reference-file
# Files with .gz or .Z extensions are assumed to be gzipped.
# Returns with retval = number of differences in reduced files
# Example :  call zcmpnfiles chk25 "6 dos-cls.$ext $testdir/dos-cls.$ext.gz"
# Creates temporary files $tmpdir/tmp_compnfile_1 $tmpdir/tmp_compnfile_2
exit
zcmpnfiles:
  set quitjob=$retcall
  set zcmpnargs = ($callarg)
  set digits = $zcmpnargs[1]
# set a = ' { for (i = NF; i > 0; --i) printf " %.'$digits'f", $i; printf "\n" }'
  set a = ' { for (i = 1; i <= NF; i++) { k = sprintf("%.'$digits'f",$i); if (k+k == 0) k = 0 ; printf "%s ", k}; printf "\n" }'

  set fn1 = $tmpdir/tmp_compnfile_1
  set fn2 = $tmpdir/tmp_compnfile_2
  if ("$zcmpnargs[2]:e" == 'gz' || "$zcmpnargs[2]:e" == 'Z') then
    set cat1 = 'gunzip -c'
  else
    set cat1 = cat
  endif
  if ("$zcmpnargs[3]:e" == 'gz' || "$zcmpnargs[3]:e" == 'Z') then
    set cat2 = 'gunzip -c'
  else
    set cat2 = cat
  endif

  $cat1  $zcmpnargs[2] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn1
  $cat2  $zcmpnargs[3] | sed s/D-/E-/g | sed s/D+/E+/g | awk "$a" > $fn2
  set ncharfile = `wc $fn1 | awk '{print $3}'`
  cmp $fn1 $fn2 >/dev/null
  set retval = $status

  if ($retval == 0) rm -f $fn1 $fn2
  if ($retval == 0) goto $quitjob

  set retval = `cmp -l $fn1 $fn2 |& grep -v EOF | wc | awk '{printf "%d", $1}'`
  if ($retval == 0) set retval = '-1'
  rm -f $fn1 $fn2
  goto $quitjob

# ---------------- diffiles --------------
exit
diffiles:
  set quitjob=$retcall
  if ($?quiet) goto $quitjob
  set files = "$callarg"
  query diff11 $quitjob "compare $files"
diff11:
  diff $files | sed -n 1,100p
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
  grep $callarg $testout
  if (`cat $testout | grep $callarg | wc | awk '{print $1}'` > 1) echo ' ---'
  zcat $refout | grep $callarg
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
  Usage: invoke with:   $testfile [switches] job-list

   Job:    tests
     1     (fccfe) static noncollinear calculation: 3k structure in fcc fe, and
                   3k+spin spiral structure, and
                   ++-- structure, through spin spiral cell or directly with 4-atom cell
     2     (fccfe) Spin statics using sdmod=1 for 4 atoms/cell Fe.
     3     (fccfe) Spin statics using sdmod=11 for 4 atoms/cell Fe.
     4     (fccfe) Spin statics using sdmod=1011 for 4 atoms/cell Fe.
     5     (fe)    Application of external B field (Zeeman term)
     6     (fe2)   Limited magnetic exchange interactions by finite rotations
     7     (er)    Noncollinear LDA+U hamiltonian
     8     (fe)    On-site orbital noncollinearity

EOF
exit

# ---------------- usage: --------------
usage:
cat <<EOF
 usage: test.nc [switches] [testcase-list | --all]
        e.g., "test.nc 3 4"
        switches:
        --list       lists the tests you can run
        --quiet      runs tests with minimal output and without prompting user
        --no-iactive runs tests without prompting user
        --all        run through a default list of test cases
        --noplot     skip any steps that generate a plot'
        --clean      clean up files generated by this script
        --add0       add suppressed or leading zeros in output for real numbers \`.nnn'
        --poszer     strips (-) sign from numbers represented as 0
        --MPIK       run lm in MPIK mode
                     Run the MPI job as
                     mpirun -n # lm ...

EOF
exit -1
