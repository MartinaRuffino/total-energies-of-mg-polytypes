### You need to set switches (1) to (6), by hand
###

# (1) ... Fortran and linker switches
# For Intel fortran
FC0 = ifort -132 -module $(HOME)/ecalj-master/fpgw/mods
# FC0 = ifort -132 -check bobunds -traceback !for check array bound
LK0 = $(FC0)

# For gnu fortran
# FC0 = gfortran
# LK0 = gfortran
# OMPF = -fopenmp
# OMPL = -openmp

# ... For Open MP, uncomment OPENMP switch and set OMPF and OMPL appropriate to your compiler
# OPENMP = 1
OMPF = -openmp
OMPL = -openmp
# OMPF = -fopenmp
# OMPL = -fopenmp
FC = $(FC0)
LK = $(LK0)
ifdef OPENMP
FC = $(FC0) $(OMPF)
LK = $(LK0) $(OMPL)
endif

# Steps for compilation with MPI
USE_MPI = 

# Steps for compilation with OPENMP
# 1. Remove the following:
# Maybe also vcoulq and x0k_sym.F?
# rm ../main/{hx0fp0.m,hx0fp0.sc.m,hsfp0.m,hsfp0.sc.m,h_uumatrix.m}.o ../nfpsrc/syscalls2.o ../gwsrc/{ppbafp.fal,ppbafp.fal_oldroutines,mkzxq,mksecp,mkseci,meltpb,wintzsg}.o
# 2. uncomment the line 'OPENMP = 1' above or set environment variable OPENMP to 1, viz:
# setenv OPENMP 1; make -j4 progs_om; unsetenv OPENMP
# make install_om

# -cm is supress all comment.
# -w95 and -w90 is to remove noisy warning related to f90 and f95 recommendations.
# See http://www.intel.com/software/products/compilers/flin/docs/ug/msg_warn.htm

# (2) ... CPP switches
CPPSWITCH_INTELLINUXIFC0 = -DEXPAND_ISWAP -DCOMMONLL -DX0KBLAS_DIV -DMbytes_X0KBLAS_DIV=2 -DEXPAND_SORTEA

#CPPSWITCH_SR8K = \
#-DEXPAND_ISWAP  -DEXPAND_VDV   -UCOMMONLL  -DDUMMY_ETIME  -DEXPAND_MELPLN2 \
#-DUSE_X0KBLAS -DX0KBLAS_DIV  -UEXPAND_SUBSTITUTION -DCOMMENTOUTfor_PARALLEL \
#-DMbytes_X0KBLAS_DIV=1024     -DNWORD_RECORDSIZE=1      -DEXPAND_SORTEA

CPPSWITCH_INTELLINUXIFC  = $(CPPSWITCH_INTELLINUXIFC0)
ifdef OPENMP
CPPSWITCH_INTELLINUXIFC  = -DOPENMP $(CPPSWITCH_INTELLINUXIFC0)
endif

# (3) ... Compiler flags
### INTEL FORTRAN PENTIUM4 LINUX ###
#FFLAGS_c0 = -O0 -Vaxlib  -cpp $(CPPSWITCH_INTELLINUXIFC)
#FFLAGS    = -Vaxlib -tpp7 -cpp $(CPPSWITCH_INTELLINUXIFC)   # for .o
FFLAGS0=-O3 -WB

# Switch with -O3 caused a stop for gfortran in ubuntu9.10 ---I needed
# older verion OR no -Osomething options. (only for sxcf.*.F)
# FFLAGS0=  -O3 -WB $(CPPSWITCH_INTELLINUXIFC) -fomit-frame-pointer -funroll-loops  -ffast-math -ffixed-line-length-132

FFLAGS0= -O3 -WB $(CPPSWITCH_INTELLINUXIFC) -fomit-frame-pointer -funroll-loops
FFLAGS0= -O3 -WB -xHost $(CPPSWITCH_INTELLINUXIFC) -fomit-frame-pointer -funroll-loops
 
# for ellery
# FFLAGS0=       $(CPPSWITCH_INTELLINUXIFC) -fomit-frame-pointer -funroll-loops
# Any intel with debug
ifdef DEBUG
FFLAGS0=  -g  $(CPPSWITCH_INTELLINUXIFC) -O0
# FFLAGS0=  -g  $(CPPSWITCH_INTELLINUXIFC) -fomit-frame-pointer -funroll-loops -heap-arrays
endif

FFLAGS= $(FFLAGS0)

#
### Don't change para_g = .o ... below (or modify it if you know how this work)
#### don't need to read here #####NoteStart
# Some groups of .f sources are compiled into .c*_o files.  (* is 1 to 4).
# The compile options are in FFLAGS_c*. The others are with .o and FFLAGS. See makefile and Search para_g or so.
# ---> It cause a problem if a source file foo.f, which compiled into foo.c*_o contains USE module,
#      because checkmodule does now just support *.o. In such a case, you have to modify checkmodule by yourself.
#      (This note is by takao. Oct.2003)
##################################NoteEnd
para_g = .o     # ppbafp.f  psi2bc1.f psi2bc.f See makefile.
sxcf_g = .o     # sxcf.f
x0kf_g = .o     # x0kf.f
hqpe_g = .o     # hqpe.f
tet5_g = .o


### alpha for older compaq compilar ###
#FFLAGS = -O4 -fast -arch host -tune host -recursive -cpp $(CPPSWITCH_DECALPHA)   # for .o
#FFLAGS  = -K -O4 -fast -arch host -tune host -recursive -cpp $(CPPSWITCH_DECALPHA)   # for .o
#para_g = .o     # ppbafp.f  psi2bc1.f psi2bc.f See makefile.
#sxcf_g = .o     # sxcf.f
#x0kf_g = .o     # x0kf.f
#hqpe_g = .o     # hqpe.f
#
#### alpha compaq compilar to override the compilar bug (rather the latest compiler requires this) #######################
#FFLAGS     = -O3 -fast -arch host -tune host -recursive -cpp $(CPPSWITCH_DECALPHA)   # for .o
#FFLAGS_c1  = -O1 -fast -arch host -tune host -recursive -cpp $(CPPSWITCH_DECALPHA)   # for .c1_o
#FFLAGS_c2  =
#FFLAGS_c3  =
#FFLAGS_c4  = -O4 -fast -arch host -tune host -recursive -cpp $(CPPSWITCH_DECALPHA)   # for .c4_o
#para_g = .o        # ppbafp.f  psi2bc1.f psi2bc.f
#sxcf_g = .c4_o     # sxcf.f
#x0kf_g = .c4_o     # x0kf.f
#hqpe_g = .c1_o     # hqpe.f
#
#
### for sr8k ###
#FFLAGS    = -Oss -loglist -Xpcomp -limit -noparallel -Xparmonitor  -nosave -64  -cpp $(CPPSWITCH_SR8K)
#FFLAGS_c1 = -Oss -loglist -Xpcomp -limit -parallel -Xparmonitor -uinline=2 -nosave -64  -cpp  $(CPPSWITCH_SR8K)
## We devide .f souces to some groups, which are compiled with the same optins to the objects with the same extentions.
#para_g = .c1_o  # ppbafp.f  psi2bc1.f psi2bc.f
#x0kf_g = .c1_o  # x0kf.f
#sxcf_g = .o     # sxcf.f
#hqpe_g = .o     # hqpe.f



#(4) BLAS + LAPACK ############
### ATLAS BLAS
#LIBMATH= $(ECAL)/BLASplusLAPACK/LAPACK_A/lapack_a.a \
#-L$(ECAL)/BLASplusLAPACK/ATLAS/ -llapack  -lcblas -lf77blas -latlas
#
### Goto's BLAS; faster than ATLAS.
#
### See http://www.cs.utexas.edu/users/flame/goto/
#LIBMATH= $(ECAL)/BLASplusLAPACK/LAPACK_A/lapack_a.a \
#$(ECAL)/BLASplusLAPACK/GotoBLAS/xerbla.o  $(ECAL)/BLASplusLAPACK/GotoBLAS/libgoto_p4_512-r0.6.so
#
## This is for multi-threaded version of GotoBlas...
##$(ECAL)/BLASplusLAPACK/GotoBLAS/xerbla.o  $(ECAL)/BLASplusLAPACK/GotoBLAS/libgoto_p4_512p-r0.6.so
#
### alpha
#LIBMATH=-lcxml #for alpha
#
# I had a problem in zgemm in pwmat.
#LIBLOC= /opt/acml4.2.0/gfortran64/lib/libacml.a -lfftw3

# it caused segmentation fault during lmf. (just after BNDFP: started).
#LIBMATH= /opt/acml4.1.0/gfortran64/lib/libacml.a -lfftw3
# LIBMATH= /opt/acml4.2.0/ifort64_mp/lib/libacml_mp.a -lfftw3 -lguide

# ubuntu
#LIBLOC= /usr/lib64/libfftw3.a /usr/lib64/liblapack.so.3gf /usr/lib64/libblas.a

#LIBLOC= /opt/acml4.4.0/gfortran64/lib/libacml.a /usr/lib64/libfftw3.so.3

LIBLOC= /usr/lib/libfftw3.so.3 /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf

# for phpdl1.ph.kcl.ac.uk
LIBLOC = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t /home/ms4/lib/libfftw3.a -lmkl_intel_lp64 -lmkl_sequential   -lmkl_lapack -lmkl_core -lpthread
LIBLOC_OPENMP = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t /opt/intel/Compiler/11.1/056/mkl/lib/em64t/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

# for samuel
LIBLOC = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t -lxc $(HOME)/lib/libfftw3.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 
LIBLOC_OPENMP = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t /opt/intel/Compiler/11.1/056/mkl/lib/em64t/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
ifeq ($(intelversion),15)
LIBLOC = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L$(LIBXCPATH)/lib
LIBLOC_OPENMP =$(LIBLOC)
endif

# For rosalind
# LIBLOC = -L${MKLROOT}/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${LIBXCPATH}/lib -lxcf90 -lxc
# # LIBLOC = -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L${LIBXCPATH}/lib -lxcf90 -lxc
# LIBLOC_OPENMP =$(LIBLOC)

# For ada
# LIBLOC = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L$(LIBXCPATH)/lib -lxcf90 -lxc 
# LIBLOC_OPENMP =$(LIBLOC)

#for ellery.eas.asu.edu
# LIBLOC = -L/opt/intel/Compiler/11.0/074/mkl/lib/em64t -lfftw3xf_intel -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core -lpthread
# LIBLOC = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t /export/svan2/lib/libfftw3.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
# LIBLOC_OPENMP = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t /opt/intel/Compiler/11.1/056/mkl/lib/em64t/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

#for henry.kcl.ac.uk, intel-11.1 compiler
# LIBLOC = -L/share/apps/intel/Compiler/11.1/080/mkl/lib/em64t  -lmkl_intel_lp64 -lmkl_sequential -lmkl_lapack -lmkl_core -lpthread
# LIBLOC_OPENMP = -L/share/apps/intel/Compiler/11.1/080/mkl/lib/em64t /share/apps/intel/Compiler/11.1/080/mkl/lib/em64t/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm

# centos yum install blas, yum install lapack
#LIBLOC= -lfftw3 /usr/lib64/liblapack.so.3.0.3 /usr/lib64/libblas.a

#LIBLOC= -lfftw3   $(HOME)/kit/numericallib/LAPACK/lapack_core2gfortran.a \
# $(HOME)/kit/numericallib/LAPACK/blas_core2gfortran.a \
# $(HOME)/kit/numericallib/LAPACK/tmglib_core2gfortran.a

#LIBLOC= -lfftw3 -L/opt/intel/mkl/10.0.2.018/lib/em64t/lib \
# -lmkl_lapack -lmkl_em64t  -lmkl_core

#for ubuntu thinkpadt61.
#LIBLOC=  /usr/lib64/libfftw3.so.3.1.2 /usr/lib64/liblapack.a

#/usr/lib64/libblas-3.a

#LIBLOC = -L/usr/lib64/atlas/ /usr/lib64/atlas/liblapack.so.3 \
#          /usr/lib64/atlas/libf77blas.so.3 /usr/lib64/atlas/libcblas.so.3 \
#          /usr/lib64/atlas/libatlas.so.3 -lfftw3
#for hpclogin
# LIBLOC= -L/usr/local/intel/2012/composerxe-2011.3.174/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64 -lpthread
# LIBLOC_OPENMP= -L/usr/local/intel/2012/composerxe-2011.3.174/mkl/lib/intel64 -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -L/usr/local/intel/2012/composerxe-2011.3.174/compiler/lib/intel64 -liomp5 -lm

#for oakley.osc.edu
# LIBLOC= -L/usr/local/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack95_lp64 -lpthread
# LIBLOC_OPENMP= -L/usr/local/intel/composer_xe_2011_sp1.6.233/mkl/lib/intel64 -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -liomp5 -lm

LIBMATH= $(LIBLOC)

# yum install atlas --> this did not work... normchk.si gave NaN
#LIBMATH = -L/usr/lib64/atlas/ /usr/lib64/atlas/liblapack.so.3 \
#            /usr/lib64/atlas/libf77blas.so.3 /usr/lib64/atlas/libcblas.so.3 \
#            /usr/lib64/atlas/libatlas.so.3 -lfftw3


#(5) Linker ####################
###for linux 586
#LKFLAGS2 = $(ECAL)/slatsm/slatsm.a  -L/usr/intel/mkl/LIB -lmkl32_lapack -lmkl32_p3  -L/usr/lib/gcc-lib/i586-redhat-linux/2.95.3 -lg2c -lpthread
#
###for sr8k
#LKFLAGS2 = $(COMMON) $(ECAL)/slatsm/slatsm.a  -lblas -llapack -lpl -parallel  -lm

LKFLAGS2 = $(LIBMATH)
LKFLAGS_OPENMP = $(LIBLOC_OPENMP)

#-Vaxlib -tpp7
#LKFLAGS2 = $(LIBMATH) -Vaxlib -tpp7 -lpthread



#(6) Root of ecal #############
ECAL   = $(HOME)/ecal
BINDIR = $(HOME)/bin/code2

