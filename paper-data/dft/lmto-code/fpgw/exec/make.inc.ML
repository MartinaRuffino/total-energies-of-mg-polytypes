### You need to set switches (1) to (6), by hand
###

#(1) Compilar ###################
# ... Fortran and linker switches for machine LINUX with intel fortran


FC = ftn  -ffixed-line-length-132
F90 = $(FC)
F90M = $(FC)


#(2) CPP SWITCHES ###################

CPPSWITCH_INTELLINUXIFC  = \
-DEXPAND_ISWAP  -DEXPAND_VDV   -DCOMMONLL  -DEXPAND_MELPLN2         \
-DUSE_X0KBLAS   -DX0KBLAS_DIV -DETIME_DUMMY  \
-DMbytes_X0KBLAS_DIV=2        -DNWORD_RECORDSIZE=1     -DEXPAND_SORTEA \
-DUSE_MPI -DSECONDS


#(3) Compilar options ###################
#
#FFLAGS = $(CPPSWITCH_INTELLINUXIFC) -O  -DUSE_LAPACK -DUSE_MPI -em -J.

FFLAGS = $(CPPSWITCH_INTELLINUXIFC) -DUSE_LAPACK -DUSE_MPI -em -J. -g -fcheck=all -finit-real=snan -Wall -Wtabs -Wno-unused -fbacktrace -ffpe-trap=invalid,zero,overflow -fsignaling-nans -fdump-core




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



#(4) BLAS + LAPACK ############
#

LIBLOC = -lfftw3 /opt/acml/default/gfortran64_int64/lib/libacml.a 


#(5) Linker ####################
#
LK = $(FC) $(FFLAGS)

LKFLAGS2 = $(LIBMATH)
LKFLAGS_OPENMP = $(LIBLOC_OPENMP)

#-Vaxlib -tpp7
#LKFLAGS2 = $(LIBMATH) -Vaxlib -tpp7 -lpthread



#(6) Root of ecal #############
ECAL   = $(HOME)/Fortran/QSGW/trunk/ecalj-master
BINDIR = $(HOME)/bin/code2

