# --- Makemakefile: machine dependent compilers and flags ---
#     machine mnemonic     : INTEL_IFORT
#     Use fortran compiler : ifort  (default = ifort)
#     Use C compiler       : gcc  (default = icc)
#

TOPDIR = $(HOME)/lm
MNEMONIC = INTEL_IFORT
# ... C compiler and flags
CC = gcc
CFLAGS =  -DINTEL_IFORT -O3 -funroll-loops -ffast-math -finline -w
#CFLAGS =  -DINTEL_IFORT -g

intelversion = 15

# ldflags = -L${MKLROOT}/lib/intel64 -lfftw3xf_intel -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread  -L$(LIBXCPATH)/lib64 -lxcf90 -lxc 
# -xCORE-AVX2X


# ... Fortran compiler and flags, and linker ... for the INTEL IA32
FC0 = ifort
# FC0 = /h/dmt/dev/ifort-d/ifort-15.0.2
LK0 = ifort -g
FC = $(FC0)
LK = $(LK0)


FFLAGS = -O3 -132 -module ../mods
FFLAGS = -O3 -132 -module ../mods
FFLAGS_LESS = -O2 -132 -module ../mods
FFLAGS_LESS2 = -O1 -132 -module ../mods
FFLAGS_LESS3 = -O0 -132 -module ../mods
FFLAGS_NONE = -g -132 -module ../mods
#FFLAGS =      -g -132 -module ../mods
#FFLAGS_LESS = -g -132 -module ../mods
#FFLAGS_LESS2 = -g -132 -module ../mods
#FFLAGS_LESS3 = -g -132 -module ../mods
#FFLAGS_NONE = -g -132 -module ../mods

# OPENMP = 1
ifdef OPENMP
FC = $(FC0) -openmp
LK = $(LK0) -openmp
endif

# DEBUG = 1
ifdef DEBUG
FFLAGS =      -g -132
FFLAGS_LESS = -g -132
FFLAGS_LESS2 = -g -132
FFLAGS_LESS3 = -g -132
FFLAGS_NONE = -g -132
endif

# Compile with special bounds checking:
# fcd -C -g -O0 -extend-source 132 -check noarg_temp_created,bounds,pointers,uninit -debug full 

# ... Path to, and MPI-specific arguments for the MPI fortran compiler
F90M = mpif90 -O3 -ip
 
# --- Libraries LIBSLA,LIBLOC ---
# LIBXC = /opt/etsf/lib/libxc.a
LIBSLA = $(HOME)/slatsm/slatsm.a
LIBSLA_MPI = $(HOME)/slatsm/slatsm-MPI.a
LIBSLA_OMP = $(HOME)/slatsm/slatsm-OMP.a 
# Intel's FFT generates wrong answers, e.g. test case lmgf, mnn, job5
# LIBLOC = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t -lfftw3xf_intel -lmkl_intel_lp64 -lmkl_sequential   -lmkl_lapack -lmkl_core -lpthread
LIBLOC = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t $(HOME)/lib/libfftw3.a -lmkl_intel_lp64 -lmkl_sequential  -lmkl_lapack -lmkl_core -lpthread
LIBLOC = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t $(LIBXC) $(HOME)/lib/libfftw3.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 
LIBLOC = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t -L$(LIBXCPATH)/lib -lxc $(HOME)/lib/libfftw3.a -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread 
LIBLOC_MPI = $(HOME)/lib/libfftw3.a $(LIBLOC) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 
ifeq ($(intelversion),15)
 LIBLOC = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L$(LIBXCPATH)/lib -lxcf90 -lxc 
 LIBLOC_MPI = $(LIBLOC) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 
 LIBLOC_MPI = $(LIBLOC)
endif
LIBLOC_OMP = -L/opt/intel/Compiler/11.1/056/mkl/lib/em64t $(HOME)/lib/libfftw3.a /opt/intel/Compiler/11.1/056/mkl/lib/em64t/libmkl_lapack95_lp64.a -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lm
LIBES = $(LIBSLA) $(LIBLOC)
LIBES_MPI = $(LIBSLA_MPI) $(LIBLOC_MPI)
LIBES_OMP = $(LIBSLA_OMP) $(LIBLOC_OMP)
ifdef OPENMP
LIBES = $(LIBSLA_OMP) $(LIBLOC_OMP)
endif
 
# --- Make rules ---
.f.o: 
	$(FC) $(FFLAGS) -c $<
 
.c.o: 
	$(CC) $(CFLAGS) -c $<
 
.f.a: 
	$(FC) $(FFLAGS) -c $<
	ar rv $@ $%
	rm -f $%
 
.c.a: 
	$(CC) $(CFLAGS) -c $<
	ar rv $@ $%
	rm -f $%
 
# Shell variables to be supplied by configure:
SHELL = /bin/sh
RANLIB = ranlib
# ... for compilers such as xlf that output numbers 0.nnn as .nnn
ADD0 = 
# ... ccomp related
CCOMPDIR = $(HOME)/bin
CCOMP = $(HOME)/bin/ccomp
INAME = subs-to-intel_ifort
# ... Compiler extensions
FC_IS_F90 = yes
