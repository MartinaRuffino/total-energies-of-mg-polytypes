# env FC=ifort CC=icc ./configure
# --- Makemakefile: machine dependent compilers and flags ---
#     machine mnemonic     : INTEL_IFORT
#     Use fortran compiler : ifort  (default = ifort)
#     Use C compiler       : icc  (default = icc)
#
MNEMONIC = INTEL_IFORT
# ... C compiler and flags
CC = icc
CFLAGS =  -DINTEL_IFORT -O
#CFLAGS =  -DINTEL_IFORT -g
# ... Fortran compiler and flags, and linker ... for the INTEL IA32
FC = ifort
LK = ifort
FFLAGS = -O3 -xHost -132
FFLAGS_LESS0 = -O3 -132
FFLAGS_LESS = -O2 -132
FFLAGS_LESS2 = -O1 -132
FFLAGS_LESS3 = -O0 -132
FFLAGS_NONE = -g -132

FFLAGS =      -g -132 -gdwarf-3 -traceback
FFLAGS_LESS = -g -132 -gdwarf-3 -traceback
FFLAGS_LESS2 = -g -132 -gdwarf-3 -traceback
FFLAGS_LESS3 = -g -132 -gdwarf-3 -traceback
FFLAGS_NONE = -g -132 -gdwarf-3 -traceback
# --- Source directory.  Substitute second line if relative path desired ---
srcdir = /h/ms4/slatsm
#srcdir = .
 
# --- Libraries LIBSLA,LIBLOC ---
LIBSLA = $(srcdir)/slatsm.a
LIBSLA_MPI = /h/ms4/slatsm/slatsm-MPI.a
LIBSLA_OPENMP = @LIBSLA_OPENMP@
LIBLOC = -L/opt/intel/composer_xe_2015.5.223/mkl -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -L/opt/libxc/3.0.0/intel/15.0.4/o/lib -lxcf90 -lxc -llapack -lblas
LIBLOC_MPI = -L/opt/intel/composer_xe_2015.5.223/mkl/lib/intel64 -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread  -L/opt/libxc/3.0.0/intel/15.0.4/o/lib -lxcf90 -lxc
LIBLOC_OPENMP = @LIBLOC_OPENMP@
LIBES = $(LIBSLA) $(LIBLOC)
LIBES_MPI = $(LIBSLA_MPI) $(LIBLOC_MPI)
LIBES_OPENMP = $(LIBSLA) $(LIBLOC_OPENMP)
 
# MPIF90 = mpif90 -I. -module -O3 -debug minimal -march=native -extend-source 132 -I/opt/intel/compilers_and_libraries_2018.0.128/linux/mkl/include -I/opt/libxc/4.2.3/intel/18.0.0/o/include -I/opt/hdf5/1.10.4/openmpi/3.1.3/intel/18.0.0/o/include -allow nofpp_comments -c

# --- Make rules ---
.f.o: 
	$(FC) $(FFLAGS) -c $<
 
.c.o: 
	$(CC) $(CMFLAGS) -c $<
 
.f.a: 
	$(FC) $(FFLAGS) -c $<
	ar rv $@ $%
	rm -f $%
 
.c.a: 
	$(CC) $(CMFLAGS) -c $<
	ar rv $@ $%
	rm -f $%
 
# --- Shell variables to be supplied by configure ---
SHELL = /bin/sh
RANLIB = ranlib
# ... ccomp related
CCOMP = $(srcdir)/../bin/ccomp
CCOMPDIR = $(srcdir)/../bin
CCOMP_SW =  -dINTEL_IFORT -dF90 -dFFTW -dIN_PLACE
# ... for compilers such as xlf that output numbers 0.nnn as .nnn
ADD0 = 
POSZER = 
# ... Compiler extensions
FC_IS_F90 = yes
FC_AUTOARRAY = 
FC_POINTER = 
# ... Fortran-C linkage
FC_UPPERCASE = 0
FC_UNDERSCORE = 1
NOPASSEDARGS = 1
CMAIN = MAIN__

CMFLAGS = $(CFLAGS) -DFC_UNDERSCORE=1 -DFC_UPPERCASE=0 -DNOPASSEDARGS=1 -DCMAIN=MAIN__ -DNARGFCALL=nargs_ -DADD_TO_NARGFCALL=0 -DARGFCALL=getarg_
