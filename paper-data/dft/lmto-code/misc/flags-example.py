import os

flags = dict(
    cc      = 'gcc',
    cflags  = '-O3 -march=corei7-avx -pipe',
    fc      = 'gfortran',
    fflags  = '-O3 -march=corei7-avx -pipe -ffixed-line-length-132',

    # The flag to set destination path for .mod files.
    modflag = '-J',

    # The flag enabling OpenMP.
    omp     = '-fopenmp',

    # The nvidia/cuda compiler. Not necessary unless you want to play with it.
    nvcc    = 'nvcc -O3 -arch=sm_35 --ftz true --prec-div false --prec-sqrt false',

    # All general linking related flags (library and run paths, libs etc...)
    ldflags = '-L{MKLROOT}/lib/intel64 -lfftw3xf_gnu -lmkl_gf_lp64 -lmkl_sequential -lmkl_core {MKLROOT}/../compiler/lib/intel64/libiomp5.a -lpthread'.format(MKLROOT=os.environ['MKLROOT']),

    # The C++ compiler to be used if necessary for glue code.
    cxx     = 'g++',

    # Preprocessing flags for the custom ccomp fortran preprocessor. Used for all .f and .F* files except for the qsgw package. The C preprocessor called from the compiler script is used there.
    ppflags = '-dunix -dAUTO-ARRAY -dF90 -dNOQUAD -dBLAS -dBLAS3 -dLAPACK -dFFTW -dFFTW3 -dINTEL_IFORT -dFORTRAN2003',

    # Kernigham & Ritchie archaic C compiler used only to compile ccomp. Beware the HECToR/cray architecture is messed up! The executables produced by this compiler NEED to run on the compilation node. This is usually not the case if the default C compiler is used.
    knrcc   = 'gcc',

    # Flags for the qsgw package (compiler flags as well as preprocessing definitions).
    gwflags = '-ffixed-line-length-132 -DEXPAND_ISWAP -DEXPAND_VDV -DCOMMONLL -DUSE_X0KBLAS -DX0KBLAS_DIV -DMbytes_X0KBLAS_DIV=2 -DNWORD_RECORDSIZE=1 -DEXPAND_SORTEA',

# Flags for routines with excessive branching which can choke a compiler but are not crucial to overall performance.
    special = {
            'v7input/m_rdctrl.f'    : '-O0 -pipe',
            'v7input/m_rdctrlchk.f' : '-O0 -pipe',
            'subs/strxops.f': '-O1 -march=corei7-avx -pipe',
            'subs/spcgrp.f' : '-O0 -march=corei7-avx -pipe'
        }
)


#print flags
