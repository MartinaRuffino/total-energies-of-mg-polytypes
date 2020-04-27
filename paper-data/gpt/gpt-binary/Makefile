################################################################################
#F90 = gfortran #mpif90
F90 = ifort
# Optimized Flags
#FLAGS = -O3 -march=native -Wall -Wfatal-errors -fbounds-check -finit-real=nan -fimplicit-none -lblas -llapack
# Debug Flags
#FLAGS = -Og -Wall -fbacktrace -fbounds-check -finit-real=nan -fimplicit-none #-lblas -llapack
################################################################################
# IF INTEL COMPILERS AND MKL USE THESE FLAGS:
FLAGS = -O3 -xHost -no-wrap-margin -mkl
#FLAGS = -O0 -g -mkl -traceback -check all -check bounds -check uninit -ftrapuv -debug all
# NO NEED FOR -L, USE -mkl INSTEAD OF -llapack ETC
################################################################################

MDGPT = ~/bin/md.gpt

#LIB = /usr/local/opt/lapack/lib

all: $(MDGPT)

.PHONY: all

$(MDGPT): main.f90
	$(F90) $(FLAGS) -o $@ $<

clean:
	$(RM) $(MDGPT)
