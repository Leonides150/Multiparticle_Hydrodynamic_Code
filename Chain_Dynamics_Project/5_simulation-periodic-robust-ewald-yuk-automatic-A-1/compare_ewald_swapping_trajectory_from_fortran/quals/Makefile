# Makefile for project hele_shaw_quadrupole
# Sagnik Singha
# Created: Fri Aug 28 15:35:41 CDT 2015
#
#
#
#
FORT=gfortran
OPT=-Os -fopenmp -march=nehalem -mcx16 -msahf -mpopcnt \
	 -msse4.2 --param l1-cache-size=32 \
	 --param l1-cache-line-size=64 \
	 --param l2-cache-size=8192 -mtune=nehalem \
	 -ftree-vectorize
LABLAS=-L/SHAREDLAB/numerical_libs/openblas/lib -lopenblas
PNL=-L/SHAREDLAB/numerical_libs/pnl -lpnl
IPNL=-I/SHAREDLAB/numerical_libs/pnl
SR=./subroutines

MOD=

OBJ=module.o \
swapping_trajectory.o
#main.o
#M=main.f90
SRC=swapping_trajectory.f90 

MODULES=module.f90 

EWALD=zgamma.f90 \
functions.f90 \
variable_making_subs.f90 \
analyt.f90 \
ewald_sum.f90 \
gradient.f90 \

EWALD_OBJ=zgamma.o \
functions.o \
variable_making_subs.o \
analyt.o \
ewald_sum.o \
gradient.o \

all:
	$(FORT) $(IPNL) $(OPT) -c $(MODULES) $(SRC)
	$(FORT) $(IPNL) $(OPT) -c $(EWALD) 
	$(FORT) -fopenmp $(LABLAS) $(PNL) $(OBJ) $(EWALD_OBJ) -o main.run

clean:
	rm main.run $(MOD) $(OBJ)
