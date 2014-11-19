#                      Makefile for ham2d

#F90     = ifort
F90      = gfortran
#CC      = icc
CC       = gcc
#FFLAGS  =  -w -O2 -r8 
FFLAGS   =  -w -O2  -lgfortran -fdefault-real-8 #-m64 -mmacosx-version-min=10.10# uncomment this for gfortran
ACCFLAGS =# -acc -ta=nvidia -Minfo=accel 
CFLAGS   =  -w  -O2 $(ACCFLAGS) #-fbacktrace -fbounds-check  #-g #-traceback  -check uninit  #-warn-all#-g #-O2 #-g
LFLAGS   = $(CFLAGS) -lm
DIR      = ./src
OBJECTS  = $(DIR)/ham2d.o \
           $(DIR)/readGrid.o \
           $(DIR)/preprocess.o \
           $(DIR)/find_faces.o\
	        $(DIR)/initflow.o \
	        $(DIR)/stepSolution.o \
	        $(DIR)/updateSoln.o \
	        $(DIR)/computeRHS.o\
           $(DIR)/computeRHSk.o \
           $(DIR)/computeRHSkv.o \
           $(DIR)/computeLinearRHS.o \
           $(DIR)/flux_roe2d.o \
	        $(DIR)/muscld.o \
	        $(DIR)/muscld_deriv.o \
	        $(DIR)/roeflx.o \
	        $(DIR)/computeForce.o \
	        $(DIR)/outputSolution.o \
	        $(DIR)/smoothGrid.o \
	        $(DIR)/jac_roe.o \
	        $(DIR)/wallFlux.o \
	        $(DIR)/flux_visc.o \
	        $(DIR)/triSolvers.o \
	        $(DIR)/mathops.o \
	        $(DIR)/findDiagonals.o\
	        $(DIR)/ADI.o \
	        $(DIR)/DADI.o \
	        $(DIR)/jac_visc.o \
	        $(DIR)/gaussSeidel.o \
	        $(DIR)/lineGaussSeidel.o \
	        $(DIR)/apply_periodic.o \
	        $(DIR)/weno.o \
	        $(DIR)/weno3.o \
	        $(DIR)/weno_deriv.o \
	        $(DIR)/periodic_bc.o \
	        $(DIR)/apply_periodic_LHS.o \

 
# Link Instruction
ham2d: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o bin/ham2d

clean:
	@rm -rf $(DIR)/*.o $(DIR)/*.mod *.*~ bin/ham2d

%.o:%.F90
	$(F90) $(FFLAGS) -c $< -o $*.o

%.o:%.f90
	$(F90) $(FFLAGS) -c $< -o $*.o
%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $*.o

%.o:%.mod
