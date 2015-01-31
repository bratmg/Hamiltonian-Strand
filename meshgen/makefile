#####################################################################
#                                                                   #
#                       MAKEFILE FOR MESHGEN                        # 
#                                                                   #
#####################################################################

CC      = gcc
CFLAGS  =  -w -g -pg #-O2# -traceback #  -fbounds-check -fbacktrace  #-g #-traceback  -check uninit  #-warn-all#-g #-O2 #-g
LFLAGS  = $(CFLAGS) -lm
SRCDIR = ./src
OBJECTS = $(SRCDIR)/globalVariables.o \
          $(SRCDIR)/loopsAndCells.o \
          $(SRCDIR)/IO.o \
          $(SRCDIR)/listOperations.o \
          $(SRCDIR)/colouringAlgorithm.o \
          $(SRCDIR)/edgeOperations.o \
          $(SRCDIR)/flushToSurface.o \
          $(SRCDIR)/smoothOperations.o \
          $(SRCDIR)/strand.o \
          $(SRCDIR)/memoryOperations.o \
          $(SRCDIR)/meshQuality.o \
          $(SRCDIR)/subdivision.o \
          $(SRCDIR)/driverMeshGen.o

# Link Instruction
meshgen: $(OBJECTS)
	$(CC) $(OBJECTS) $(LFLAGS) -o bin/meshgen

clean:
	@rm  -f $(SRCDIR)/*.o $(SRCDIR)/*.mod bin/meshgen

%.o:%.c
	$(CC) $(CFLAGS) -c $< -o $*.o

%.o:%.mod

#####################################################################
#                           END OF FILE                             #
#####################################################################
