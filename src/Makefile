#----------------------------------------------------------------------
# If you modify any of these options, make sure     
# to recompile the code by typing "make clean; make".      
#                                                                      
#----------------------------------------------------------------------


# Select ODE solver

OPT += -DRKDOPRI5
#OPT += -DROSENBROCK4
#OPT += -DEULER
#OPT += -DRK78
#OPT += -DBS


# Only touch if you know what you are doing
#OPT += -DSTROEMGRENTEST
#OPT += -DFIREHOSEDEBUG
OPT += -DMAKETABLES



#----------------------------------------------------------------------
# Here, select compile environment for the target machine.  
# Follow example below.
#----------------------------------------------------------------------

CC       = gcc               # sets the C-compiler

CPP      = g++

OPTIMIZE = -O2 -Wall -g   # sets optimization and warning flags

OPTIONS  = $(OPTIMIZE) $(OPT) 
#OPTIONS  = $(OPTIMIZE) $(OPT) -Wunused-but-set-variable

EXEC     = STARDUST

OBJS     = main.o utils.o log.o config.o allvars.o \
           functions.o sed.o table_read.o  \
           table_ion.o table_temp.o table_compton.o \
           memory.o interpolation.o rt.o

INCL     = allvars.h prototype.h  Makefile

# Libraries
GSL_INCL = -I${HOME}/local/include
GSL_LIB  = -L${HOME}/local/bin

CONF_INCL= -I${HOME}/local/include
CONF_LIB = -L${HOME}/local/lib

BOOST_INCL = -I${HOME}/local/boost-current

CFLAGS   = $(OPTIONS) $(GSL_INCL) $(BOOST_INCL)
CPPFLAGS = 


LIBS     =  -lconfig $(CONF_LIB) $(CONF_INCL) -lgsl -lgslcblas -lm  $(GSL_LIB) $(GSL_INCL)  $(BOOST_INCL)


# Implicit Rules
.SUFFIXES: .o .c .cpp


.c.o:
	$(CPP) -c $< -o $*.o $(CFLAGS)

.cpp.o:
	$(CPP) -c $< -o $@ $(CPPFLAGS)

$(EXEC): $(OBJS) 
	$(CPP) $(OBJS) $(LIBS)  -o $(EXEC)  

$(OBJS): $(INCL) 


clean:
	rm -f $(OBJS) $(EXEC)

