.SUFFIXES: # remove all implicit rules
MAKEFLAGS += --check-symlink-times

# Default target
default_target = $(main)

# Provide runnam as command line argument as
#   make bin runnam=<runnam>

# Default runnam in the absence of the command-line arg.
#runnam_def	=	test-007
runnam_def	=
runnam	:=	$(if $(runnam),$(runnam),$(runnam_def))

ifneq ($(strip $(runnam)),)
suffix := -$(runnam)
endif

#Source file of the main program
#main_src	=	ctwmodes-007.f90
main_src	=	ctwmodes.f90

#Name of the executable
main		=	ctwmodes$(suffix)

#--- Complier flags
FC=gfortran
debug_flags	= -g -ffpe-trap=invalid -Wall -Wno-unused-function
opt_flags	=
#opt_flags	= -O
lapack_flags	= -framework Accelerate # use macOS's builtin lapack
FCFLAGS	=	$(lapack_flags) $(opt_flags) $(debug_flags)

#--- Dependencies
dep_mod	=	ctwmodes_pars.mod
dep_obj	=	ctwmodes_pars.o

%.o: %.f90
	$(FC) $(FCFLAGS) -c $?

%.mod: %.f90
	$(FC) $(FCFLAGS) -c $?

$(main): $(main_src) $(dep_mod) $(dep_obj)
	$(FC) $(FCFLAGS) -o $@ $(main_src) $(dep_obj)

# Default target
.DEFAULT_GOAL := $(default_target)

# ===
.PHONY: clean bin run all
clean:
	-rm -f *.o *.mod *.smod

run: bin
	./$(main)

all: run
	julia ./fortbin2nc.jl
