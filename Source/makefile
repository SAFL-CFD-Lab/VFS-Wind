##################################################################
# Copyright (C) by Regents of the University of Minnesota.       #
#                                                                #
# This Software is released under GNU General Public License 2.0 #
# http://www.gnu.org/licenses/gpl-2.0.html                       #
#                                                                #
##################################################################

###################
# Makefile for VWiS
###################
#
# Do not modify this file
#
# Set the following variables in makefile.local if the default is not applicable
#  ACML ......... path to your acml installation
#  PETSC ........ path to your PETSC installation
#  HYPRE ........ path to your HYPRE installation
#  TEC360HOME ... path to your tecplot installation
#  MPICXX ....... MPI c++ compiler
#  USE_TECPLOT .. set this variable if tecplot is used

sinclude makefile.local

# Set default values
ACML ?= /safl/software/x86_64/acml/4.4.0/gfortran64
PETSC ?= /mnt/lustre2/fsotiropoulos/cald0075/system/petsc-3.1
HYPRE ?= /mnt/lustre2/fsotiropoulos/cald0075/system/petsc-3.1
TEC360HOME ?= /usr
MPICXX ?= mpicxx

# Petsc library and include paths
PETSCLIB=$(PETSC)/lib
PETSCINC=$(PETSC)/include
PETSC_FOUND=$(shell for line in '\#include "petsc.h"' '\#include "petscversion.h"' 'main(){}' ; do echo $$line ; done | $(MPICXX) -xc++ -I $(PETSCINC) -L $(PETSCLIB) -lpetsc -o /dev/null - 2> /dev/null && echo 1 || echo -1)
PETSC_VERSION=$(shell for line in '\#include <stdio.h>' '\#include "petscversion.h"' 'main(){printf("%d.%d", PETSC_VERSION_MAJOR, PETSC_VERSION_MINOR);}' ; do echo $$line ; done | $(MPICXX) -xc++ -I $(PETSCINC) -o petsc_version - && ./petsc_version && rm -f petsc_version)

# Tecplot library and include paths
TECINC=$(TEC360HOME)/include
TECLIB=$(TEC360HOME)/lib
TECPLOT_FOUND=$(shell for line in '\#include "TECIO.h"' 'main(){}' ; do echo $$line ; done | $(MPICXX) -xc++ -I $(TECINC) -L $(TECLIB) -L./ -ltecio64 -o /dev/null - 2> /dev/null && echo 1 || echo -1)

# ACML library paths
ACMLLIB=$(ACML)/lib
ACML_FOUND=$(shell for line in 'main(){}' ; do echo $$line ; done | $(MPICXX) -xc++ -L $(ACMLLIB) -lacml -o /dev/null - 2> /dev/null && echo 1 || echo -1)

# Hypre library and include paths
HYPRELIB = $(HYPRE)/lib
HYPREINC = $(HYPRE)/include
HYPRE_FOUND=$(shell for line in '\#include "HYPRE.h"' 'main(){}' ; do echo $$line ; done | $(MPICXX) -xc++ -I $(HYPREINC) -L $(HYPRELIB) -lHYPRE -o /dev/null - 2> /dev/null && echo 1 || echo -1)

# If any library is set to /usr/lib and we are on a 64 bit system, then change that to /usr/lib64
libdir = /usr/lib
ifeq ($(shell uname -m), x86_64)
	libdir = /usr/lib64
ifeq ($(PETSCLIB), /usr/lib)
	PETSCLIB=$(libdir)
endif
ifeq ($(TECLIB), /usr/lib)
	TECLIB=$(libdir)
endif
ifeq ($(HYPRELIB), /usr/lib)
	HYPRELIB=$(libdir)
endif
ifeq ($(ACMLLIB), /usr/lib)
	ACMLLIB=$(libdir)
endif
endif

#########

LIBDIR=-L$(ACMLLIB) -L$(PETSCLIB) -L$(TECLIB) -L$(HYPRELIB) -L./

LIBFLAG	=	-lpthread -lrt -ldl -lstdc++ -lX11 \
		-lpetsc -lHYPRE -lacml -lgfortran -g

SOURCEC	=	bcs.c bmv.c compgeom.c ibm.c ibm_io.c init.c \
		main.c metrics.c poisson.c rhs.c timeadvancing.c \
		timeadvancing1.c variables.c fsi.c implicitsolver.c\
		fsi_move.c solvers.c rhs2.c wallfunction.c \
		les.c k-omega.c distance.c level.c momentum.c poisson_hypre.c\
                wave.c rotor_model.c

OBJSC	=	$(SOURCEC:%.c=%.o)

ALLSRC = $(wildcard *.c)

CPPFLAGS =	-DNDEBUG -I$(PETSCINC) -I$(TECINC) -I$(HYPREINC) -O3 -MMD -MP

ifdef USE_TECPLOT
	CPPFLAGS += -DTECIO
	LIBTECPLOT = -ltecio64
endif

# Compilation directive for .c files
%.o: %.c
	$(MPICXX) -c -o $@ $(CPPFLAGS) $<

all: vwis

tests:
ifeq ($(PETSC_FOUND), -1)
	$(error PETSC not found. Please set your PETSC variable correctly. It is currently set to: "$(PETSC)")
endif
ifneq ($(PETSC_VERSION), 3.1)
	$(error PETSC version not supported. Your version is $(PETSC_VERSION) whereas it should be 3.1)
endif
ifeq ($(ACML_FOUND), -1)
	$(error ACML not found. Please set your ACML variable correctly. It is currently set to: "$(ACML)")
endif
ifeq ($(HYPRE_FOUND), -1)
	$(error HYPRE not found. Please set your HYPRE variable correctly. It is currently set to: "$(HYPRE)")
endif

tecplottest:
ifeq ($(TECPLOT_FOUND), -1)
ifdef USE_TECPLOT
	$(error TECPLOT not found. Please set your TEC360HOME variable correctly or do not set USE_TECPLOT. The variable is currently set to: "$(TEC360HOME)")
endif
endif

# Main executable
vwis: tests $(OBJSC)
	$(MPICXX) -o vwis $(OBJSC) $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

data: tests tecplottest data.o
	$(MPICXX) -o data data.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) $(LIBTECPLOT)

data05: tests tecplottest data05.o
	$(MPICXX) -o data05 data05.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) $(LIBTECPLOT)

itfcsearch: tests itfcsearch.o
	$(MPICXX) -o itfcsearch itfcsearch.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG)

shear: tests tecplottest shear.o
	$(MPICXX) -o csashear shear.o $(CPPFLAGS) $(LIBDIR) $(LIBFLAG) $(LIBTECPLOT)

clean:
	rm -f *.o *\~  Momt* Force* surface* Distance*.dat* *_0.dat test1

.PHONY: clean

-include $(ALLSRC:%.c=%.d)
