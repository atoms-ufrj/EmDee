# Define DEBUG or FAST mode:
#   Build the "fast" version with: `make` or `make DEBUG=0`
#   Build the "debug" version with: `make DEBUG=1`
DEBUG ?= 0

# Installation prefix:
PREFIX ?= /usr/local

# Compilers and their basic options:
FORT ?= gfortran
CC ?= gcc

BASIC_F_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1
BASIC_C_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1

# Warning-related options:
BASIC_F_OPTS += -Wall -Wno-maybe-uninitialized
BASIC_C_OPTS += -Wall -Wno-maybe-uninitialized

# Option FAST (default):
FAST_F_OPTS = -Ofast
FAST_C_OPTS = -Ofast

# Option DEBUG:
DEBUG_F_OPTS = --coverage -g -Og -fcheck=all -Ddebug
DEBUG_C_OPTS = --coverage -g -Og -fstack-check -fsanitize=null -fbounds-check -Ddebug

# Checks chosen option:
ifeq ($(DEBUG), 1)
  F_OPTS = $(BASIC_F_OPTS) $(DEBUG_F_OPTS)
  C_OPTS = $(BASIC_C_OPTS) $(DEBUG_C_OPTS)
else
  F_OPTS = $(BASIC_F_OPTS) $(FAST_F_OPTS)
  C_OPTS = $(BASIC_C_OPTS) $(FAST_C_OPTS)
endif

LN_INC_OPT = -I$(INCDIR)
LN_SO_OPT = -Wl,-rpath,'$$ORIGIN/../lib'
LN_OPTS = $(LN_INC_OPT) $(LN_SO_OPT)

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./bin
LIBDIR = ./lib
INCDIR = ./include
TSTDIR = ./test

LIBS = -lgfortran -lm -lgomp -lnfft3 -lfftw3
EMDEELIB = -L$(LIBDIR) -lemdee

obj = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(1)))
src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))

# Force-field models (from source files starting with pair_, coul_, bond_, angle_, and dihedral_):
MODELTERMS  = pair coul bond angle dihedral
PAIRMODELS  = $(patsubst $(SRCDIR)/%.f90,%,$(wildcard $(SRCDIR)/pair_*.f90))
COULMODELS  = $(patsubst $(SRCDIR)/%.f90,%,$(wildcard $(SRCDIR)/coul_*.f90))
BONDMODELS  = $(patsubst $(SRCDIR)/%.f90,%,$(wildcard $(SRCDIR)/bond_*.f90))
ANGLEMODELS = $(patsubst $(SRCDIR)/%.f90,%,$(wildcard $(SRCDIR)/angle_*.f90))
DIHEDMODELS = $(patsubst $(SRCDIR)/%.f90,%,$(wildcard $(SRCDIR)/dihedral_*.f90))

ALLMODELS   = $(shell bash $(SRCDIR)/make_pair_list.sh $(PAIRMODELS)) \
              $(COULMODELS) $(BONDMODELS) $(ANGLEMODELS) $(DIHEDMODELS)

KSPACE = $(patsubst $(SRCDIR)/%.f90,%,$(wildcard $(SRCDIR)/kspace_*.f90))

OBJECTS = $(call obj,EmDeeCode EmDeeData ArBee math structs models \
                     $(ALLMODELS) $(KSPACE) $(addprefix modelClass_,$(MODELTERMS) kspace) \
                     modelClass lists nfft global)

COMPUTES = $(addprefix compute_,$(MODELTERMS))
ENERGYCOMPUTES = $(addprefix energy_compute_,pair coul)
VIRIALCOMPUTES = $(addprefix virial_compute_,pair coul)

TESTS = $(patsubst %.f90,%,$(wildcard $(TSTDIR)/*.f90))

.PHONY: all test clean install uninstall lib include

.DEFAULT_GOAL := all

# Phony targets:

all: lib include $(TSTDIR)/testfortran
# DELETE $(TSTDIR)/testfortran above

test: $(addprefix $(BINDIR)/,testc testjulia) $(TESTS)
	cd $(TSTDIR) && bash run_tests.sh

clean:
	rm -rf $(OBJDIR) $(LIBDIR) $(BINDIR) $(INCDIR)
	rm -rf $(call src,$(COMPUTES) $(ENERGYCOMPUTES) $(VIRIALCOMPUTES) models)
	rm -rf $(TESTS) *.gcda *.gcno

install:
	cp $(LIBDIR)/libemdee.* $(PREFIX)/lib/
	cp $(INCDIR)/*.* $(PREFIX)/include/
	ldconfig

uninstall:
	rm -f $(PREFIX)/lib/libemdee.so
	rm -f $(addprefix $(PREFIX)/include/,emdee.h emdee.f03 libemdee.jl)
	ldconfig

lib: $(LIBDIR)/libemdee.so

include: $(INCDIR)/emdee.f03 $(INCDIR)/emdee.h $(INCDIR)/libemdee.jl

# Executables:

$(TSTDIR)/%: $(TSTDIR)/%.f90 $(INCDIR)/emdee.f03 $(TSTDIR)/common/contained.f90 $(OBJDIR)/mConfig.o
	$(FORT) $(F_OPTS) -o $@ $(LN_OPTS) -J$(OBJDIR) $< $(OBJDIR)/mConfig.o $(EMDEELIB)

$(BINDIR)/testc: $(SRCDIR)/testc.c $(INCDIR)/emdee.h $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	$(CC) $(C_OPTS) -o $@ $(LN_OPTS) $< $(EMDEELIB) -lm

$(BINDIR)/testjulia: $(SRCDIR)/testjulia.jl $(INCDIR)/libemdee.jl $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	(echo '#!/usr/bin/env julia' && echo 'DIR="${CURDIR}"' && cat $<) > $@
	chmod +x $@

$(OBJDIR)/mConfig.o: $(TSTDIR)/common/mConfig.f90 $(LIBDIR)/libemdee.so
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

# Shared library and includes:

$(LIBDIR)/libemdee.so: $(OBJECTS)
	mkdir -p $(LIBDIR)
	$(FORT) $(F_OPTS) -shared -fPIC -o $@ $(OBJECTS) $(LIBS)

$(INCDIR)/emdee.f03: $(SRCDIR)/emdee_header.f03 $(SRCDIR)/make_f_header.sh
	mkdir -p $(INCDIR) $(LIBDIR)
	bash $(SRCDIR)/make_f_header.sh $(ALLMODELS) $(KSPACE) > $(INCDIR)/emdee.f03

$(INCDIR)/emdee.h: $(SRCDIR)/emdee_header.h $(SRCDIR)/make_c_header.sh
	mkdir -p $(INCDIR) $(LIBDIR)
	bash $(SRCDIR)/make_c_header.sh $(ALLMODELS) $(KSPACE) > $(INCDIR)/emdee.h

$(INCDIR)/libemdee.jl: $(SRCDIR)/emdee_header.jl $(SRCDIR)/make_j_header.sh
	mkdir -p $(INCDIR) $(LIBDIR)
	bash $(SRCDIR)/make_j_header.sh $(ALLMODELS) $(KSPACE) > $(INCDIR)/libemdee.jl

# Object files:

$(OBJDIR)/EmDeeCode.o: $(call src,EmDeeCode inner_loop) \
                       $(call src,$(COMPUTES) $(ENERGYCOMPUTES) $(VIRIALCOMPUTES)) \
                       $(call obj,EmDeeData ArBee structs models lists global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/EmDeeData.o: $(call src,EmDeeData inner_loop) \
                       $(call obj,ArBee structs models lists math global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/ArBee.o: $(SRCDIR)/ArBee.f90 $(call obj,math global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/structs.o: $(SRCDIR)/structs.f90 $(OBJDIR)/models.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(SRCDIR)/compute_pair.f90: $(call src,$(PAIRMODELS)) $(SRCDIR)/make_compute.sh
	bash $(SRCDIR)/make_compute.sh pair $(PAIRMODELS) > $@

$(SRCDIR)/compute_coul.f90: $(call src,$(COULMODELS)) $(SRCDIR)/make_compute.sh
	bash $(SRCDIR)/make_compute.sh coul $(COULMODELS) > $@

$(SRCDIR)/compute_bond.f90: $(call src,$(BONDMODELS)) $(SRCDIR)/make_compute.sh
	bash $(SRCDIR)/make_compute.sh bond $(BONDMODELS) > $@

$(SRCDIR)/compute_angle.f90: $(call src,$(ANGLEMODELS)) $(SRCDIR)/make_compute.sh
	bash $(SRCDIR)/make_compute.sh angle $(ANGLEMODELS) > $@

$(SRCDIR)/compute_dihedral.f90: $(call src,$(DIHEDMODELS)) $(SRCDIR)/make_compute.sh
	bash $(SRCDIR)/make_compute.sh dihedral $(DIHEDMODELS) > $@

$(SRCDIR)/energy_compute_pair.f90: $(call src,$(PAIRMODELS)) $(SRCDIR)/make_energy_compute.sh
	bash $(SRCDIR)/make_energy_compute.sh pair $(PAIRMODELS) > $@

$(SRCDIR)/energy_compute_coul.f90: $(call src,$(COULMODELS)) $(SRCDIR)/make_energy_compute.sh
	bash $(SRCDIR)/make_energy_compute.sh coul $(COULMODELS) > $@

$(SRCDIR)/virial_compute_pair.f90: $(call src,$(PAIRMODELS)) $(SRCDIR)/make_virial_compute.sh
	bash $(SRCDIR)/make_virial_compute.sh pair $(PAIRMODELS) > $@

$(SRCDIR)/virial_compute_coul.f90: $(call src,$(COULMODELS)) $(SRCDIR)/make_virial_compute.sh
	bash $(SRCDIR)/make_virial_compute.sh coul $(COULMODELS) > $@

$(OBJDIR)/models.o: $(call obj,$(ALLMODELS) $(addprefix modelClass_,$(MODELTERMS))) \
                    $(call obj,$(KSPACE) modelClass_kspace) $(SRCDIR)/make_models_module.sh
	bash $(SRCDIR)/make_models_module.sh $(ALLMODELS) $(KSPACE) > $(SRCDIR)/models.f90
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $(SRCDIR)/models.f90

$(OBJDIR)/pair_%.o: $(SRCDIR)/pair_%.f90 $(OBJDIR)/modelClass_pair.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/coul_%.o: $(SRCDIR)/coul_%.f90 $(OBJDIR)/modelClass_coul.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/bond_%.o: $(SRCDIR)/bond_%.f90 $(OBJDIR)/modelClass_bond.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/angle_%.o: $(SRCDIR)/angle_%.f90 $(OBJDIR)/modelClass_angle.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/dihedral_%.o: $(SRCDIR)/dihedral_%.f90 $(OBJDIR)/modelClass_dihedral.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/kspace_%.o: $(SRCDIR)/kspace_%.f90  $(OBJDIR)/modelClass_kspace.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/modelClass_%.o: $(SRCDIR)/modelClass_%.f90 $(OBJDIR)/modelClass.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/modelClass.o: $(SRCDIR)/modelClass.f90 $(call obj,lists math nfft)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/nfft.o: $(SRCDIR)/nfft.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/lists.o: $(SRCDIR)/lists.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/math.o: $(SRCDIR)/math.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/global.o: $(SRCDIR)/global.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<
