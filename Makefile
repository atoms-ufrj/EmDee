# Define DEBUG or FAST mode:
#   Build the "fast" version with: `make` or `make DEBUG=0`
#   Build the "debug" version with: `make DEBUG=1`
DEBUG?=0

# Define models:
PAIR = lj_cut lj_cut_coul_cut lj_sf lj_sf_coul_sf
BOND = harmonic

# Add prefixes:
PAIRMODELS = $(addprefix pair_,$(PAIR))
BONDMODELS = $(addprefix bond_,$(BOND))

# Compilers and their basic options:
FORT = gfortran
CC   = gcc

BASIC_F_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1
BASIC_C_OPTS = -march=native -m64 -fPIC -fopenmp -cpp -fmax-errors=1

# Warning-related options:
BASIC_F_OPTS += -Wall -Wno-maybe-uninitialized
BASIC_C_OPTS += -Wall -Wno-maybe-uninitialized

# Option FAST (default):
FAST_F_OPTS = -Ofast
FAST_C_OPTS = -Ofast

# Option DEBUG:
DEBUG_F_OPTS = -g -Og -fopenmp -fcheck=all -Ddebug
DEBUG_C_OPTS = -g -Og -fopenmp -fstack-check -fsanitize=null -fbounds-check -Ddebug

# Checks chosen option:
ifeq ($(DEBUG), 1)
  F_OPTS = $(BASIC_F_OPTS) $(DEBUG_F_OPTS)
  C_OPTS = $(BASIC_C_OPTS) $(DEBUG_C_OPTS)
else
  F_OPTS = $(BASIC_F_OPTS) $(FAST_F_OPTS)
  C_OPTS = $(BASIC_C_OPTS) $(FAST_C_OPTS)
endif

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./test
LIBDIR = ./lib
INCDIR = ./include

LIBS = -lgfortran -lm -lgomp
EMDEELIB = -L$(LIBDIR) -lemdee

LN_INC_OPT = -I$(INCDIR)
LN_SO_OPT = -Wl,-rpath,'$$ORIGIN/../lib'
LN_OPTS = $(LN_INC_OPT) $(LN_SO_OPT)

obj = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(1)))
src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))

OBJECTS = $(call obj,EmDeeCode ArBee math structs models \
                     $(PAIRMODELS) pairModelClass \
                     $(BONDMODELS) bondModelClass \
                     modelClass lists global)

.PHONY: all test clean install uninstall lib

.DEFAULT_GOAL := all

all: lib

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)
	rm -rf $(call src,compute_pair compute_bond models)

install:
	cp $(LIBDIR)/libemdee.* /usr/local/lib/
	cp $(INCDIR)/emdee.* /usr/local/include/
	ldconfig

uninstall:
	rm -f /usr/local/lib/libemdee.so
	rm -f /usr/local/include/emdee.h /usr/local/include/emdee.f03
	ldconfig

# Executables:

test: $(addprefix $(BINDIR)/,testfortran testc testjulia)

$(BINDIR)/testfortran: $(SRCDIR)/testfortran.f90 $(INCDIR)/emdee.f03 $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	$(FORT) $(F_OPTS) -o $@ $(LN_OPTS) -J$(OBJDIR) $< $(EMDEELIB)

$(BINDIR)/testc: $(SRCDIR)/testc.c $(INCDIR)/emdee.h $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	$(CC) $(C_OPTS) -o $@ $(LN_OPTS) $< $(EMDEELIB) -lm

$(BINDIR)/testjulia: $(SRCDIR)/testjulia.jl
	mkdir -p $(BINDIR)
	(echo '#!/usr/bin/env julia' && echo 'DIR="${CURDIR}"' && cat $<) > $@
	chmod +x $@

# Static and shared libraries:

lib: $(LIBDIR)/libemdee.so

$(LIBDIR)/libemdee.so: $(OBJECTS)
	mkdir -p $(LIBDIR)
	$(FORT) -shared -fPIC -o $@ $^ $(LIBS)

# Object files:

$(OBJDIR)/EmDeeCode.o: $(call src,EmDeeCode compute_pair compute_bond compute_angle compute_dihedral) \
                       $(call obj,ArBee structs models lists global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/ArBee.o: $(SRCDIR)/ArBee.f90 $(call obj,math global)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/math.o: $(SRCDIR)/math.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/structs.o: $(SRCDIR)/structs.f90 $(OBJDIR)/models.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(SRCDIR)/compute_pair.f90: $(call src,$(addprefix compute_,$(PAIRMODELS)))
	bash $(SRCDIR)/make_compute.sh $(PAIRMODELS) > $@

$(SRCDIR)/compute_bond.f90: $(call src,$(addprefix compute_,$(BONDMODELS)))
	bash $(SRCDIR)/make_compute.sh $(BONDMODELS) > $@

$(OBJDIR)/models.o: $(SRCDIR)/models.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(SRCDIR)/models.f90: $(call obj,$(PAIRMODELS)) $(call obj,$(BONDMODELS))
	bash $(SRCDIR)/make_models_module.sh $(PAIRMODELS) $(BONDMODELS) > $@

$(OBJDIR)/pair_%.o: $(SRCDIR)/pair_%.f90 $(SRCDIR)/compute_pair_%.f90 $(OBJDIR)/pairModelClass.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/pairModelClass.o: $(SRCDIR)/pairModelClass.f90 $(OBJDIR)/modelClass.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/bond_%.o: $(SRCDIR)/bond_%.f90 $(SRCDIR)/compute_bond_%.f90 $(OBJDIR)/bondModelClass.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/bondModelClass.o: $(SRCDIR)/bondModelClass.f90 $(OBJDIR)/modelClass.o
	$(FORT) $(F_OPTS) -Wno-unused-dummy-argument -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/modelClass.o: $(SRCDIR)/modelClass.f90 $(OBJDIR)/global.o
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/lists.o: $(SRCDIR)/lists.f90
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/global.o: $(SRCDIR)/global.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(F_OPTS) -J$(OBJDIR) -c -o $@ $<
