FORT = gfortran
CC   = gcc
OPTS = -march=native -ffast-math -fstrict-aliasing -Ofast -fPIC -m64 -fopenmp -Wunused -cpp

SRCDIR = ./src
OBJDIR = $(SRCDIR)/obj
BINDIR = ./test
LIBDIR = ./lib
INCDIR = ./include

EMDEELIB = -L$(LIBDIR) -lemdee

obj = $(addprefix $(OBJDIR)/, $(addsuffix .o, $(1)))
src = $(addprefix $(SRCDIR)/, $(addsuffix .f90, $(1)))

OBJECTS = $(call obj,EmDeeCode ArBee math structs models lists global)

.PHONY: all test clean install uninstall lib

all: lib

clean:
	rm -rf $(OBJDIR)
	rm -rf $(LIBDIR)
	rm -rf $(BINDIR)

install:
	cp $(LIBDIR)/libemdee.* /usr/local/lib/
	cp $(INCDIR)/emdee.* /usr/local/include/
	ldconfig

uninstall:
	rm -f /usr/local/lib/libemdee.a /usr/local/lib/libemdee.so
	rm -f /usr/local/include/emdee.h /usr/local/include/emdee.f03
	ldconfig

# Executables:

test: $(addprefix $(BINDIR)/,testfortran testc testjulia)

$(BINDIR)/testfortran: $(SRCDIR)/testfortran.f90 $(INCDIR)/emdee.f03 $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	$(FORT) $(OPTS) -static-libgfortran -o $@ -J$(OBJDIR) $< $(EMDEELIB)

$(BINDIR)/testc: $(SRCDIR)/testc.c $(INCDIR)/emdee.h $(LIBDIR)/libemdee.so
	mkdir -p $(BINDIR)
	$(CC) -static-libgfortran $(OPTS) -o $@ $< $(EMDEELIB) -lm

$(BINDIR)/testjulia: $(SRCDIR)/testjulia.jl
	mkdir -p $(BINDIR)
	cp $< $@ && chmod +x $@

# Static and shared libraries:

lib: $(LIBDIR)/libemdee.so $(LIBDIR)/libemdee.a

$(LIBDIR)/libemdee.so: $(OBJECTS)
	mkdir -p $(LIBDIR)
	$(FORT) -shared -fPIC -o $@ $^ -lgfortran -lm

$(LIBDIR)/libemdee.a: $(OBJECTS)
	mkdir -p $(LIBDIR)
	ar cr $@ $^

# Object files:

$(OBJDIR)/EmDeeCode.o: $(call src,EmDeeCode compute_pair compute_bond compute_angle compute_dihedral) \
                       $(call obj,ArBee structs models lists global)
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/ArBee.o: $(SRCDIR)/ArBee.f90 $(call obj,math global)
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/math.o: $(SRCDIR)/math.f90 $(OBJDIR)/global.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/structs.o: $(SRCDIR)/structs.f90 $(OBJDIR)/models.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/models.o: $(SRCDIR)/models.f90 $(OBJDIR)/global.o
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/lists.o: $(SRCDIR)/lists.f90
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

$(OBJDIR)/global.o: $(SRCDIR)/global.f90
	mkdir -p $(OBJDIR)
	$(FORT) $(OPTS) -J$(OBJDIR) -c -o $@ $<

