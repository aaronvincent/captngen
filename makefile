SRCDIR = src
NUMDIR = numerical
QAGDIR = $(NUMDIR)/dqag
WDIR = Wfunctions
RDIR = Rfunctions
OBJDIR = obj
BINDIR = bin

# The files must be sorted in their module call order so that they compile in order
CAPTNSRCS = $(addprefix $(SRCDIR)/, \
				sharedcap.f90 \
				gencap.f90 \
				opercap.f90 \
				alphakappamod.f90 spergelpressmod.f90 \
				transgen.f90 fastevap.f90 \
			)
MAINSRC = $(SRCDIR)/main.f90
# Grab the source files via wildcards
WRSRCS = $(wildcard $(SRCDIR)/$(WDIR)/*.f $(SRCDIR)/$(RDIR)/*.f)
NUMSRCS = $(wildcard $(SRCDIR)/$(NUMDIR)/*.f*)
QAGSRCS = $(wildcard $(SRCDIR)/$(QAGDIR)/*.f)

# Use a string replace to get target names and directories for each object file
CAPTNOBJS = $(CAPTNSRCS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
MAINOBJ = $(MAINSRC:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
WROBJS = $(WRSRCS:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
temp = $(NUMSRCS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
NUMOBJS = $(temp:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
QAGOBJS = $(QAGSRCS:$(SRCDIR)/%.f=$(OBJDIR)/%.o)

CAPTNGEN_LIBNAME = gencap
TESTING_EXE = gentest.x

FC=gfortran
#legacy is required if you are running gcc 10 or later
FFLAGS=-fPIC -std=legacy -J $(OBJDIR)
ifeq ($(debug_mode),true)
	FFLAGS+= -g -O0 -Wall -fbounds-check
else
	FFLAGS+= -O3 -fopenmp
endif

LINKER=$(FC)
# The testing executable needs the rpath set by the linker such that at runtime it can find the library inside the binary folder  
LDFLAGS=-fopenmp -L $(BINDIR) -I $(OBJDIR) -Wl,-rpath,"$(realpath $(BINDIR))"
#If the library follows the lib[name].so naming convention, then -l[name] can be used instead of -l:[name]lib.so
LSLIBS=-l$(CAPTNGEN_LIBNAME)


# Phony targets to rename the functional binary_directory/file.out targets.
lib$(CAPTNGEN_LIBNAME).so: $(BINDIR)/lib$(CAPTNGEN_LIBNAME).so
$(TESTING_EXE): $(BINDIR)/$(TESTING_EXE)


# Targets to put the library and test executable in the binary folder
$(BINDIR)/lib$(CAPTNGEN_LIBNAME).so: $(NUMOBJS) $(QAGOBJS) $(CAPTNOBJS) $(WROBJS) | $(BINDIR)
	$(FC) -shared $^ -o $@

$(BINDIR)/$(TESTING_EXE): $(MAINOBJ) lib$(CAPTNGEN_LIBNAME).so | $(BINDIR)
	${LINKER} $(LDFLAGS) $< $(LSLIBS) -o $@


# Targets for each object file
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(WDIR)/%.o: $(SRCDIR)/$(WDIR)/%.f | $(OBJDIR)/$(WDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(RDIR)/%.o: $(SRCDIR)/$(RDIR)/%.f | $(OBJDIR)/$(RDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Both pchip.f90 and fftpack5.f90 raise a large number of 'argument-mismatch' errors
$(OBJDIR)/$(NUMDIR)/%.o: $(SRCDIR)/$(NUMDIR)/%.f90 | $(OBJDIR)/$(NUMDIR)
	$(FC) $(FFLAGS) -Wno-argument-mismatch -c  $< -o $@

$(OBJDIR)/$(NUMDIR)/%.o: $(SRCDIR)/$(NUMDIR)/%.f | $(OBJDIR)/$(NUMDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(QAGDIR)/%.o: $(SRCDIR)/$(QAGDIR)/%.f | $(OBJDIR)/$(QAGDIR)
	$(FC) $(FFLAGS) -c $< -o $@


# Targets to inform the makefile how to create the directories if they don't exist yet
$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/$(WDIR):
	mkdir -p $(OBJDIR)/$(WDIR)

$(OBJDIR)/$(RDIR):
	mkdir -p $(OBJDIR)/$(RDIR)

$(OBJDIR)/$(NUMDIR):
	mkdir -p $(OBJDIR)/$(NUMDIR)

$(OBJDIR)/$(QAGDIR):
	mkdir -p $(OBJDIR)/$(QAGDIR)

$(BINDIR):
	mkdir -p $(BINDIR)


# clean clears all objects and modules
clean:
	rm -f *.mod $(OBJDIR)/*.mod $(OBJDIR)/*/*.mod $(OBJDIR)/*/*/*.mod
	rm -f *.o $(OBJDIR)/*.o $(OBJDIR)/*/*.o $(OBJDIR)/*/*/*.o

# nuke invokes clean and also clears the testing executable and library
nuke: clean
	rm -f $(BINDIR)/lib$(CAPTNGEN_LIBNAME).so
	rm -f $(BINDIR)/$(TESTING_EXE)

.PHONY: clean nuke lib$(CAPTNGEN_LIBNAME).so $(TESTING_EXE)
