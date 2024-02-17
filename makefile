# -------------------------------- Directories ---------------------------------
SRCDIR = src
NUMDIR = numerical
QAGDIR = $(NUMDIR)/dqag
WDIR = Wfunctions
RDIR = Rfunctions
OBJDIR = obj
BINDIR = bin


# ----------------------- Source Files and their Targets -----------------------
# The files must be sorted in their module call order so they compile in order
CAPTNSRCS = $(addprefix $(SRCDIR)/, \
				sharedcap.f90 \
				gencap.f90 \
				opercap.f90 \
				alphakappamod.f90 spergelpressmod.f90 \
				transgen.f90 fastevap.f90 \
			)
MAINSRC = $(SRCDIR)/main.f90
# Grab the f and f90 source files via wildcards
WRSRCS = $(wildcard $(SRCDIR)/$(WDIR)/*.f $(SRCDIR)/$(RDIR)/*.f)
NUMSRCS = $(wildcard $(SRCDIR)/$(NUMDIR)/*.f*)
QAGSRCS = $(wildcard $(SRCDIR)/$(QAGDIR)/*.f)

# Use a string replace to get target directory/filename.o for each source file
CAPTNOBJS = $(CAPTNSRCS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
MAINOBJ = $(MAINSRC:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
WROBJS = $(WRSRCS:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
temp = $(NUMSRCS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
NUMOBJS = $(temp:$(SRCDIR)/%.f=$(OBJDIR)/%.o)
QAGOBJS = $(QAGSRCS:$(SRCDIR)/%.f=$(OBJDIR)/%.o)

# Name of the library and testing executable
CAPTNGEN_LIBNAME = gencap
TESTING_EXE = gentest.x


# ----------------------------- Compiler and Flags -----------------------------
FC=gfortran
#legacy is required if you are running gcc 10 or later
FFLAGS=-fopenmp -fPIC -std=legacy -J $(OBJDIR)
MISMATCH=-Wno-argument-mismatch # add mismatch flag to some compilations
ifeq ($(debug_mode),true)
	FFLAGS+= -g -O0 -Wall -fbounds-check
else
	FFLAGS+= -O3
endif

# -L tells where the linker to look at compile time
# -Wl sends a comma separated list of arguments to the linker
# -rpath tells the exe where to look at runtime (hence the use of the full path)
LDFLAGS=-L $(BINDIR) -Wl,-rpath,"$(realpath $(BINDIR))"
LDLIBS=-l $(CAPTNGEN_LIBNAME)


# ------------------------------- Phony Targets --------------------------------
.PHONY: lib$(CAPTNGEN_LIBNAME).so $(TESTING_EXE) clean nuke

lib$(CAPTNGEN_LIBNAME).so: $(BINDIR)/lib$(CAPTNGEN_LIBNAME).so
$(TESTING_EXE): $(BINDIR)/$(TESTING_EXE)
clean: # clears all objects and modules
	rm -f *.mod $(OBJDIR)/*.mod $(OBJDIR)/*/*.mod $(OBJDIR)/*/*/*.mod
	rm -f *.o $(OBJDIR)/*.o $(OBJDIR)/*/*.o $(OBJDIR)/*/*/*.o
nuke: clean # and also clears the testing executable and library
	rm -f $(BINDIR)/lib$(CAPTNGEN_LIBNAME).so
	rm -f $(BINDIR)/$(TESTING_EXE)


# -------------------------------- Main Targets --------------------------------
# Targets with recipes to put the library and executable in the correct folders
$(BINDIR)/lib$(CAPTNGEN_LIBNAME).so: $(NUMOBJS) $(QAGOBJS) $(CAPTNOBJS) $(WROBJS) | $(BINDIR)
	$(FC) $(FFLAGS) -shared $^ -o $@

$(BINDIR)/$(TESTING_EXE): $(MAINOBJ) lib$(CAPTNGEN_LIBNAME).so | $(BINDIR)
	$(FC) $(FFLAGS) $(LDFLAGS) $< $(LDLIBS) -o $@


# ------------------------------- Object Targets -------------------------------
# Targets with recipes for each object file in directory structure
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(WDIR)/%.o: $(SRCDIR)/$(WDIR)/%.f | $(OBJDIR)/$(WDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(RDIR)/%.o: $(SRCDIR)/$(RDIR)/%.f | $(OBJDIR)/$(RDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Both pchip.f90 and fftpack5.f90 raise a large number of 'argument-mismatch' errors
$(OBJDIR)/$(NUMDIR)/%.o: $(SRCDIR)/$(NUMDIR)/%.f90 | $(OBJDIR)/$(NUMDIR)
	$(FC) $(FFLAGS) $(MISMATCH) -c $< -o $@

$(OBJDIR)/$(NUMDIR)/%.o: $(SRCDIR)/$(NUMDIR)/%.f | $(OBJDIR)/$(NUMDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(QAGDIR)/%.o: $(SRCDIR)/$(QAGDIR)/%.f | $(OBJDIR)/$(QAGDIR)
	$(FC) $(FFLAGS) -c $< -o $@


# ----------------------------- Directory Targets ------------------------------
# Targets with recipes to create the output directories if they don't exist yet
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

