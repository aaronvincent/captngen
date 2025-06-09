# Default target is the library
CAPTNGEN_LIBNAME = gencap
lib$(CAPTNGEN_LIBNAME).so:
# -------------------------------- Directories ---------------------------------
SRCDIR = src
NUMDIR = numerical
QAGDIR = dqag
WDIR = Wfunctions
RDIR = Rfunctions
OBJDIR_name = obj
BINDIR_name = bin
LIBDIR_name = lib
DEBUG_name = debug

# Put debug-compiled files in a seperate directory
ifneq ($(strip $(debug)),)
	OBJDIR = $(OBJDIR_name)-$(DEBUG_name)
	BINDIR = $(BINDIR_name)-$(DEBUG_name)
	LIBDIR = $(LIBDIR_name)-$(DEBUG_name)
else
	OBJDIR = $(OBJDIR_name)
	BINDIR = $(BINDIR_name)
	LIBDIR = $(LIBDIR_name)
endif


# ----------------------- Source Files and their Targets -----------------------
# The file module call dependencies are defined here:
# gencap.f90 and opercap.f90 use sharedcap.f90
$(addprefix $(OBJDIR)/, \
	gencap.o \
	opercap.o \
): $(addprefix $(OBJDIR)/, \
	sharedcap.o \
)

# spergelpressmod.f90 and fastevap.f90 use gencap.f90
$(addprefix $(OBJDIR)/, \
	spergelpressmod.o \
	fastevap.o \
): $(addprefix $(OBJDIR)/, \
	gencap.o \
)

# transgen.f90 uses gencap.f90, spergelpressmod.f90, and alphakappamod.f90
$(addprefix $(OBJDIR)/, \
	transgen.o \
): $(addprefix $(OBJDIR)/, \
	gencap.o \
	spergelpressmod.o \
	alphakappamod.o \
)

# Grab the f and f90 source files via wildcards
CAPTNSRCS = $(wildcard $(SRCDIR)/*.f90)
MAINSRC = $(SRCDIR)/main.f90
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

# Name of the testing executable
TESTING_EXE = gentest.x


# ----------------------------- Compiler and Flags -----------------------------
FC=gfortran
#legacy is required if you are running gcc 10 or later due to the arguement-mismatch warning being promoted to error
FFLAGS=-fopenmp -fPIC -std=legacy -J $(OBJDIR)
ifneq ($(strip $(debug)),)
	FFLAGS+= -g -O0 -Wall -Wextra -Wconversion
	FFLAGS+= -fbacktrace -fcheck=all -ffpe-trap=zero,overflow,underflow,denormal
	FFLAGS+= -fdebug-aux-vars# -fimplicit-none --- numerical/dgamic.f misbehaves with the implicit-none restriction!
else
	FFLAGS+= -O3
endif

# -L tells where the linker to look at compile time
# -Wl sends a comma separated list of arguments to the linker
# -rpath tells the exe where to look at runtime (hence the use of the full path)
LDFLAGS=-L $(LIBDIR) -Wl,-rpath,"$(realpath $(LIBDIR))"
LDLIBS=-l $(CAPTNGEN_LIBNAME)


# ------------------------------- Phony Targets --------------------------------
# debug=true isn't an actual target, but PHONY makes it appear in tab completion
.PHONY: lib$(CAPTNGEN_LIBNAME).so $(TESTING_EXE) clean nuke debug=true

lib$(CAPTNGEN_LIBNAME).so: $(LIBDIR)/lib$(CAPTNGEN_LIBNAME).so
$(TESTING_EXE): $(BINDIR)/$(TESTING_EXE)
clean: # clears all objects and modules
	rm -f $(SRCDIR)/*.mod $(OBJDIR)/*.mod $(OBJDIR)/*/*.mod
	rm -f $(OBJDIR)/*.o $(OBJDIR)/*/*.o
nuke: clean # and also clears the testing executable and library
	rm -f $(LIBDIR)/lib$(CAPTNGEN_LIBNAME).so
	rm -f $(BINDIR)/$(TESTING_EXE)


# -------------------------------- Main Targets --------------------------------
# Targets with recipes to put the library and executable in the correct folders
$(LIBDIR)/lib$(CAPTNGEN_LIBNAME).so: $(NUMOBJS) $(QAGOBJS) $(CAPTNOBJS) $(WROBJS) | $(LIBDIR)
	$(FC) $(FFLAGS) -shared $^ -o $@

$(BINDIR)/$(TESTING_EXE): $(MAINOBJ) $(LIBDIR)/lib$(CAPTNGEN_LIBNAME).so | $(BINDIR)
	$(FC) $(FFLAGS) $(LDFLAGS) $< $(LDLIBS) -o $@


# ------------------------------- Object Targets -------------------------------
# Targets with recipes for each object file in directory structure
$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(WDIR)/%.o: $(SRCDIR)/$(WDIR)/%.f | $(OBJDIR)/$(WDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(RDIR)/%.o: $(SRCDIR)/$(RDIR)/%.f | $(OBJDIR)/$(RDIR)
	$(FC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/$(NUMDIR)/%.o: $(SRCDIR)/$(NUMDIR)/%.f* | $(OBJDIR)/$(NUMDIR)
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

$(LIBDIR):
	mkdir -p $(LIBDIR)
