FC=gfortran
FOPT= -O3 -fPIC -std=legacy -fopenmp# -Wall -fbounds-check -g  #legacy is required if you are running gcc 10 or later 
NUMDIR = ./numerical
QAGDIR = ./numerical/dqag
# TSDIR = ./numerical/TSPACK
WDIR = ./Wfunctions
RDIR = ./Rfunctions

MAIN = main.o
CAPTURE = calcCaptures.o
MFSHR = sharedcap.o
MFOBJ = gencap.o
MFCAP = opercap.o
TRGOBJ = alphakappamod.o spergelpressmod.o transgen.o fastevap.o
NUMFOBJ =  dgamic.o d1mach.o
NUMF90OBJ = sgolay.o spline.o pchip.o fftpack5.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o
WFUNC = WM.o WS2.o WS1.o WP2.o WMP2.o WP1.o WD.o WS1D.o
RFUNC = RM.o RS2.o RS1.o RP2.o RMP2.o RP1.o RD.o RS1D.o


# TSOBJ = ENDSLP.o SIGS.o SNHCSH.o STORE.o YPCOEF.o YPC1.o YPC1P.o YPC2.o YPC2P.o TSPSI.o \
 INTRVL.o HVAL.o HPVAL.o


gencaplib.so: $(MFSHR) $(MFOBJ) $(MFCAP) $(TRGOBJ) $(NUMFOBJ) $(NUMF90OBJ) $(QAG) $(WFUNC) $(RFUNC)
	$(FC) $(FOPT) -shared -o $@ $(MFSHR) $(MFOBJ) $(MFCAP) $(TRGOBJ) $(NUMFOBJ) $(NUMF90OBJ) $(QAG) $(WFUNC) $(RFUNC)

# -L tells the linker where to look for shared libraries
# -rpath puts the location of the libraries in the executable so the load can find them at runtime
# -Wl lets us send options to the linker (which are comma seperated)
gentest.x: $(MAIN) gencaplib.so
	${FC} $(FOPT) -L. -Wl,-rpath,. -o gentest.x $(MAIN) gencaplib.so
#	rm $(MFOBJ) $(NUMFOBJ) $(QAG)

calcCaps.x: $(CAPTURE) gencaplib.so
	${FC} $(FOPT) -L. -Wl,-rpath,. -o calcCaps.x $(CAPTURE) gencaplib.so


$(NUMFOBJ): %.o : $(NUMDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(NUMF90OBJ): %.o : $(NUMDIR)/%.f90
	$(FC) $(FOPT) -c  $<

$(TSOBJ): %.o : $(TSDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(MFSHR): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFOBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFCAP): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(TRGOBJ): %.o: %.f90
	$(FC) $(FOPT) -c $<

$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(CAPTURE): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(NUMOBJ): %.o: $(NUMDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(QAG): %.o: $(QAGDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(WFUNC): %.o: $(WDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(RFUNC): %.o: $(RDIR)/%.f
	$(FC) $(FOPT) -c  $<


clean:
	rm -f *.o *.mod *.so gentest.x calcCaps.x
