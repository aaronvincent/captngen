FC=gfortran
FOPT= -O3 -fPIC -std=legacy #latter is required if you are running gcc 10 or later #-Wall -fbounds-check
NUMDIR = ./numerical
QAGDIR = ./numerical/dqag
# TSDIR = ./numerical/TSPACK
WDIR = ./Wfunctions
RDIR = ./Rfunctions

MAIN = main.o
MFSHR = sharedcap.o
MFOBJ = gencap.o
MFCAP = opercap.o
TRGOBJ = alphakappamod.o transgen.o fastevap.o
NUMFOBJ =  dgamic.o d1mach.o
NUMF90OBJ = sgolay.o spline.o pchip.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o
WFUNC = WM.o WS2.o WS1.o WP2.o WMP2.o WP1.o WD.o WS1D.o
RFUNC = RM.o RS2.o RS1.o RP2.o RMP2.o RP1.o RD.o RS1D.o


# TSOBJ = ENDSLP.o SIGS.o SNHCSH.o STORE.o YPCOEF.o YPC1.o YPC1P.o YPC2.o YPC2P.o TSPSI.o \
 INTRVL.o HVAL.o HPVAL.o


gentest.x: $(MAIN) $(MFSHR) $(MFOBJ)  $(TRGOBJ) $(NUMFOBJ) $(NUMF90OBJ) $(QAG) $(WFUNC) $(RFUNC)
	${FC} -o gentest.x $(MAIN) $(MFSHR) $(MFOBJ)  $(TRGOBJ) $(NUMFOBJ) $(NUMF90OBJ) $(QAG) $(WFUNC) $(RFUNC)
#	rm $(MFOBJ) $(NUMFOBJ) $(QAG)


gencaplib.so: $(MFSHR) $(MFOBJ)  $(TRGOBJ) $(NUMFOBJ) $(NUMF90OBJ) $(QAG) $(WFUNC) $(RFUNC)
	$(FC) $(FOPT) -shared -o $@ $(MFSHR) $(MFOBJ)  $(TRGOBJ) $(NUMFOBJ) $(NUMF90OBJ) $(QAG) $(WFUNC) $(RFUNC)

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

$(NUMOBJ): %.o: $(NUMDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(QAG): %.o: $(QAGDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(WFUNC): %.o: $(WDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(RFUNC): %.o: $(RDIR)/%.f
	$(FC) $(FOPT) -c  $<


clean:
	rm -f *.o *.mod *.so gentest.x
