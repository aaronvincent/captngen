FC=gfortran
FOPT= -O3 -fPIC
NUMDIR = ./numerical
QAGDIR = ./numerical/dqag
WDIR = ./Wfunctions
RDIR = ./Rfunctions

MAIN = mainTest.o #mainOper.o
MFSHR = sharedcap.o
MFOBJ = gencap.o
MFCAP = opercap.o
NUMOBJ =  dgamic.o d1mach.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o
WFUNC = WM.o WS2.o WS1.o WP2.o WMP2.o WP1.o WD.o WS1D.o
RFUNC = RM.o RS2.o RS1.o RP2.o RMP2.o RP1.o RD.o RS1D.o


gentest.x: $(MAIN) $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)
	${FC} $(FOPT) -o gentest.x $(MAIN) $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)


gencaplib.so: $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)
	$(FC) $(FOPT) -shared -o $@ $(MFSHR) $(MFOBJ) $(MFCAP) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)


$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFSHR): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFOBJ): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFCAP): %.o: %.f90
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
	rm -f *.o *.so *.mod gentest.x
