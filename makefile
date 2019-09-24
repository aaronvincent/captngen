FC=gfortran
FOPT= -O3
NUMDIR = ./numerical
QAGDIR = ./numerical/dqag
WDIR = ./Wfunctions
RDIR = ./Rfunctions

MAIN = mainOper.o
MFOBJ = gencap.o
NUMOBJ =  dgamic.o d1mach.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o
WFUNC = WM.o WS2.o WS1.o WP2.o WMP2.o WP1.o WD.o WS1D.o
RFUNC = RM.o RS2.o RS1.o RP2.o RMP2.o RP1.o RD.o RS1D.o


gentest.x: $(MAIN) $(MFOBJ) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)
	${FC} -o gentest.x $(MAIN) $(MFOBJ) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)
#	rm $(MFOBJ) $(NUMOBJ) $(QAG)


gencaplib.so: $(MFOBJ) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)
	$(FC) $(FOPT) -shared -o $@ $(MFOBJ) $(NUMOBJ) $(QAG) $(WFUNC) $(RFUNC)


$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<

$(MFOBJ): %.o: %.f90
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
	rm -f *.o *.so gentest.x
