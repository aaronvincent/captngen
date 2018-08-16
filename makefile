FC=gfortran
FOPT= -O3
AUXDIR = ./aux
QAGDIR = ./aux/dqag

MAIN = main3.o
MFOBJ = gencap.o
NUMOBJ =  dgamic.o d1mach.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o

WAddress = ./Wfunctions#./Ws_sun/Ws

gentest.x: $(MAIN) $(MFOBJ) $(NUMOBJ) $(QAG) $(WAddress)/WM.f $(WAddress)/WS2.f $(WAddress)/WS1.f $(WAddress)/WP2.f $(WAddress)/WMP2.f $(WAddress)/WP1.f $(WAddress)/WD.f $(WAddress)/WS1D.f
	${FC} -o gentest.x $(MAIN) $(MFOBJ) $(NUMOBJ) $(QAG) $(WAddress)/WM.f $(WAddress)/WS2.f $(WAddress)/WS1.f $(WAddress)/WP2.f $(WAddress)/WMP2.f $(WAddress)/WP1.f $(WAddress)/WD.f $(WAddress)/WS1D.f
#	rm $(MFOBJ) $(NUMOBJ) $(QAG)


gencaplib.so: $(MFOBJ) $(NUMOBJ) $(QAG)
	$(FC) $(FOPT) -shared -o $@ $(MFOBJ) $(NUMOBJ) $(QAG)

$(NUMOBJ): %.o : $(AUXDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(QAG): %.o : $(QAGDIR)/%.f
	$(FC) $(FOPT) -c  $<


$(MFOBJ): %.o: %.f90
	$(FC) $(FOPT) -c $<

$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<




clean:
	rm -f *.o *.so gentest.x
