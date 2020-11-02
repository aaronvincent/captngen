FC=gfortran
FOPT= -O3 -fPIC #-Wall -fbounds-check
AUXDIR = ./numerical
QAGDIR = ./numerical/dqag
# TSDIR = ./numerical/TSPACK

MAIN = main.o
MFOBJ = gencap.o
TRGOBJ = alphakappamod.o transgen.o fastevap.o
NUMOBJ =  dgamic.o d1mach.o sgolay.o spline.o pchip.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o


# TSOBJ = ENDSLP.o SIGS.o SNHCSH.o STORE.o YPCOEF.o YPC1.o YPC1P.o YPC2.o YPC2P.o TSPSI.o \
 INTRVL.o HVAL.o HPVAL.o


gentest.x: $(MAIN) $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(QAG)
	${FC} -o gentest.x $(MAIN) $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(QAG)
#	rm $(MFOBJ) $(NUMOBJ) $(QAG)


gencaplib.so: $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(QAG)
	$(FC) $(FOPT) -shared -o $@ $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(QAG)

$(NUMOBJ): %.o : $(AUXDIR)/%.f $(AUXDIR)/%.f90
	$(FC) $(FOPT) -c  $<

$(TSOBJ): %.o : $(TSDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(QAG): %.o : $(QAGDIR)/%.f
	$(FC) $(FOPT) -c  $<


$(MFOBJ): %.o: %.f90
	$(FC) $(FOPT) -c $<

$(TRGOBJ): %.o: %.f90
	$(FC) $(FOPT) -c $<

$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c  $<




clean:
	rm -f *.o *.mod *.so gentest.x
