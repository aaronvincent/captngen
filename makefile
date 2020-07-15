FC=gfortran
FOPT= -O3 -fPIC #-Wall -fbounds-check
AUXDIR = ./numerical
QAGDIR = ./numerical/dqag
# TSDIR = ./numerical/TSPACK

MAIN = main.o
MFOBJ = gencap.o
TRGOBJ = alphakappamod.o nonlocalmod.o transgen.o fastevap.o
NUMOBJ =  dgamic.o d1mach.o
AUXOBJ = sgolay.o spline.o pchip.o fftpack5.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o


# TSOBJ = ENDSLP.o SIGS.o SNHCSH.o STORE.o YPCOEF.o YPC1.o YPC1P.o YPC2.o YPC2P.o TSPSI.o \
 INTRVL.o HVAL.o HPVAL.o


gentest.x: $(MAIN) $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG)
	${FC} -o gentest.x $(MFOBJ) $(MAIN) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG)
#	rm $(MFOBJ) $(NUMOBJ) $(QAG)


gencaplib.so: $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG)
	$(FC) $(FOPT) -shared -o $@ $(MFOBJ) $(TRGOBJ) $(NUMOBJ) $(AUXOBJ) $(QAG)

$(NUMOBJ): %.o : $(AUXDIR)/%.f
	$(FC) $(FOPT) -c  $<

$(AUXOBJ): %.o : $(AUXDIR)/%.f90
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
	rm -f *.o *.so gentest.x akmod.mod capmod.mod
