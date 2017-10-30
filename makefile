FC=gfortran
switch= -O3
AUXDIR = ./aux
QAGDIR = ./aux/dqag

MAIN = main.o
MFOBJ = gencap.o
NUMOBJ =  dgamic.o d1mach.o
QAG=  dsntdqagse.o dqelg.o dqk21.o dqpsrt.o dsntdqk21.o


gentest.x: $(MAIN) $(MFOBJ) $(NUMOBJ) $(QAG)
	${FC} -o gentest.x $(MAIN) $(MFOBJ) $(NUMOBJ) $(QAG)
#	rm $(MFOBJ) $(NUMOBJ) $(QAG)


gencaplib.so: $(MFOBJ) $(NUMOBJ) $(QAG)
	$(FC) $(FOPT) -shared -o $@ $(MFOBJ) $(NUMOBJ) $(QAG)

$(NUMOBJ): %.o : $(AUXDIR)/%.f
	$(FC) $(FOPT) -c $(switch) $<

$(QAG): %.o : $(QAGDIR)/%.f
	$(FC) $(FOPT) -c $(switch) $<


$(MFOBJ): %.o: %.f90
	$(FC) $(FOPT) -c $(switch) $<

$(MAIN): %.o: %.f90
	$(FC) $(FOPT) -c $(switch) $<


#%.o: %.f
#	$(FC) -c $(switch) $<


clean:
	rm -f *.o *.so gentest.x