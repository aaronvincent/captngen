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
	$(FC) -shared -o $@ $(MFOBJ) $(NUMOBJ) $(QAG)

$(NUMOBJ): %.o : $(AUXDIR)/%.f
	$(FC) -c $(switch) $<

$(QAG): %.o : $(QAGDIR)/%.f
	$(FC) -c $(switch) $<


$(MFOBJ): %.o: %.f90
	$(FC) -c $(switch) $<

$(MAIN): %.o: %.f90
	$(FC) -c $(switch) $<


#%.o: %.f
#	$(FC) -c $(switch) $<


clean:
	rm $(MFOBJ) $(NUMOBJ) $(QAG) gentest.x