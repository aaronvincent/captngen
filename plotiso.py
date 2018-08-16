import matplotlib.pyplot as plot

def plotiso(filename, isotope, otherfilename, title="Title", savefigname=None):

	isotopes = ["H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24","Al27","Si28","S32","Ar40","Ca40","Fe56","Ni58"]	
	isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{58}Ni$"]	
	
	# read in the file's data
	file = open(filename,'r')
	currentLine = file.readline()
	Xs=[]
	Ys=[]
	while currentLine != "":
		theLine = currentLine.split()
		Xs.append(float(theLine[0]))
		Ys.append(float(theLine[1+isotopes.index(isotope)]))
		currentLine = file.readline()
	file.close()

	# read in the file's data
	file = open("Catena_data/"+otherfilename+"_data.dat",'r')
	currentLine = file.readline()
	otherXs=[]
	otherYs=[]
	while currentLine != "":
		theLine = currentLine.split()
		otherXs.append(float(theLine[0]))
		otherYs.append(float(theLine[1]))
		currentLine = file.readline()
	file.close()

	plot.clf()
	ax = plot.subplot(111)
	plot.plot(Xs, Ys, color="blue", marker='.', linestyle='--', linewidth=0.4, markersize=3, label=isotopeLabels[isotopes.index(isotope)])
	plot.plot(otherXs, otherYs, color="Black",  linestyle='-', linewidth=0.8, label="Catena Total")
	# Shrink current axis by 20%
	box = ax.get_position()
	ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
	# Put a legend to the right of the current axis
	ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	plot.grid()
	plot.yscale('log')
	plot.xscale('log')
	plot.xlabel(r"Dark Matter Mass $[GeV]$")
	plot.ylabel(r"Capture Rate $[\frac{1}{s}]$")
	plot.title(title)
    
	if savefigname != None:
		plot.savefig(savefigname, dpi=500, bbox_inches="tight")
	else:
		plot.show()
	print()


isotopes = ["H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24","Al27","Si28","S32","Ar40","Ca40","Fe56","Ni58"]	
for i in range(len(isotopes)):
	plotiso("captest_oper_c1-0_alliso.dat", isotopes[i], "c1-0", r"$c_{1}^{0}$", "c1-0_plots_"+str(i+1)+"-"+isotopes[i]+".png")


