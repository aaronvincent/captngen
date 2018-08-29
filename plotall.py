import matplotlib.pyplot as plot

def plotall(filename, otherfilename, title="Title", savefigname=None):
	
	isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{58}Ni$"]
	colours = ["#ef1a1a", "#2fef19", "#0055ff", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670","#ef1a1a", "#0055ff", "#2fef19", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670"]
	
	# read in the file's data
	file = open(filename,'r')
	currentLine = file.readline()
	Ms=[]
	ISOs=[]
	while currentLine != "":
		theLine = currentLine.split()
		Ms.append(float(theLine[0]))
		tempList = []
		for i in range(1,17):
			tempList.append(float(theLine[i]))
		ISOs.append(tempList)
		currentLine = file.readline()
	file.close()

	totalCap = []
	for currentMass in ISOs:
		totalCap.append(sum(currentMass))

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

	# make each sub list a list of only the isotope changing with DM mass
	ISOs = [list(temp) for temp in zip(*ISOs)]

	plot.clf()
	ax = plot.subplot(111)
	plot.plot(Ms, totalCap, label="Total", color="Black", marker='.', linestyle='--', linewidth=0.8, markersize=3)
	plot.plot(otherXs, otherYs, color="Black", linestyle='-', linewidth=0.8, label="Catena Total")
	for i in range(0,16):
		if i<8:
			plot.plot(Ms, ISOs[i], color=colours[i],  label=isotopeLabels[i], marker='.', linestyle='-', linewidth=0.4, markersize=3)
		else:
			plot.plot(Ms, ISOs[i], color=colours[i], label=isotopeLabels[i], marker='.', linestyle='--', linewidth=0.4, markersize=3)
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

plotall("captest_oper_c3-0_testR.dat", "c3-0", r"$c_{3}^{0}$", "c3-0_plots_alliso_testR.png")

