import matplotlib.pyplot as plot

# global list to reference the order of the isotopes, and their colours
isotopeList = ["H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24","Al27","Si28","S32","Ar40","Ca40","Fe56","Ni59"]
isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{59}Ni$"]
colours = ["#ef1a1a", "#2fef19", "#0055ff", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670","#ef1a1a", "#0055ff", "#2fef19", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670"]


# plotcatena takes the capture rates from the Catena paper (https://arxiv.org/abs/1501.03729) for a
# specified coupling constant (eg c1-0, c5-1, ect...). This is designed to mimic the plots found in the paper
def plotcatena(couplingConstant, title="Catena Plot", savefigname=None):
	
	
	folder = "Catena_data_updated/Catena_data_"+couplingConstant+"/"

	# read in the file's data from catena paper
	catenaMs = []
	catenaCs = []
	for i in range(len(isotopeList)+1):
		if i == len(isotopeList):
			filename = folder+"Total_"+couplingConstant+".dat"
		else:
			filename = folder+isotopeList[i]+".dat"
		file = open(filename,'r')
		currentLine = file.readline()
		Ms=[]
		Cs=[]
		while currentLine != "":
			theLine = currentLine.split()
			if len(theLine) > 1:
				Ms.append(float(theLine[0]))
				Cs.append(float(theLine[1]))
			currentLine = file.readline()
		catenaMs.append(Ms)
		catenaCs.append(Cs)
		file.close()

	plot.clf()
	ax = plot.subplot(111)
	for i in range(len(isotopeList)+1):
		if i == len(isotopeList):
			plot.plot(catenaMs[i], catenaCs[i], color="Black", linestyle='-', linewidth=0.8, label="Total")
		elif 7<i<len(isotopeList):
			plot.plot(catenaMs[i], catenaCs[i], color=colours[i], label=isotopeLabels[i], linestyle='--', linewidth=0.4)
		else:
			plot.plot(catenaMs[i], catenaCs[i], color=colours[i], label=isotopeLabels[i], linestyle='-', linewidth=0.4)
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

couplingConstants = ["c1-0", "c3-0", "c4-0", "c5-0", "c6-0", "c7-0", "c8-0", "c9-0", "c10-0", "c11-0", "c12-0", "c13-0", "c14-0", "c15-0"]

for c in couplingConstants:
	plotcatena(c, "Copy of Catena Plot of "+c, "Catena_Plot_Copies/"+c+"_catena_copy.png")
