import matplotlib.pyplot as plot

# global list to reference the order of the isotopes, and their colours
isotopeList = ["H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24","Al27","Si28","S32","Ar40","Ca40","Fe56","Ni59"]
isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{59}Ni$"]
colours = ["#ef1a1a", "#2fef19", "#0055ff", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670","#ef1a1a", "#0055ff", "#2fef19", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670"]


# plotiso takes the isotope-specific capture from the Catena paper (https://arxiv.org/abs/1501.03729) for a
# specified coupling constant (eg c1_0, c5_1, ect...) and isotope. This is then plotted against the 
# isotope-specific capture rate the code captn_oper finds
def plotiso(isotopePick, couplingConstant, title="Isotope Plot Comparison", savefigname=None):

	isoIndex = isotopeList.index(isotopePick)

	filename = "Oper_data/captest_oper_"+couplingConstant+"_alliso-gs98.dat"
	# read in the file's data outputed from captn
	# organized into columns of DM mass (x axis), then 17 of isotope specific capture rates
	file = open(filename,'r')
	currentLine = file.readline()
	Ms=[]
	ISOs=[]
	while currentLine != "":
		theLine = currentLine.split()
		Ms.append(float(theLine[0]))
		ISOs.append(float(theLine[1+isoIndex])) # offset to pickout the chosen isotope to plot
		currentLine = file.readline()
	file.close()


	filename = "Catena_data_updated/Catena_data_"+couplingConstant+"/"+isotopeList[isoIndex]+".dat"
	# read in the file's data from catena paper
	catenaMs = []
	catenaCs = []
	file = open(filename,'r')
	currentLine = file.readline()
	while currentLine != "":
		theLine = currentLine.split()
		if len(theLine) > 1:
			catenaMs.append(float(theLine[0]))
			catenaCs.append(float(theLine[1]))
		currentLine = file.readline()
	file.close()


	plot.clf()
	ax = plot.subplot(111)

	# plotting the Catena and CaptnOper Data
	plot.plot(catenaMs, catenaCs, color=colours[isoIndex], label="Catena "+isotopeLabels[isoIndex], linestyle='--', linewidth=0.3)
	plot.plot(Ms, ISOs, color=colours[isoIndex], label="CaptnOper "+isotopeLabels[isoIndex], marker='.', linestyle='-', linewidth=0.4, markersize=3)
	
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

iso = isotopeList[0]
c = couplingConstants[0]
plotiso(iso, c, title="Compare "+iso+" in coupling constant "+c, savefigname=c+"_"+iso+"_comparison_plot.png")
