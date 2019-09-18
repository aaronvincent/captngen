import matplotlib.pyplot as plot

# global list to reference the order of the isotopes
isotopeList = ["H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24","Al27","Si28","S32","Ar40","Ca40","Fe56","Ni59"]

# plotiso takes the isotope-specific capture from the Catena paper (https://arxiv.org/abs/1501.03729) for a
# specified coupling constant (eg c1_0, c5_1, ect...) adn isotope. This is then plotted against the 
# isotope-specific capture rate the code captn_oper finds
def plotiso(filename, isotopePick, catenafilename, title="Title", savefigname=None):

	isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{59}Ni$"]
	colours = ["#ef1a1a", "#2fef19", "#0055ff", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670","#ef1a1a", "#0055ff", "#2fef19", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670"]
	
	isoIndex = isotopeList.index(isotopePick)

	# read in the file's data outputed from captn
	# organized into columns of DM mass (x axis), then 17 of isotope specific capture rates
	file = open(filename,'r')
	currentLine = file.readline()
	dmMasses=[]
	capRates=[]
	while currentLine != "":
		theLine = currentLine.split()
		dmMasses.append(float(theLine[0]))
		capRates.append(float(theLine[1+isoIndex])) # offset to pickout the chosen isotope to plot
		currentLine = file.readline()
	file.close()

	# read in the file's data from catena paper (isotope-specific capture rate data)
	file = open("Catena_data/Catena_data_"+catenafilename+"/"+isotopePick+".dat",'r')
	currentLine = file.readline()
	catenaDMms=[]
	catenaCaps=[]
	while currentLine != "":
		theLine = currentLine.split()
		catenaDMms.append(float(theLine[0]))
		catenaCaps.append(float(theLine[1]))
		currentLine = file.readline()
	file.close()

	plot.clf()
	ax = plot.subplot(111)
	plot.plot(dmMasses, capRates, color=colours[isoIndex], marker='.', linestyle='--', linewidth=0.4, markersize=3, label=isotopeLabels[isoIndex])
	plot.plot(catenaDMms, catenaCaps, color=colours[isoIndex],  linestyle='-', linewidth=0.8, label="Catena "+isotopeLabels[isoIndex])
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


for iso in isotopes:
	plotiso("captest_oper_c7-0_alliso-gs98.dat", iso, "c7-0", r"$c_{7}^{0}$", "c7-0_plots_"+str(i+1)+"-"+isotopes[i]+".png")

