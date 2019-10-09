import matplotlib.pyplot as plot

# global list to reference the order of the isotopes, and their colours
isotopeList = ["H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24","Al27","Si28","S32","Ar40","Ca40","Fe56","Ni58"]
isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{58}Ni$"]
colours = ["#ef1a1a", "#2fef19", "#0055ff", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670","#ef1a1a", "#0055ff", "#2fef19", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670"]

# plotoper takes the capture rate data for each isotope for a given coupling constant
# that captnoper calculates. This is designed to mimic Catena's plots
def plotoper(couplingConstant, title="Operator Plot", savefigname=None):

	# filename = "Oper_data/captest_oper_"+couplingConstant+"_alliso-gs98.dat"
	# filename = "Oper_factor_of_two_data/factortwotest_oper_"+couplingConstant+"_alliso-gs98.dat"
	# filename = "Oper_factor_of_two_data/factortwotest_oper_"+couplingConstant+"_alliso-b16.dat"
	filename = "Oper_data/captest_oper_"+couplingConstant+"_alliso-ags05.dat"
	# filename = "Oper_data/captest_oper_"+couplingConstant+"_alliso-agss09.dat"
	# filename = "Oper_data/captest_oper_"+couplingConstant+"_alliso-agss09ph.dat"

	# read in the file's data outputed from captn
	# organized into columns of DM mass (x axis), then 17 of isotope specific capture rates
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
	

	# ISOs is now a list of lists: each sub-list has one capture rate entry for each isotope, given the DM mass
	# now create a list of the total capture rate at each DM mass (sum over each sub list)

	totalCap = []
	for subList in ISOs:
		totalCap.append(sum(subList))

	# flip the list of lists around
	# (I need each sub-list to be of one isotope walking through the DM mass to plot it correctly)
	# eg. [[1,2,3],[4,5,6]] --> [[1,4],[2,5],[3,6]]
	ISOs = [list(temp) for temp in zip(*ISOs)]


	plot.clf()
	ax = plot.subplot(111)
	for i in range(len(isotopeList)+1):
		if i == len(isotopeList):
			plot.plot(Ms, totalCap, label="Total", color="Black", marker='.', linestyle='-', linewidth=0.8, markersize=3)
		elif 7<i<len(isotopeList):
			plot.plot(Ms, ISOs[i], color=colours[i], label=isotopeLabels[i], marker='.', linestyle='--', linewidth=0.4, markersize=3)
		else:
			plot.plot(Ms, ISOs[i], color=colours[i],  label=isotopeLabels[i], marker='.', linestyle='-', linewidth=0.4, markersize=3)

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
	# plotoper(c, "Oper results from "+c+" gs98", "Oper_plots/"+c+"_gs98_oper_plot.pdf")
	# plotoper(c, "Oper results from "+c+" gs98", "Oper_factor_of_two_plots/"+c+"_gs98_oper_plot.pdf")
	# plotoper(c, "Oper results from "+c+" b16", "Oper_factor_of_two_plots/"+c+"_b16_oper_plot.pdf")
	plotoper(c, "Oper results from "+c+" ags05", "Oper_plots/"+c+"_ags05_oper_plot.pdf")
	# plotoper(c, "Oper results from "+c+" agss09", "Oper_plots/"+c+"_agss09_oper_plot.pdf")
	# plotoper(c, "Oper results from "+c+" agss09ph", "Oper_plots/"+c+"_agss09ph_oper_plot.pdf")
