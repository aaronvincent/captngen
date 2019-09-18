import matplotlib.pyplot as plot

# plotall takes the total capture from the Catena paper (https://arxiv.org/abs/1501.03729) for a
# specified coupling constant (eg c1_0, c5_1, ect...). This is then plotted against the total
# capture rate the code captn_oper finds, and also plots captn's individual isotope capture rates
def plotall(filename, catenafilename, title="Title", savefigname=None):
	
	isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{59}Ni$"]
	colours = ["#ef1a1a", "#2fef19", "#0055ff", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670","#ef1a1a", "#0055ff", "#2fef19", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670"]
	
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
	
	# ISOs is now a list of lists: each sub-list has an isotope-specific capture rate for the given DM mass
	# now create a list of the total capture rate at each DM mass (sum over each sub list)

	totalCap = []
	for subList in ISOs:
		totalCap.append(sum(subList))

	# read in the file's data from catena paper (total capture rate data)
	file = open("Catena_data/"+catenafilename+"_data.dat",'r')
	currentLine = file.readline()
	catenaXs=[]
	catenaYs=[]
	while currentLine != "":
		theLine = currentLine.split()
		catenaXs.append(float(theLine[0]))
		catenaYs.append(float(theLine[1]))
		currentLine = file.readline()
	file.close()

	# flip the list of lists around
	# (I need each sub-list to be of one isotope walking through the DM mass to plot it correctly)
	# eg. [[1,2,3],[4,5,6]] --> [[1,4],[2,5],[3,6]]
	ISOs = [list(temp) for temp in zip(*ISOs)]

	plot.clf()
	ax = plot.subplot(111)
	plot.plot(Ms, totalCap, label="Total", color="Black", marker='.', linestyle='--', linewidth=0.8, markersize=3)
	plot.plot(catenaXs, catenaYs, color="Black", linestyle='-', linewidth=0.8, label="Catena Total")
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

plotall("captest_oper_c1-0_alliso-gs98.dat", "c1-0", r"$c_{1}^{0}$", "c1-0_plots_alliso.png")

