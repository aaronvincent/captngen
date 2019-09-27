import matplotlib.pyplot as plot

# global list to reference the order of the isotopes, and their colours
isotopeList = ["H","He3","He4","C12","N14","O16","Ne20","Na23","Mg24","Al27","Si28","S32","Ar40","Ca40","Fe56","Ni58"]
isotopeLabels = [r"$H$",r"$^{3}He$",r"$^{4}He$",r"$^{12}C$",r"$^{14}N$",r"$^{16}O$",r"$^{20}Ne$",r"$^{23}Na$",r"$^{24}Mg$",r"$^{27}Al$",r"$^{28}Si$",r"$^{32}S$",r"$^{40}Ar$",r"$^{40}Ca$",r"$^{56}Fe$",r"$^{58}Ni$"]
colours = ["#ef1a1a", "#2fef19", "#0055ff", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670","#ef1a1a", "#0055ff", "#2fef19", "#137708", "#00fff2", "#ff5efc", "#9b9b9b", "#080670"]


# plottogether takes the capture rates from the Catena paper (https://arxiv.org/abs/1501.03729) and 
# from CaptnOper for a specified coupling constant (eg c1-0, c5-1, ect...). This is designed to 
# directly compare the results between Catena and CaptnOper
def plottogether(couplingConstant, title="Double Plot", savefigname=None):
	
	# Read in the Catena Data
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


	# Read in the CaptnOper data
	filename = "Oper_data/captest_oper_"+couplingConstant+"_alliso-gs98.dat"
	#filename = "../factortwotest_oper_"+couplingConstant+"_alliso-gs98.dat"#
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
	# plotting all of the Catena Data
	for i in range(len(isotopeList)+1):
		if i == len(isotopeList):
			plot.plot(catenaMs[i], catenaCs[i], color="Black", linestyle='-.', linewidth=0.8, label="Catena Total")
		elif 7<i<len(isotopeList):
			plot.plot(catenaMs[i], catenaCs[i], color=colours[i], label="Catena "+isotopeLabels[i], linestyle=':', linewidth=0.4)
		else:
			plot.plot(catenaMs[i], catenaCs[i], color=colours[i], label="Catena "+isotopeLabels[i], linestyle='-.', linewidth=0.4)
	# plotting all of the CaptnOper data
	for i in range(len(isotopeList)+1):
		if i == len(isotopeList):
			plot.plot(Ms, totalCap, label="CaptnOper Total", color="Black", marker='.', linestyle='-', linewidth=0.8, markersize=3)
		elif 7<i<len(isotopeList):
			plot.plot(Ms, ISOs[i], color=colours[i], label="CaptnOper "+isotopeLabels[i], marker='.', linestyle='--', linewidth=0.4, markersize=3)
		else:
			plot.plot(Ms, ISOs[i], color=colours[i],  label="CaptnOper "+isotopeLabels[i], marker='.', linestyle='-', linewidth=0.4, markersize=3)

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


# same as plottogether, but uses the csv from WebPlotDigitizer
def plottogether_csv(couplingConstant, title="Catena Plot", savefigname=None):
	
	# Read in the Catena Data
	filename = "Catena_data_updated/Catena_"+couplingConstant+".csv"
	file = open(filename,'r')
	# skip two header lines
	currentLine = file.readline()
	currentLine = file.readline()
	currentLine = file.readline()
	catenaMs = []
	catenaCs = []
	while currentLine != "":
		theLine = currentLine.split(",")
		if len(theLine) > 1:
			Ms = []
			Cs = []
			for i in range(len(isotopeList)+1): # loop through each isotope + total capture rate
				# the csv file WebPlotDigitizer outputs is X,Y,X,Y,X,Y... ect for each isotope
				# so I need to pull the DM mass and capture rate for each isotope on each line I read
				# Python doesn't like converting the empty string into a None, so I did it manually
				M = theLine[0+2*i]
				C = theLine[1+2*i]
				C = C.rstrip("\n") # the csv has newline characters on the final capture rate of each line

				if M == "":
					M = None
				else:
					M = float(M)

				if C == "":
					C = None
				else:
					C = float(C)

				Ms.append(M)
				Cs.append(C)

		catenaMs.append(Ms)
		catenaCs.append(Cs)
		currentLine = file.readline()
	file.close()

	# now I need to flip the lists around (the first entry of each sublist in catenaMs and catenaCs corresponds to Total,
	# the next is H, then He3, ect ..., Ni58). I want each sublist to be one isotope
	# eg. [[1,2,3],[4,5,6]] --> [[1,4],[2,5],[3,6]]
	catenaMs = [list(temp) for temp in zip(*catenaMs)]
	catenaCs = [list(temp) for temp in zip(*catenaCs)]



	# Read in the CaptnOper data
	#filename = "Oper_data/captest_oper_"+couplingConstant+"_alliso-gs98.dat"
	filename = "Oper_factor_of_two_data/factortwotest_oper_"+couplingConstant+"_alliso-gs98.dat"
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
	# plotting all of the Catena Data
	plot.plot(catenaMs[0], catenaCs[0], color="Black", linestyle='-.', linewidth=0.8, label="Total")
	for i in range(len(isotopeList)):
		if i<8:
			plot.plot(catenaMs[i+1], catenaCs[i+1], color=colours[i], label=isotopeLabels[i], linestyle='-.', linewidth=0.4)
		else:
			plot.plot(catenaMs[i+1], catenaCs[i+1], color=colours[i], label=isotopeLabels[i], linestyle=':', linewidth=0.4)
	# plotting all of the CaptnOper data
	for i in range(len(isotopeList)+1):
		if i == len(isotopeList):
			plot.plot(Ms, totalCap, label="CaptnOper Total", color="Black", marker='.', linestyle='-', linewidth=0.8, markersize=3)
		elif 7<i<len(isotopeList):
			plot.plot(Ms, ISOs[i], color=colours[i], label="CaptnOper "+isotopeLabels[i], marker='.', linestyle='--', linewidth=0.4, markersize=3)
		else:
			plot.plot(Ms, ISOs[i], color=colours[i],  label="CaptnOper "+isotopeLabels[i], marker='.', linestyle='-', linewidth=0.4, markersize=3)

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
	#plottogether_csv(c, "Comparison plot of "+c, "Comparison_Plots/"+c+"_comparison.png")
	plottogether_csv(c, "Comparison plot of "+c, "New_Comparisons/"+c+"_comparison.png")
