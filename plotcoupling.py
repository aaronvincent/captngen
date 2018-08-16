import matplotlib.pyplot as plot

def plotcoupling(filename, otherfilename, title="Title", savefigname=None):
	# read in the file's data
	file = open(filename,'r')
	currentLine = file.readline()
	Xs=[]
	Ys=[]
	while currentLine != "":
		theLine = currentLine.split()
		Xs.append(float(theLine[0]))
		Ys.append(float(theLine[1]))
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
	plot.plot(Xs, Ys, color="blue", marker='.', linestyle='--', linewidth=0.4, markersize=3, label="Captn_oper")
	plot.plot(otherXs, otherYs, color="Black", linestyle='-', linewidth=0.8, markersize=3, label="Catena Total")
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

plotcoupling("captest_oper_c1-0_test.dat", "c1-0", r"$c_{1}^{0}$", "c1-0_plots_test.png")

