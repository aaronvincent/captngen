import numpy as np
import matplotlib.pyplot as plot
import os
os.chdir(os.path.dirname(__file__)) # cd into this script's directory

def readtimingfromcaptn(dirandfilename):
	'''
	reads in the timing test files from captnoper, where all isotopes are already summed over
	it ran with 5 gev and 50 gev DM masses for each of the isoscalar coupling constants
	it also has the output from the linux 'time' command to determine run speed
	funtion returns a list of the capture rates for each coupling constant at both masses,
	 the coupling constants themselves, and the time to run gentest.x in seconds
	'''
	couplinglist = []
	fiveGeVlist = []
	fiftyGeVlist = []
	time = 0
	with open(dirandfilename) as f:
		couplingline = ""
		while couplingline != "c15-0":
			couplingline = f.readline().strip()
			fiveGeVline = f.readline().strip().split()
			fiftyGeVline = f.readline().strip().split()
			f.readline() # blank line to skip
		
			couplinglist.append(couplingline)
			fiveGeVlist.append(float(fiveGeVline[1]))
			fiftyGeVlist.append(float(fiftyGeVline[1]))
		
		f.readline() # blank line to skip
		timeline = f.readline().strip().split()
		timesplit = timeline[1].split('m')
		minutes = int(timesplit[0])
		seconds = float(timesplit[1].strip('s'))
		time = 60*minutes + seconds
		f.close()
	return [time, np.array(couplinglist), np.array(fiveGeVlist), np.array(fiftyGeVlist)]

def computeErrors(cutchoice, cutnormal):
	timepercall = cutchoice[0]/(14*2)
	# print(timepercall)
	cutchoice.append(timepercall)

	errorpercent5GeV = 100 * abs(cutchoice[2]-cutnormal[2])/cutnormal[2]
	# print(errorpercent5GeV)
	cutchoice.append(errorpercent5GeV)
	
	errorpercent50GeV = 100 * abs(cutchoice[3]-cutnormal[3])/cutnormal[3]
	# print(errorpercent50GeV)
	cutchoice.append(errorpercent50GeV)


dirname = "./"
# filename = "testacc_cut1-ags05.dat"

# cut1 = readtimingfromcaptn(dirname+"testacc_cut1-ags05.dat")
# cut2 = readtimingfromcaptn(dirname+"testacc_cut2-ags05.dat")
# cut5 = readtimingfromcaptn(dirname+"testacc_cut5-ags05.dat")
# cut10 = readtimingfromcaptn(dirname+"testacc_cut10-ags05.dat")
# cut100 = readtimingfromcaptn(dirname+"testacc_cut100-ags05.dat")

# computeErrors(cut1,cut1)
# computeErrors(cut2,cut1)
# computeErrors(cut5,cut1)
# computeErrors(cut10,cut1)
# computeErrors(cut100,cut1)

# xAxis = [1,2,5,10,100]

# plot.clf()
# plot.xlabel("Save only every Nth line")
# plot.ylabel("Call time")
# plot.title("CaptnOper time when Solar Model Files are Cut in Size")
# for i in range(len(cut1[1])):
# 	plot.plot(xAxis, [cut1[4], cut2[4], cut5[4], cut10[4], cut100[4]], label=cut1[1][i])#, color="", marker=".", linestyle="")
# plot.legend()
# plot.savefig("figTime.pdf")
# plot.show()

# plot.clf()
# plot.xlabel("Save only every Nth line")
# plot.ylabel("Error Percent")
# plot.title("Error % from CaptnOper when Solar Model Files are Cut in Size for DM 5 GeV")
# for i in range(len(cut1[1])):
# 	plot.plot(xAxis, [cut1[5][i], cut2[5][i], cut5[5][i], cut10[5][i], cut100[5][i]], label=cut1[1][i])#, color="", marker=".", linestyle="")
# plot.legend()
# plot.savefig("fig5GeV.pdf")
# plot.show()

# plot.clf()
# plot.xlabel("Save only every Nth line")
# plot.ylabel("Error Percent")
# plot.title("Error % from CaptnOper when Solar Model Files are Cut in Size for DM 50 GeV")
# for i in range(len(cut1[1])):
# 	plot.plot(xAxis, [cut1[6][i], cut2[6][i], cut5[6][i], cut10[6][i], cut100[6][i]], label=cut1[1][i])#, color="", marker=".", linestyle="")
# plot.legend()
# plot.savefig("fig50GeV.pdf")
# plot.show()

cut1 = readtimingfromcaptn(dirname+"testacc-agss09.dat")
cut10lin = readtimingfromcaptn(dirname+"testacc_10lin-agss09.dat")
cut10log = readtimingfromcaptn(dirname+"testacc_10log-agss09.dat")
cut20lin = readtimingfromcaptn(dirname+"testacc_20lin-agss09.dat")
cut20log = readtimingfromcaptn(dirname+"testacc_20log-agss09.dat")

computeErrors(cut1,cut1)
computeErrors(cut10lin,cut1)
computeErrors(cut10log,cut1)
computeErrors(cut20lin,cut1)
computeErrors(cut20log,cut1)

xAxis = [1,10,20]

plot.clf()
plot.xlabel("Save only every Nth line")
plot.ylabel("Call time")
plot.title("CaptnOper time when Solar Model File agss09 is Cut in Lin Size")
for i in range(len(cut1[1])):
	plot.plot(xAxis, [cut1[4], cut10lin[4], cut20lin[4]], label=cut1[1][i], color=colours[i], marker="x", linestyle="--")
plot.legend()
plot.savefig("figTime_agss09_lin.pdf")
plot.show()

plot.clf()
plot.xlabel("Save only every Nth line")
plot.ylabel("Call time")
plot.title("CaptnOper time when Solar Model File agss09 is Cut in Log Size")
for i in range(len(cut1[1])):
	plot.plot(xAxis, [cut1[4], cut10log[4], cut20log[4]], label=cut1[1][i], color=colours[i], marker="x", linestyle="--")
plot.legend()
plot.savefig("figTime_agss09_log.pdf")
plot.show()

plot.clf()
plot.xlabel("Save only every Nth line")
plot.ylabel("Error Percent")
plot.title("Error % from CaptnOper when Solar Model File agss09 is Cut in Lin Size for DM 5 GeV")
for i in range(len(cut1[1])):
	plot.plot(xAxis, [cut1[5][i], cut10lin[5][i], cut20lin[5][i]], label=cut1[1][i], color=colours[i], marker="x", linestyle="--")
plot.legend()
plot.savefig("fig5GeV_agss09_lin.pdf")
plot.show()

plot.clf()
plot.xlabel("Save only every Nth line")
plot.ylabel("Error Percent")
plot.title("Error % from CaptnOper when Solar Model File agss09 is Cut in Log Size for DM 5 GeV")
for i in range(len(cut1[1])):
	plot.plot(xAxis, [cut1[5][i], cut10log[5][i], cut20log[5][i]], label=cut1[1][i], color=colours[i], marker="x", linestyle="--")
plot.legend()
plot.savefig("fig5GeV_agss09_log.pdf")
plot.show()

plot.clf()
plot.xlabel("Save only every Nth line")
plot.ylabel("Error Percent")
plot.title("Error % from CaptnOper when Solar Model File agss09 is Cut in Lin Size for DM 50 GeV")
for i in range(len(cut1[1])):
	plot.plot(xAxis, [cut1[6][i], cut10lin[6][i], cut20lin[6][i]], label=cut1[1][i], color=colours[i], marker="x", linestyle="--")
plot.legend()
plot.savefig("fig50GeV_agss09_lin.pdf")
plot.show()

plot.clf()
plot.xlabel("Save only every Nth line")
plot.ylabel("Error Percent")
plot.title("Error % from CaptnOper when Solar Model File agss09 is Cut in Log Size for DM 50 GeV")
for i in range(len(cut1[1])):
	plot.plot(xAxis, [cut1[6][i], cut10log[6][i], cut20log[6][i]], label=cut1[1][i], color=colours[i], marker="x", linestyle="--")
plot.legend()
plot.savefig("fig50GeV_agss09_log.pdf")
plot.show()