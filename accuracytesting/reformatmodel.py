import numpy as np

cutnumber = 10

foldername_model = "solarmodels/"
filename_model = "struct_b16_gs98.dat"
myfilename = "struct_b16_gs98_reduce10.dat"

with open(myfilename, 'w') as myfile:
	with open(foldername_model+filename_model, 'r') as modelfile:
		# Skip header of length head
		print("header:")
		head = 9
		for i in range(head):
			line = modelfile.readline()
			# print(i, line)
			myfile.write(line)
		
		print("body:")
		i=0
		line = modelfile.readline()
		while line != "":
			if i%cutnumber == 0: 
				# print(i+head, line)
				myfile.write(line)
			i += 1
			line = modelfile.readline()
