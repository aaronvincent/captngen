import numpy as np

cutnumber = 100

foldername_model = "solarmodels/Models_From_DarkSUSY/"
filename_model = "Serenelli-model_ags05.dat"
myfilename = "accuracytesting/Cut"+str(cutnumber)+"_"+filename_model

with open(myfilename, 'w') as myfile:
	with open(foldername_model+filename_model, 'r') as modelfile:
		# Skip header of length head
		print("header:")
		head = 20
		for i in range(head):
			line = modelfile.readline()
			# print(i, line)
			# myfile.write(line)
		
		print("body:")
		i=0
		line = modelfile.readline()
		while line != "":
			if i%cutnumber == 0: 
				print(i+head, line)
				myfile.write(line)
			i += 1
			line = modelfile.readline()
