import os 
import numpy as np
import matplotlib.pyplot as plt 

file_name = ["Output_NonUniSpace/"+fn for fn in os.listdir(os.getcwd()+'\Output_NonUniSpace')]
def ReadInData(file):
	content = []
	file = open(file,'r')
	file = file.read().splitlines()
	for line in file:
		line = line.split(' ')
		row  = [float(item) for item in line if item != '']
		content.append(row)
	return content 
 
center = ReadInData(file_name[0])[0]
r_12_i = ReadInData(file_name[1])[0]
dr     = ReadInData(file_name[2])[0]

index_c = np.linspace(1,len(center)+1,len(center))
index_e = np.linspace(1,len(r_12_i)+1,len(r_12_i))

fig = plt.figure()
ax1 = fig.add_subplot(3,1,1)
ax2 = fig.add_subplot(3,1,2)
ax3 = fig.add_subplot(3,1,3)
ax1.set_ylabel('Cell Centers')
ax1.set_xlabel('Index')
ax2.set_ylabel('Cell Edges')
ax2.set_xlabel('Index')
ax3.set_ylabel(' Log Cell Widths')
ax3.set_xlabel('Index')
ax1.scatter(index_c,center)
ax2.scatter(index_e,r_12_i)
ax3.scatter(index_c,np.log(dr))
plt.show()
