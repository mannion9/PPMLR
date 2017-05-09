Output = '0'            # Output [0 = create .png files        , 1 = create animation pop up  , 2 = create .mp4 file]
height,width = 12.,15.  # Height width of plot in inches (12 is good for a video)
floor  = [-3000.,3000.] # Limits on plotting axes
pad    = 4
pad    = [ i for i in range(0,pad)]

import os
import sys
from matplotlib import pyplot as plt
from matplotlib import animation

def ReadInData(file):
	content = []
	file = open(file,'r')
	file = file.read().splitlines()
	for line in file:
		line = line.split(' ')
		row  = [float(item) for item in line if item != '']
		content.append(row)
	return content

def limitFinder(array,frac):
	mini = min([min(element) for element in array])
	maxi = max([max(element) for element in array])
	mini -= frac*abs(mini)
	maxi += frac*abs(maxi)
	limits = [max(floor[0],mini),min(floor[1],maxi)]
	return limits

def Plot_set(ax1,ax2,ax3,ax4,i=0):
	ax1.set_xlabel('r'),ax1.set_ylabel(r'$\rho$')
	ax2.set_xlabel('r'),ax2.set_ylabel('u')
	ax3.set_xlabel('r'),ax3.set_ylabel('P')
	ax4.set_xlabel('r'),ax4.set_ylabel('e')
	plt.suptitle('Time %f ' % t[i][0] )
	ax1.set_xlim(x_lim)   , ax2.set_xlim(x_lim)
	ax3.set_xlim(x_lim)   , ax4.set_xlim(x_lim)
	ax1.set_ylim(rho_lim) , ax2.set_ylim(u_lim)
	ax3.set_ylim(P_lim)   , ax4.set_ylim(e_lim)
	ax1.grid()            ,ax2.grid()
	ax3.grid()            ,ax4.grid()
	return

def Plot_set_single(ax1,i=0):
	ax1.set_xlabel('r',fontsize=20),ax1.set_ylabel(r'$\rho$',fontsize=20)
	plt.suptitle(' Time %f  ' %( t[i][0] ))
	ax1.set_xlim(x_lim)  , ax1.set_ylim(rho_lim) , ax1.grid()
	return

def SavePng(plt,ax1,ax2,ax3,ax4,j):
	for i in pad:
		if 10**(i) > j and j<10**(i+1):
			plt.savefig('Output/Pictures/Image'+(pad[-1]-i)*'0'+'%i.png' % j)
			break
	ax1.cla(),ax2.cla(),ax3.cla(),ax4.cla()
	return

def SavePng_single(plt,ax1,j):
	for i in pad:
		if 10**(i) > j and j<10**(i+1):
			plt.savefig('Output/Pictures/Image'+(pad[-1]-i)*'0'+'%i.png' % j)
			break
	ax1.cla()
	return

# Read in the data
file_name = sorted([ 'Output/'+fn for fn in os.listdir(os.getcwd()+'/Output') if fn.endswith('.txt')],key=str.lower)
energy    = ReadInData(file_name[0])  # Exact internal energy
lrc       = ReadInData(file_name[1])  # Exact positions
press     = ReadInData(file_name[2])  # Exact pressure
rho	      = ReadInData(file_name[3])  # Exact rho
vel  	  = ReadInData(file_name[4])  # Exact velocity
file_name = sorted([ 'Inputs/'+fn for fn in os.listdir(os.getcwd()+'/Inputs') if fn.endswith('.txt')],key=str.lower)
t         = ReadInData(file_name[0])

# Determine plotting limits
x_lim   = limitFinder(lrc,0.)
rho_lim = limitFinder(rho,.1)
u_lim   = limitFinder(vel,.1)
P_lim   = limitFinder(press,.1)
e_lim   = limitFinder(energy,.1)
# x_lim   = [0.4,1.0]
# rho_lim = [0.,4.2]

# Create figure to plot to
if Output!='3':
	fig = plt.figure(1)
	fig.set_size_inches(width,height, True)
	ax1 = fig.add_subplot(2,2,1)
	ax2 = fig.add_subplot(2,2,2)
	ax3 = fig.add_subplot(2,2,3)
	ax4 = fig.add_subplot(2,2,4)
else:
	fig = plt.figure(1)
	fig.set_size_inches(width,height, True)
	ax1 = fig.add_subplot(111)

# Create images
if Output == '0':
	for i in range(len(rho)):
		# Create line plot objects
		Plot_set(ax1,ax2,ax3,ax4,i=i)
		ax1.plot(lrc[i],rho[i])
		ax2.plot(lrc[i],vel[i])
		ax3.plot(lrc[i],press[i])
		ax4.plot(lrc[i],energy[i])
		j = int(i+1)
		SavePng(plt,ax1,ax2,ax3,ax4,j)

# Create image of only Density
# Create images
if Output == '3':
	for i in range(len(rho)):
		# Create line plot objects
		Plot_set_single(ax1,i=i)
		sol = ax1.plot(lrc[i],rho[i],label='Exact')
		j = int(i+1)
		SavePng_single(plt,ax1,j)

# Create video of animation
if Output == '1' or Output == '2':
		# Create intial figure
		Plot_set(ax1,ax2,ax3,ax4,i=0)
		plt1, = ax1.plot(lrc[0],rho[0])
		plt2, = ax2.plot(lrc[0],vel[0])
		plt3, = ax3.plot(lrc[0],press[0])
		plt4, = ax4.plot(lrc[0],energy[0])


		def updateExact(i,fig,plt1,plt2,plt3,plt4):
			plt1.set_data(lrc[i],rho[i])
			plt2.set_data(lrc[i],vel[i])
			plt3.set_data(lrc[i],press[i])
			plt4.set_data(lrc[i],energy[i])
			plt.suptitle(' Time %f' % t[i][0])
			Plot_set(ax1,ax2,ax3,ax4,i=i)
			return plt1,plt2,plt3,plt4

		anim = animation.FuncAnimation(fig, updateExact,fargs=(fig,plt1,plt2,plt3,plt4),frames=len(rho),interval=500)

		if Output == '1':
			plt.show()
		elif Output == '2':
			anim.save('Output/out.mp4', fps=2)
