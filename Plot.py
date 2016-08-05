animate   = 1    # Set to 1 to animate set to 0 to get individual plots
step      = 1  	# When animation == 0; Plot steps size 
plotter   = 1 	# When animation == 0; plotter == 0 -> plot total E , plotter = 1 -> plot primatives
size      = 6.    # Size of edge of square image in inchest (12 is good for a video)
# Put a limit on plotting domain, even if solution goes crazy at late time 
rho_floor = [-0.5,10.]
u_floor   = [-1.,10.]
P_floor   = [-0.5,10.]
e_floor   = [-0.5,10.]  
#rho_floor = [-0.5,1.2]
#u_floor   = [-1.,3.]
#P_floor   = [-0.5,1.2]
#e_floor   = [-0.5,4.0] 

import os
import numpy as np
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
	
file_name_m	= [ 'Output/'+fn for fn in os.listdir(os.getcwd()+'/Output') if fn.endswith('.txt')]
file_name_e	= [ 'RP-Exact/Output/'+fn for fn in os.listdir(os.getcwd()+'/RP-Exact/Output') if fn.endswith('.txt')]

r         = ReadInData(file_name_m[0])[0]    # CellCenter
rho       = ReadInData(file_name_m[1])       # Density
dt        = ReadInData(file_name_m[2])	     # dt
Energy    = ReadInData(file_name_m[3])	     # Energy
energy    = ReadInData(file_name_m[4])       # Internal Energy
lr        = ReadInData(file_name_m[5])       # Lagrange cell center
press     = ReadInData(file_name_m[6]) 	     # Pressure
Lagrange  = ReadInData(file_name_m[7])[0][0] # RemapChoice
mass      = ReadInData(file_name_m[8])       # Total mass
vel       = ReadInData(file_name_m[9])	     # Velocity

r_e       = ReadInData(file_name_e[0])[0]     # Domain
energy_e  = ReadInData(file_name_e[1])		# Exact internal energy
press_e   = ReadInData(file_name_e[2])		# Exact pressure
rho_e	    = ReadInData(file_name_e[3])		# Exact rho
vel_e	    = ReadInData(file_name_e[4])		# Exact velocity
t         = np.cumsum(dt)				# Time

# Determine limits of plots
## x limits
if Lagrange == 0:
    x_max = max([max(element) for element in lr])  # For each Lagrange step, find the max cell center and find the max of those maximum
    x_min = min([min(element) for element in lr])  # For each Lagrange step, find the min cell center and find the min of those minimum
else:
    x_max = max(r) # Find the max stationary cell center
    x_min = min(r) # Find the min stationary cell center
x_lim = [x_min,x_max]
## y limits
### Density
rho_max = max([max(element) for element in rho])
rho_min = min([min(element) for element in rho])
rho_max += .1*rho_max # Asthetic for plotting 
rho_min -= .1*rho_min 
rho_lim = [max(rho_floor[0],rho_min),min(rho_floor[1],rho_max)]
### Velocity
u_max = max([max(element) for element in vel])
u_min = min([min(element) for element in vel])
u_max += .1*u_max # Asthetic for plotting 
u_min -= .1*u_min 
u_lim = [max(u_floor[0],u_min),min(u_floor[1],u_max)]
### Pressure
P_max = max([max(element) for element in press])
P_min = min([min(element) for element in press])
P_max += .1*P_max # Asthetic for plotting 
P_min -= .1*P_min 
P_lim = [max(P_floor[0],P_min),min(P_floor[1],P_max)]
### energy
e_max = max([max(element) for element in energy])
e_min = min([min(element) for element in energy])
e_max += .1*e_max # Asthetic for plotting 
e_min -= .1*e_min 
e_lim = [max(e_floor[0],e_min),min(e_floor[1],e_max)]

if animate == 1:
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(1)
    fig.set_size_inches(size,size, True)
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    ax1.set_xlabel('r'),ax1.set_ylabel(r'$\rho$')
    ax2.set_xlabel('r'),ax2.set_ylabel('u')
    ax3.set_xlabel('r'),ax3.set_ylabel('P')
    ax4.set_xlabel('r'),ax4.set_ylabel('e')
    ax1.set_xlim(x_lim) , ax2.set_xlim(x_lim)
    ax3.set_xlim(x_lim) , ax4.set_xlim(x_lim)
    ax1.set_ylim(rho_lim)  # rho
    ax2.set_ylim(u_lim)    # velocity
    ax3.set_ylim(P_lim)    # pressure
    ax4.set_ylim(e_lim)    # Energy
    if Lagrange == 0:
        x = lr[0]
    else:
        x = r
    scat1 = ax1.scatter(x,rho[0],marker="+")
    scat2 = ax2.scatter(x,vel[0],marker='+')
    scat3 = ax3.scatter(x,press[0],marker='+')
    scat4 = ax4.scatter(x,energy[0],marker='+')
    plt1, = ax1.plot(r_e,rho_e[0])
    plt2, = ax2.plot(r_e,vel_e[0])
    plt3, = ax3.plot(r_e,press_e[0])
    plt4, = ax4.plot(r_e,energy_e[0])
    
    # animation function.  This is called sequentially
    def update(i,fig,scat1,scat2,scat3,scat4,plt1,plt2,plt3,plt4):
    #def update(i,fig,scat1,scat2,scat3,scat4):    
        if Lagrange == 0:
            x = lr[i]
        else:
            x = r
        scat1.set_offsets([[x[j],rho[i][j]]    for j in range(len(rho[i])-1)])
        scat2.set_offsets([[x[j],vel[i][j]]    for j in range(len(vel[i])-1)])
        scat3.set_offsets([[x[j],press[i][j]]  for j in range(len(press[i])-1)])
        scat4.set_offsets([[x[j],energy[i][j]] for j in range(len(energy[i])-1)])
        plt1.set_data(r_e,rho_e[i])
        plt2.set_data(r_e,vel_e[i])
        plt3.set_data(r_e,press_e[i])
        plt4.set_data(r_e,energy_e[i])
        plt.suptitle(' Step %i Time %f , Total Mass %f ' %( i , t[i] , mass[i][0]))
        return scat1,scat2,scat3,scat4
    #anim = animation.FuncAnimation(fig, update,fargs=(fig,scat1,scat2,scat3,scat4),frames=len(rho),interval=500)
    anim = animation.FuncAnimation(fig, update,fargs=(fig,scat1,scat2,scat3,scat4,plt1,plt2,plt3,plt4),frames=len(rho),interval=500)
    anim.save('Output/out.mp4', fps=1)
    #plt.show(1)
    
else:
    for i in range(0,int(len(rho)),step):
        fig = plt.figure()
        if plotter == 0:  # Determines number colums of plots in subplot
            j = 3 
        else:
            j = 2
        ax1 = fig.add_subplot(2,j,1)
        ax2 = fig.add_subplot(2,j,2)
        ax3 = fig.add_subplot(2,j,3)
        ax4 = fig.add_subplot(2,j,4)
        
        ax1.scatter(r,rho[i]   ,facecolors='none',edgecolor='black')
        ax2.scatter(r,vel[i]   ,facecolors='none',edgecolor='black')
        ax3.scatter(r,press[i] ,facecolors='none',edgecolor='black')
        ax4.scatter(r,energy[i],facecolors='none',edgecolor='black')

        ax1.plot(r_e,rho_e[i])
        ax2.plot(r_e,vel_e[i])
        ax3.plot(r_e,press_e[i])
        ax4.plot(r_e,energy_e[i])
        
        ax1.set_xlabel('r') , ax1.set_ylabel('Density')
        ax2.set_xlabel('r') , ax2.set_ylabel('Velocity')
        ax3.set_xlabel('r') , ax3.set_ylabel('Pressure')
        ax4.set_xlabel('r') , ax4.set_ylabel('Internal Energy')
        plt.suptitle('Time %f' % t[i])
        
        if plotter == 0:
            ax5 = fig.add_subplot(2,j,5)
            ax5.scatter(r,Energy[i],facecolors='none',edgecolor='black')
            ax5.set_xlabel('r'), ax5.set_ylabel('Total Energy')
        else:
            plt.show()

        t += dt[i][0]
    	
