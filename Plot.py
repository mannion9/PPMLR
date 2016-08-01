animate   = 1    # Set to 1 to animate set to 0 to get individual plots
step      = 1  	# When animation == 0; Plot steps size 
plotter   = 1 	# When animation == 0; plotter == 0 -> plot total E , plotter = 1 -> plot primatives

import os
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
file_name_e	= [ 'Exact-Output/'+fn for fn in os.listdir(os.getcwd()+'/Exact-Output') if fn.endswith('.txt')]

r      = ReadInData(file_name_m[0])[0]  # CellCenter
rho    = ReadInData(file_name_m[1])     # Density
dt     = ReadInData(file_name_m[2])	    # dt
Energy = ReadInData(file_name_m[3])	    # Energy
energy = ReadInData(file_name_m[4])    	# Internal Energy
press  = ReadInData(file_name_m[5]) 	# Pressure
vel    = ReadInData(file_name_m[6])	    # Velocity

r_e       = ReadInData(file_name_e[0])[0]   # Domain
energy_e  = ReadInData(file_name_e[1])		# Exact internal energy
press_e   = ReadInData(file_name_e[2])		# Exact pressure
rho_e	  = ReadInData(file_name_e[3])		# Exact rho
vel_e	  = ReadInData(file_name_e[4])		# Exact velocity
t       = 0.

if animate == 1:
    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure(1)
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)
    ax1.set_xlabel('r'),ax1.set_ylabel('rho')
    ax2.set_xlabel('r'),ax2.set_ylabel('u')
    ax3.set_xlabel('r'),ax3.set_ylabel('P')
    ax4.set_xlabel('r'),ax4.set_ylabel('E')
    ax1.set_xlim([0.,1.]) , ax2.set_xlim([0.,1.])
    ax3.set_xlim([0.,1.]) , ax4.set_xlim([0.,1.])
    ax1.set_ylim([0.,1.2])  # rho
    ax2.set_ylim([0,2.])    # velocity
    ax3.set_ylim([0,1.])    # pressure
    ax4.set_ylim([1.6,2.8]) # Energy
    
    scat1 = ax1.scatter(r,rho[0])
    scat2 = ax2.scatter(r,vel[0])
    scat3 = ax3.scatter(r,press[0])
    scat4 = ax4.scatter(r,energy[0])
    plt1, = ax1.plot(r_e,rho_e[0])
    plt2, = ax2.plot(r_e,vel_e[0])
    plt3, = ax3.plot(r_e,press_e[0])
    plt4, = ax4.plot(r_e,energy_e[0])
    
    t=0
    # animation function.  This is called sequentially
    def update(i,fig,scat1,scat2,scat3,scat4,plt1,plt2,plt3,plt4,t):
        t += dt[i][0]
        scat1.set_offsets([[r[j],rho[i][j]]    for j in range(len(rho[i])-1)])
        scat2.set_offsets([[r[j],vel[i][j]]    for j in range(len(vel[i])-1)])
        scat3.set_offsets([[r[j],press[i][j]]  for j in range(len(press[i])-1)])
        scat4.set_offsets([[r[j],energy[i][j]] for j in range(len(energy[i])-1)])
        plt1.set_data(r_e,rho_e[i])
        plt2.set_data(r_e,vel_e[i])
        plt3.set_data(r_e,press_e[i])
        plt4.set_data(r_e,energy_e[i])
        plt.suptitle('Time %f' % t)
        return scat1,scat2,scat3,scat4
    
    anim = animation.FuncAnimation(fig, update,fargs=(fig,scat1,scat2,scat3,scat4,plt1,plt2,plt3,plt4,t),frames=len(rho),interval=500)
    anim.save('animation.mp4', fps=10)
    plt.show(1)
    
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
        plt.suptitle('Time %f' % t)
        
        if plotter == 0:
            ax5 = fig.add_subplot(2,j,5)
            ax5.scatter(r,Energy[i],facecolors='none',edgecolor='black')
            ax5.set_xlabel('r'), ax5.set_ylabel('Total Energy')
        else:
            plt.show()

        t += dt[i][0]
    	
