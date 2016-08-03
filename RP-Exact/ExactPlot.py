import matplotlib.pyplot as plt 
domain_f  = open('Output/domain.txt','r')
pressure_f= open('Output/pressure.txt','r')
rho_f     = open('Output/rho.txt','r')
velocity_f= open('Output/velocity.txt','r')
energy_f  = open('Output/energy.txt','r')

def ReadInData(files):
    ''' Input is a list of the files that have been opened.
        This will break up the content of the files by row and column'''
    content =[] 
    for data in files:
        data = data.read().splitlines()
        events = []
        for line in data:
            row = []
            line = line.split(' ')
            for i in line:
                if i != '': 
                    row.append(float(i))
            events.append(row)
        content.append(events)
    return content
Data = ReadInData([domain_f,pressure_f,rho_f,velocity_f,energy_f])

step = 10 # only plot ever 10 itterations

x       = Data[0][0]
pressure= Data[1]
rho     = Data[2]
velocity= Data[3]
energy  = Data[4]

fig = plt.figure()
ax1 = fig.add_subplot(2,2,1)
ax1.set_xlabel('Position')
ax1.set_ylabel('Density')
ax2 = fig.add_subplot(2,2,2)
ax2.set_xlabel('Position')
ax2.set_ylabel('Velocity')
ax3 = fig.add_subplot(2,2,3)
ax3.set_xlabel('Position')
ax3.set_ylabel('Pressure')
ax4 = fig.add_subplot(2,2,4)
ax4.set_xlabel('Position')
ax4.set_ylabel('Energy')

ax1.plot(x,rho)
ax2.plot(x,velocity)
ax3.plot(x,pressure)
ax4.plot(x,energy)

plt.show()
