import matplotlib.pyplot as plt
import numpy as np

print('Uniform Spacing:')
x_min = 0.
x_max = 1.
N = 100

dx_uni = (x_max-x_min)/(N-1)

r_12_i = [x_min+i*dx_uni             for i in range(N+1)]
dr     = [r_12_i[i+1]-r_12_i[i]      for i in range(N)]
center = [.5*(r_12_i[i+1]+r_12_i[i]) for i in range(N)]

index_c = np.linspace(1,N,N)
index_e = np.linspace(1,N+1,N+1)

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
#plt.show()


print('Non-Uniform Spacing:')
x_min = 0.
x_max = 1.
N = 100

# Shock region
f  = .90 
xL = .4
xR = .6
NS = int(f*N)
dxS= (xR-xL)/(NS-1)

# Regular region
Nremain = N - NS
NL = int(.5*Nremain)
x_max_L = xL - dxS
x_min_L = x_min
NR = int(.5*Nremain)
x_min_R = xR + dxS
x_max_R = x_max

dx_uni_L = (x_max_L-x_min_L)/(NL-1)
dx_uni_R = (x_max_R-x_min_R)/(NR-1)

r_12_i_L = [x_min_L+i*dx_uni_L         for i in range(NL)]
r_12_i_S = [xL+i*dxS                   for i in range(NS)]
r_12_i_R = [x_min_R+i*dx_uni_R         for i in range(NR+1)]

r_12_i = r_12_i_L + r_12_i_S + r_12_i_R
    
dr     = [r_12_i[i+1]-r_12_i[i]      for i in range(len(r_12_i)-1)]
center = [.5*(r_12_i[i+1]+r_12_i[i]) for i in range(len(r_12_i)-1)]

index_c = np.linspace(1,len(center),len(center))
index_e = np.linspace(1,len(r_12_i),len(r_12_i))

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
#plt.show()

print('NS: %i' %NS)
print('dxS: %f ' %dxS)
print('dx_uni_R: %f ' %dx_uni_R)
print('dx_uni_l %f ' %dx_uni_L)