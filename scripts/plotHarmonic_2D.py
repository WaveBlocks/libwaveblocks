from numpy import *
from numpy.linalg import det
from matplotlib.pyplot import *
from matplotlib import gridspec
import h5py as hdf
import sys
sys.path.insert(0, 'IMPORTANT/WaveBlocksND/src/WaveBlocksND')
from Plot import stemcf

F = hdf.File("harmonic_2D.hdf5", "r")
dt = 0.01
T = 12
time = linspace(0, T, T/dt+1)
f = lambda x: ('%f' % x).rstrip('0').rstrip('.')

def read_data(base, time):
	L = [F[base + "@" + f(t)][:] for t in time]
	return array(L)

q = read_data("q", time)
p = read_data("p", time)
Q = read_data("Q", time)
P = read_data("P", time)
EPot = read_data("Epot",time)
EKin = read_data("Ekin",time)
Etot = EPot + EKin

figure()
Etot = squeeze(Etot);
semilogy(time, abs(Etot - Etot[0]), 'm-', label='$|E_{kin}(0) - E_{kin}(t)|$')
xlabel('time')
legend(loc="upper left", bbox_to_anchor=(1,1))
grid(True)


fig = figure()
gs = gridspec.GridSpec(2,2)
ax1 = fig.add_subplot(gs[0,0])
ax1.plot(q[:,0,:],q[:,1,:], 'b-', label = '$(q_x,q_y)$')
legend()
grid(True)
ax2 =  fig.add_subplot(gs[0,1])

ax2.plot(p[:,0,:],p[:,1,:], 'r-', label= '$(p_x,p_y)$')
legend()
grid(True)




# In[13]:

detQ = det(Q)
detP = det(P)


# In[14]:

fig.add_subplot(gs[1,0])
plot(time, abs(detQ), 'b-', label= '$abs(det(Q))$')
xlabel('time')
legend()
grid(True)
fig.add_subplot(gs[1,1])
plot(time, abs(detP), 'r-', label = '$abs(det(P))$')
xlabel('time')
legend()
grid(True)



# In[15]:
Ekin = read_data("Ekin", time)
Epot = read_data("Epot", time)
Etot = read_data("Etot", time)

# In[16]:

Ekin.shape


# In[17]:

figure()
plot(time, squeeze(Ekin), 'r-', label='$E_{kin}$')
plot(time, squeeze(Epot), 'b-', label='$E_{pot}$' )
plot(time, squeeze(Etot), 'm-', label='$E_{pot}$')
xlabel('time')
legend(loc="upper left", bbox_to_anchor=(1,1))
grid(True)


# In[18]:

C = read_data("coefficients", time)


# In[19]:

C.shape


#~ # In[28]:
K = len(C[0,:,0]);

figure()
subplot(2,3,1)
title('$t=0$')
c = C[0, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)


# In[29]:


subplot(2,3,2)
title('$t=4$')
c = C[400, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)



# In[30]:

subplot(2,3,3)
title('$t=6$')
c = C[600, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)

# In[31]:

subplot(2,3,4)
title('$t=8$')
c = C[800, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)
# In[27]:

subplot(2,3,5)
title('$t=10$')
c = C[1000, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)

subplot(2,3,6)

title('$t=12$')
c = C[1200, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)

show()

