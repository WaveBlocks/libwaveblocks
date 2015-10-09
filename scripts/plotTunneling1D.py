
# coding: utf-8

# In[1]:


# In[2]:

from numpy import *
from numpy.linalg import det
from matplotlib.pyplot import *
from matplotlib import gridspec
import sys
sys.path.insert(0, 'WaveBlocksND/src/WaveBlocksND')
from Plot import stemcf

# In[3]:

import h5py as hdf

K = 512
# In[6]:

Fs = [hdf.File("tunneling_1D_" +str(x) + ".hdf5", "r") for x in xrange(0,14)]



# In[7]:

dt = 0.005
T = 70
time = linspace(0, T, T/dt+1)[:-1]


# In[8]:$

f = lambda x: ('%f' % x).rstrip('0').rstrip('.')


# In[9]:

def read_data(base, time):
  L = []
  for x in xrange(0,14):
    L.extend([Fs[x][base + "@" + f(t)][:] for t in time[1000*x:1000*(x+1)]])
  return array(L)

# In[10]:
q = read_data("q", time)
p = read_data("p", time)
Q = read_data("Q", time)
P = read_data("P", time)


# In[11]:

q.shape


# In[12]:
fig = figure()
gs = gridspec.GridSpec(2,2)
ax1 = fig.add_subplot(gs[0,0])
ax1.plot(time, q[:,0,:], 'b-', label = '$q$')
xlabel('time')
legend()
grid(True)
ax2 =  fig.add_subplot(gs[0,1])

ax2.plot(time, p[:,0,:], 'r-', label= '$p$')
xlabel('time')
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

figure()
Etot = squeeze(Etot);
semilogy(time, abs(Etot - Etot[0]), 'm-', label='$|E_{kin}(0) - E_{kin}(t)|$')
xlabel('time')
legend(loc="upper left", bbox_to_anchor=(1,1))
grid(True)




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

figure()
subplot(2,3,1)
title('$t=0$')
c = C[200*0, :,0];
stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)


# In[29]:

subplot(2,3,2)
title('$t=20$')
c = C[200*20, :,0];
stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)


# In[30]:

subplot(2,3,3)
title('$t=30$')
c = C[200*30, :,0];
stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)


# In[31]:

subplot(2,3,4)
title('$t=40$')
c = C[200*40, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)


# In[27]:

subplot(2,3,5)
title('$t=50$')
c = C[200*50, :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)

subplot(2,3,6)
title('$t=70$')
c = C[200*70 -1 , :,0];

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(-1, K + 1)
ylim(-0.1, 1.1)


# 51.63
c = C[10326,:,0]
figure()
subplot(2,1,1)
title('$t=51.63$')

stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(50, K + 1)
ylim(0.0, 0.1)

# 61.96
c = C[12392,:,0]
subplot(2,1,2)
title('$t=61.96$')
stemcf(arange(K), angle(c), abs(c))
grid(True)
xlim(50, K + 1)
ylim(0.0, 0.1)

show()
