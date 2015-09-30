from numpy import *
from numpy.linalg import det
from matplotlib.pyplot import *
import h5py as hdf
F = hdf.File("anharmonic_1D.hdf5", "r")
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

figure()
plot(q[:,0,:], 'b-')
plot(p[:,0,:], 'r-')
grid(True)
detQ = det(Q)
detP = det(P)
figure()
plot(Q[:,0,0], 'b-')
plot(P[:,0,0], 'r-')
grid(True)
figure()
plot(EPot[:,0,0], 'b-')
plot(EKin[:,0,0],'r-')
plot(EKin[:,0,0] + EPot[:,0,0], 'm-+')
grid(True)
show()

