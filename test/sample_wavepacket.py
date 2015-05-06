import numpy
import numpy.linalg
import math
from WaveBlocksND import *

def createSampleWavepacket(D, shape):
    packet = HagedornWavepacket(D,1,0.9)
    
    # create sample parameters
    p = numpy.zeros((D))
    q = numpy.zeros((D))
    Q = numpy.eye(D, dtype=complex)
    P = 1j*numpy.eye(D, dtype=complex)
    for i in range(0,D):
        q[i] = math.cos(i+1)
        p[i] = math.sin(i+1)
        for j in range(0,D):
            Q[i,j] += 0.3*math.sin(i+j+1) + 0.3j*math.cos(i+j+1)
            P[i,j] += 0.3*math.cos(i+j+1) + 0.3j*math.sin(i+j+1)
    packet.set_parameters( (q,p,Q,P), key=("q","p","Q","P"))
    
    # set shape
    packet.set_basis_shapes((shape),0)
    
    # create sample coefficients
    for node in shape.get_node_iterator():
        falloff = 0.1

        x = 0.0
        y = 0.0
        s = 0
        for d in range(0,D):
            x += math.sin(node[d] + 1.0*(d+1)/D)
            y += math.cos(node[d] + 1.5*(d+1)/D)
            s += node[d]

        x *= math.exp(-falloff*s)
        y *= math.exp(-falloff*s)

        packet.set_coefficient(0, node, x + y*1j)
    
    # normalize coefficients
    norm = numpy.linalg.norm(packet.get_coefficients(0))**2
    
    packet.set_coefficient_vector(packet.get_coefficient_vector() / norm)
    
    return packet