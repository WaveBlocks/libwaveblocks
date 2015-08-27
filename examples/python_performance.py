from WaveBlocksND import *
import numpy.random
import time

def run(D):
    # Number of dimensions
    shape = LimitedHyperbolicCutShape(D,2**D,[4]*D)
    
    packet = HagedornWavepacket(D,1,0.9)
    
    # PARAMETERS (All components share same parameter set)
    p = numpy.random.randn(D)
    q = numpy.random.randn(D)
    P = numpy.random.randn(D,D) + 1.0j*numpy.random.randn(D,D)
    Q = numpy.random.randn(D,D) + 1.0j*numpy.random.randn(D,D)
    packet.set_parameters((q,p,Q,P), key=('q','p','Q','P'))
    
    # SHAPE
    packet.set_basis_shapes((shape),0)
    
    # COEFFICIENTS
    basis_size = shape.get_basis_size()
    coeffs = numpy.random.randn(basis_size) + 1.0j*numpy.random.randn(basis_size)
    coeffs /= numpy.linalg.norm(coeffs)
    packet.set_coefficients((coeffs,))
    
    # EVALUATE
    grid = numpy.linspace(-1.0, 1.0,D)
    
    print "Evaluate wavepacket"
    start = time.clock()
    result = packet.slim_recursion(grid,0)
    end = time.clock()
    print "   shape: ", shape
    print "   size of shape: ", shape.get_basis_size()
    print "   value: ", result
    print "   time:  ", str((end - start)*1000), "[ms]"
    
    #print "Evaluate gradient: "
    #nabla = packet.get_gradient_operator()
    #for index, gradwp in enumerate(nabla.apply_gradient(packet)):
        #print "   ", gradwp.slim_recursion(grid,0)

def main():
    run(5)
    run(8)

if __name__ == "__main__":
    numpy.random.seed(0)
    
    main()