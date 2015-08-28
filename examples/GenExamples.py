import numpy.random
import numpy.linalg

from WaveBlocksND import *

from HaWp2YAML import *
from HaWp2CSV import *

def main():
    # Number of dimensions
    D = 3

    shape = LimitedHyperbolicCutShape(D,7.0,(5,5,5))

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
    
    print "Grid: ", grid
    print "Value of wavepacket on grid: "
    print "   ", packet.slim_recursion(grid,0)
    
    print "Values of basis functions on grid: "
    for mindex, entry in zip(list(packet.get_basis_shapes(0).get_node_iterator()), 
                     packet.evaluate_basis_at(grid,0)):
        print "   ", mindex, entry
    
    
    nabla = packet.get_gradient_operator()
    
    print "Values of wavepacket gradient on grid: "
    for index, gradwp in enumerate(nabla.apply_gradient(packet)):
        print "   ", gradwp.slim_recursion(grid,0)
        
        '''
        if index == 2:
            for mindex, entry in zip(list(gradwp.get_basis_shapes(0).get_node_iterator()), 
                     gradwp.get_coefficients(0)):
                print "   ", mindex, entry
        '''
    
    WriteHaWpAsYAML(packet, "3d_scalar_hawp_basis.yaml")
    WriteHaWpCoefficientsAsCSV(packet.get_coefficients()[0], list(packet.get_basis_shapes(0).get_node_iterator()), "3d_scalar_hawp_coeffs.csv")

if __name__ == "__main__":
    numpy.random.seed(0)
    
    main()