from numpy import *
from time import time

from WaveBlocksND import *

def run1D():
    order = 8
    n_runs = 100

    print
    print "1D homogeneous quadrature of order", str(order) + "."
    print "number of coefficients vs. evaluation time [ms]", \
          "averaged over", n_runs, "runs:"

    for N in [ order * 2**i for i in range(3, 9) ]:
        K = HyperCubicShape([N])
        HAWP = HagedornWavepacket(1, 1, 0.2)
        HAWP.set_basis_shapes([K])
        HAWP.set_coefficients([ones((N, 1))])

        QR = GaussHermiteQR(order)
        IP = HomogeneousInnerProduct(DirectHomogeneousQuadrature(QR))

        # Time many runs.
        start_time = time()
        for run in range(n_runs):
            IP.quadrature(HAWP)
        end_time = time()

        # Print number of coefficients and average time in ms.
        print N, 1000 * (end_time - start_time) / n_runs

def runMultiD():
    N = 10
    order = 8

    print
    print "Homogeneous quadrature of order", order, "with", N, \
          "coefficients per dimension."
    print "number of dimensions vs. evaluation time [ms]:"

    for D, n_runs in zip([1, 2, 3, 4], [5000, 500, 10, 1]):
        K = HyperCubicShape(D * [N])
        HAWP = HagedornWavepacket(D, 1, 0.2)
        HAWP.set_basis_shapes([K])
        HAWP.set_coefficients([ones((N**D, 1))])

        TQR = TensorProductQR(D * [GaussHermiteQR(order)])
        IP = HomogeneousInnerProduct(DirectHomogeneousQuadrature(TQR))

        start_time = time()
        for run in range(n_runs):
            IP.quadrature(HAWP)
        end_time = time()

        # Print number of dimensions and average time in ms.
        print D, 1000 * (end_time - start_time) / n_runs

def runMultiComponent():
    D = 2
    N = 10
    order = 8
    K = HyperCubicShape(D * [N])
    n_runs = 10

    print
    print "Multi-component", str(D) + "D homogeneous quadrature of order", \
          order, "with", N, "coefficients per dimension."
    print "number of components vs. evaluation time [ms] averaged over", \
          n_runs, "runs:"

    for n_components in range(1, 9):
        HAWP = HagedornWavepacket(D, n_components, 0.2)
        HAWP.set_basis_shapes(n_components * [K])
        HAWP.set_coefficients(n_components * [ones((N**D, 1))])

        TQR = TensorProductQR(D * [GaussHermiteQR(order)])
        IP = InhomogeneousInnerProduct(DirectInhomogeneousQuadrature(TQR))

        start_time = time()
        for run in range(n_runs):
            IP.quadrature(HAWP)
        end_time = time()

        # Print number of dimensions and average time in ms.
        print n_components, 1000 * (end_time - start_time) / n_runs

def runGenzKeister():
    N = 10
    level = 6

    print
    print "Homogeneous Genz-Keister quadrature of level", level, "with", N, \
          "coefficients per dimension."
    print "number of dimensions vs. evaluation time [ms]:"

    for D, n_runs in zip([1, 2, 3, 4], [5000, 500, 10, 1]):
        K = HyperCubicShape(D * [N])
        HAWP = HagedornWavepacket(D, 1, 0.2)
        HAWP.set_basis_shapes([K])
        HAWP.set_coefficients([ones((N**D, 1))])

        QR = GenzKeisterQR(D, level)
        IP = HomogeneousInnerProduct(DirectHomogeneousQuadrature(QR))

        start_time = time()
        for run in range(n_runs):
            IP.quadrature(HAWP)
        end_time = time()

        # Print number of dimensions and average time in ms.
        print D, 1000 * (end_time - start_time) / n_runs

run1D()
runMultiD()
runMultiComponent()
runGenzKeister()
