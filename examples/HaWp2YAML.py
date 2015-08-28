import yaml
import csv
import numpy

from WaveBlocksND import *
from HaWp2CSV import *

def convertComplexNumpyArrayToYAML(array):
    return array.astype(str).tolist()

def ParametersetAsYAML(paramset):
    (q,p,Q,P) = paramset
    return {
        "q": q.astype(float).T[0].tolist(), 
        "p": p.astype(float).T[0].tolist(),
        "Q": convertComplexNumpyArrayToYAML(Q),
        "P": convertComplexNumpyArrayToYAML(P)
    }

def BasisshapeAsYAML(shape):
    if type(shape) is HyperCubicShape:
        return {
            "type": "hypercubic",
            "dimensionality": shape.get_dimension(),
            "limits": shape.get_limits()
        }
    elif type(shape) is HyperbolicCutShape:
        return {
            "type": "hyperbolic",
            "dimensionality": shape.get_dimension(),
            "sparsity": shape._sparsity
        }
    elif type(shape) is LimitedHyperbolicCutShape:
        return {
            "type": "hyperbolic",
            "dimensionality": shape.get_dimension(),
            "sparsity": shape._sparsity,
            "limits": shape.get_limits()
        }
    else:
        raise Exception("I dont know how to serialize " + str(shape));

def WriteHaWpAsYAML(packet, name):
    dom_params = []
    if type(packet) is HagedornWavepacketInhomogeneous:
        for item in packet.get_parameters(key=('q', 'p', 'Q', 'P')):
            dom_params.append(ParametersetAsYAML(item))
    else:
        dom_params.append(ParametersetAsYAML(packet.get_parameters(key=('q', 'p', 'Q', 'P'))))
    
    dom_shapes = []
    for item in packet.get_basis_shapes():
        dom_shapes.append(BasisshapeAsYAML(item))
    
    '''
    dom_coeffs = []
    for index, item in enumerate(packet.get_coefficients()):
        filename = name+"_coefficients_"+str(index)+".csv"
        WriteHaWpCoefficientsAsCSV(item, list(packet.get_basis_shapes(index).get_node_iterator()), filename)
        dom_coeffs.append(filename)
    '''
    
    dom = {
        "eps": packet.get_eps(),
        "dimensionality": packet.get_dimension(),
        "n_components": packet.get_number_components(),
        #"coefficients": dom_coeffs,
        "parameters": dom_params,
        "shapes": dom_shapes
    }
    
    yaml.safe_dump(dom, file(name, "w"))