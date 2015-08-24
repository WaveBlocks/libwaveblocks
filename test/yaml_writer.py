import yaml
import numpy

from WaveBlocksND import *

def convertComplexNumpyArrayToYAML(array):
    return array.astype(str).tolist()

def yaml_write_hawp(packet, name):
    (q,p,Q,P) = packet.get_parameters(key=("q","p","Q","P"))
    
    dom = {
        "eps": packet.get_eps(),
        "dimensionality": packet.get_dimension(),
        "coefficients": name+"_coefficients.csv",
        "type": "scalar",
        "parameterset": {
            "q": q.astype(float).T[0].tolist(), 
            "p": p.astype(float).T[0].tolist(),
            "Q": convertComplexNumpyArrayToYAML(Q),
            "P": convertComplexNumpyArrayToYAML(P)
        }
    }
    
    yaml.safe_dump(dom, file(name+".yaml", "w"))