import numpy as np

def jacob(vmag, vang, G, B, z, zloc, nbus, nP, nq, neq):

    
    J = np.zeros((neq,neq))
    
    return J