import math

def flowres(vmag, vang, Ybarra, G, B, Z, line, nbus, nlin):

    IPA=np.zeros(nbus)
    IPR=np.zeros(nbus)
    ICA=np.zeros(nbus)
    ICR=np.zeros(nbus)
    FPA=np.zeros(nbus)
    FPR=np.zeros(nbus)
    FPA1=np.zeros(nbus)
    FPR1=np.zeros(nbus)
    FCA1=np.zeros(nbus)
    FCR1=np.zeros(nbus)

    Ibarra=np.zeros(nbus)
    FC=np.zeros(nbus) 
    FC1=np.zeros(nbus)

    for i in range(nbus):
        for j in range(nbus):
            IPA[i]=IPA[i]+vmag[j]*(G[i,j]* math.cos())

    return FPA, FPA1, FPR, FPR1, IPA, IPR, Ibarra, FC, FC1