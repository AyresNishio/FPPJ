import numpy as np

def delta_z(vmag, vang, Z, G, B, nbus):

    zest=np.zeros(nbus, np.float64)
    zestp=np.zeros(nbus, np.float64)
    zestq=np.zeros(nbus, np.float64)

    zloc=np.zeros(nbus, np.float64)
    zlocp=np.zeros(nbus, np.float64)
    zlocq=np.zeros(nbus, np.float64)
    
    res=np.zeros(nbus, np.float64)
    resp=np.zeros(nbus, np.float64)
    resq=np.zeros(nbus, np.float64)
    
    nP=-1
    nQ=-1
    
    # Calculo das injecoes de potencia e dos residuos deltaP e deltaQ
    for i in range(nbus):

        if Z['Tipo'][i+1] == 0: # Barras PQ
            nP=nP+1
            nQ=nQ+1

            # injecao de potencia ativa
            for j in range(nbus):
                zestp[nP]=zestp[nP]+vmag[j]*(G[i,j]*np.cos(vang[i]-vang[j])+B[i,j]*np.sin(vang[i]-vang[j]))
            
            zestp[nP]=vmag[i]*zestp[nP]
            resp[nP]=(Z['Pg'][i+1]-Z['Pl'][i+1])-zestp[nP] # vetor contendo os residuos de potencia ativa (mismatches deltaP = Pesp - Pcal)
            zlocp[nP]=i

            # injecao de potencia reativa
            for j in range(nbus):
                zestq[nQ]=zestq[nQ]+vmag[j]*(G[i,j]*np.sin(vang[i]-vang[j])-B[i,j]*np.cos(vang[i]-vang[j]))
            
            zestq[nQ]=vmag[i]*zestq[nQ]
            #resq[nQ]=(-Z.at[i,'Ql'])-zestq[nQ] # vetor contendo os residuos de potencia reativa (mismatches deltaQ = Qesp - Qcal)
            resq[nQ]=(-Z['Ql'][i+1])-zestq[nQ] # vetor contendo os residuos de potencia reativa (mismatches deltaQ = Qesp - Qcal)
            zlocq[nQ]=i
        
        elif Z['Tipo'][i+1] == 1: # Barras PV
            nP=nP+1

            # injecao de potencia ativa
            for j in range(nbus):
                zestp[nP]=zestp[nP]+vmag[j]*(G[i,j]*np.cos(vang[i]-vang[j])+B[i,j]*np.sin(vang[i]-vang[j]))
            
            zestp[nP]=vmag[i]*zestp[nP]
            #resp[nP]=(Z.at[i,'Pg']-Z.at[i,'Pl'])-zestp[nP]
            resp[nP]=(Z['Pg'][i+1]-Z['Pl'][i+1])-zestp[nP]
            zlocp[nP]=i

    nEQ=(nP+1)+(nQ+1)
    zest=np.concatenate((zestp,zestq),axis=0) # numero total de equacoes (deltaP e deltaQ)
    res=np.concatenate((resp[0:nP+1],resq[0:nQ+1]),axis=0) # res - vetor que agrega os subvetores resp e resq
    zloc=np.concatenate((zlocp[0:nP+1],zlocq[0:nQ+1]),axis=0)


    return zest, res, zloc, nP, nQ, nEQ