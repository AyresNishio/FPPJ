import numpy as np

def flowres(vmag, vang, Ybus, G, B, Z, line, nbus, nlin):

    IPA=np.zeros(nbus, np.float64)
    IPR=np.zeros(nbus, np.float64)
    FPA=np.zeros(nlin, np.float64)
    FPA1=np.zeros(nlin, np.float64)
    FPR=np.zeros(nlin, np.float64)
    FPR1=np.zeros(nlin, np.float64)
    FC=np.zeros(nlin, np.complex64)
    FC1=np.zeros(nlin, np.complex64)

    Vbus=np.zeros(nbus, np.complex64)
    Ibus=np.zeros(nbus, np.complex64)

    # Calculo das injecoes de potencia ativa e reativa nas barras
    for i in range(nbus):
        for j in range(nbus):
            IPA[i]=IPA[i]+vmag[j]*(G[i,j]*np.cos(vang[i]-vang[j])+B[i,j]*np.sin(vang[i]-vang[j]))
            IPR[i]=IPR[i]+vmag[j]*(G[i,j]*np.sin(vang[i]-vang[j])-B[i,j]*np.cos(vang[i]-vang[j]))
        IPA[i]=vmag[i]*IPA[i]
        IPR[i]=vmag[i]*IPR[i]

    # Calculo das injecoes de corrente nas barras (forma retangular)
    for i in range(nbus):
        Vbus[i]=complex(vmag[i]*np.cos(vang[i]),vmag[i]*np.sin(vang[i]))

    Ibus=np.dot(Ybus,Vbus)


    # Calculo dos fluxos de potencia ativa e reativa nos ramos
    for index, row in line.iterrows():
        i = index - 1
        From=int(row['De']-1)
        To=int(row['Para']-1)
        zs=complex(row['R'], row['X'])
        ys=pow(zs,-1)
        gs=ys.real
        bs=ys.imag
        bsh=row['Bsh']/2
        FPA[i]=pow(vmag[From],2)*gs-vmag[From]*vmag[To]*(gs*np.cos(vang[From]-vang[To])+bs*np.sin(vang[From]-vang[To]))
        FPA1[i]=pow(vmag[To],2)*gs-vmag[From]*vmag[To]*(gs*np.cos(vang[To]-vang[From])+bs*np.sin(vang[To]-vang[From]))
        FPR[i]=-pow(vmag[From],2)*(bsh+bs)-vmag[From]*vmag[To]*(gs*np.sin(vang[From]-vang[To])-bs*np.cos(vang[From]-vang[To]))
        FPR1[i]=pow(vmag[To],2)*(bsh+bs)-vmag[From]*vmag[To]*(gs*np.sin(vang[To]-vang[From])-bs*np.cos(vang[To]-vang[From]))
        
        # Calculo dos fluxos de corrente nos ramos (forma retangular)
        FC[i]=ys*(Vbus[From]-Vbus[To])+complex(0,bsh)*Vbus[From]
        FC1[i]=ys*(Vbus[To]-Vbus[From])+complex(0,bsh)*Vbus[To]



    return FPA, FPA1, FPR, FPR1, IPA, IPR, Ibus, FC, FC1