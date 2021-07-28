import math

def flowres(vmag, vang, Ybus, G, B, Z, line, nbus, nlin):

    IPA=np.zeros(nbus)
    IPR=np.zeros(nbus)
    FPA=np.zeros(nbus)
    FPA1=np.zeros(nbus)
    FPR=np.zeros(nbus)
    FPR1=np.zeros(nbus)
    FC=np.zeros(nbus, np.complex64)
    FC1=np.zeros(nbus, np.complex64)

    Vbus=np.zeros(nbus, np.complex64)
    Ibus=np.zeros(nbus, np.complex64)

    # Calculo das injecoes de potencia ativa e reativa nas barras
    for i in range(nbus):
        for j in range(nbus):
            IPA[i]=IPA[i]+vmag[j]*(G[i,j]*math.cos(vang[i]-vang[j])+B[i,j]*math.sin(vang[i]-vang[j]))
            IPR[i]=IPR[i]+vmag[j]*(G[i,j]*math.sin(vang[i]-vang[j])-B[i,j]*math.cos(vang[i]-vang[j]))
        IPA=vmag[i]*IPA[i]
        IPR=vmag[i]*IPR[i]

    # Calculo das injecoes de corrente nas barras (forma retangular)
    for i in range(nbus):
        Vbus[i]=complex(vmag[i]*math.cos(vang[i]),vmag[i]*math.sin(vang[i]))

    Ibus=Ybus*Vbus


    # Calculo dos fluxos de potencia ativa e reativa nos ramos
    for index, row in line.iterrows():
        From=row['De']-1
        To=row['Para']-1
        zs=complex(row['R'])
        ys=pow(zs,-1)
        gs=ys.real
        bs=ys.imag
        bsh=row['Bsh']
        FPA[i]=pow(vmag[From],2)*gs-vmag[From]*vmag[To]*(gs*math.cos(vang[From]-vang[To])+bs*math.sin(vang[From]-vang[To]))
        FPA1[i]=pow(vmag[To],2)*gs-vmag[From]*vmag[To]*(gs*math.cos(vang[To]-vang[From])+bs*math.sin(vang[To]-vang[From]))
        FPR[i]=-pow(vmag[From],2)*(bsh+bs)-vmag[From]*vmag[To]*(gs*math.sin(vang[From]-vang[To])-bs*math.cos(vang[From]-vang[To]))
        FPR1[i]=pow(vmag[To],2)*(bsh+bs)-vmag[From]*vmag[To]*(gs*math.sin(vang[To]-vang[From])-bs*math.cos(vang[To]-vang[From]))
        
        # Calculo dos fluxos de corrente nos ramos (forma retangular)
        FC[i]=ys*(Vbus[From]-Vbus[To])+complex(0,bsh)*Vbus[From]
        FC1[i]=ys*(Vbus[To]-Vbus[From])+complex(0,bsh)*Vbus[To]



    return FPA, FPA1, FPR, FPR1, IPA, IPR, Ibus, FC, FC1