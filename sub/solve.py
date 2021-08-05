
from sub.deltaz import *
from sub.flowres import *
from sub.jacob import *
import numpy as np

def solver(Lines, Z, nbus, nlin):
    # Determinacao do estado (Metodo de Newton Raphson)
    # Inicializacoes

    zest = []
    tol = 0.0001 # tolerancia para convergencia (a ser confrontada com os mismatches deltaP e deltaQ)
    maxiter = 10 # numero maximo de iteracoes
    indconv = 0 # indicador de convergencia (assume valor 1 quando o processo iterativo converge)
    iter=0;

    Ybarra = np.zeros((nbus, nbus), np.complex64)
    
    #
    # Montagem da matriz Ybarra
    for index, line in Lines.iterrows():
        zs=complex(line['R'],line['X']);
        y= 1/zs
        bsh=complex(0,line['Bsh'])/2;
        de = int(line['De'])-1
        para = int(line['Para'])-1
        tap = np.float64(line['Tap'])
        Ybarra[de][de]=Ybarra[de][de]+(y/pow(tap,2))+bsh  # elemento diagonal (i,i)
        Ybarra[para][para]=Ybarra[para][para]+y+bsh       # elemento diagonal (j,j)
        Ybarra[de][para]= - y/tap                         # elemento fora da diagonal (i,j)
        Ybarra[para][de]= - y/tap                         # elemento fora da diagonal (j,i)

    G = Ybarra.real
    B = Ybarra.imag

    for i in range(nbus):
        B[i,i]=B[i,i] + Z['CS'][i+1] 

    # Inicializacao do estado (flat start)
    # vmag = nbus*[1]
    # vang = nbus*[0]

    vmag = np.ones(nbus)
    vang = np.zeros(nbus)

    for i in range(nbus):
        if(Z['Tipo'][i+1] > 0):
            vmag[i]=Z['V'][i+1]
        else:
            vmag[i]=1

    #Processo Iterativo
    while(iter <= maxiter and indconv !=1):
        
        # Calculo dos residuos - mismatches deltaP e deltaQ (subrotina deltaz)
        zest, res, zloc, nP, nq, neq = delta_z(vmag, vang, Z, G, B, nbus)

        # Verificacao da convergencia
        absres=abs(res)
        maxres=max(res)  # maior mismatch observado
        
        if(maxres <= tol):
            indconv=1
        
        # Formacao da matriz Jacobiano (subrotina jacob)
        J = jacob(vmag, vang, G, B, Z, zloc, nbus, nP, nq, neq)

        # Atualizacao do estado 
        # deltax = np.dot(np.linalg.inv(J),res.transpose())
        deltax = np.inv(J) * res[:, None]

        # Solucao do problema linear para determinar a modificacao a ser realizada no estado atual 
        for i in range(neq):
            izloc = int(zloc[i])
            if(i <= nP):
                vang[izloc]=vang[izloc]+deltax[i] # atualizacao do angulo da tensao
            else:
                vmag[izloc]=vmag[izloc]+deltax[i] # atualizacao da magnitude da tensao
            iter=iter+1; #incrementa contador de iteracoes

        vang*(180/3.1416)

        # Calculo dos fluxos/injecoes de potencia e corrente para o estado convergido (subrotina flowres)
        FPA, FPA1, FPR, FPR1, IPA, IPR, Ibarra, FC, FC1=flowres(vmag, vang, Ybarra, G, B, Z, Lines, nbus, nlin)