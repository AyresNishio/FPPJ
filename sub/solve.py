import numpy as np
import cmath 

def solver(Lines, z, nbus, nlin):
    # Determinacao do estado (Metodo de Newton Raphson)
    # Inicializacoes
    vang = []
    vmag = []
    zest = []
    tol = 0.0001 # tolerancia para convergencia (a ser confrontada com os mismatches deltaP e deltaQ)
    maxiter = 10 # numero maximo de iteracoes
    indconv = 0 # indicador de convergencia (assume valor 1 quando o processo iterativo converge)
    iter=0;

    array =  np.arange(nbus*nbus).reshape(nbus,nbus)
    Ybarra = np.zeros_like(array,np.complex64)
    
    # Montagem da matriz Ybarra
    for line in Lines.iterrows():
        zs=complex(line['R'],line['X']);
        y= 1/zs
        bsh=complex(0,line['Bsh'])/2;
        de = int(line['De'])-1
        para = int(line['Para'])-1
        tap = np.float64(line['Tap'])
        Ybarra[de][de]=Ybarra[de][de]+(y/pow(tap,2))+bsh  # elemento diagonal (i,i)
        Ybarra[para][para]=Ybarra[para][para]+y+bsh       # elemento diagonal (j,j)
        Ybarra[de][para]= - y/tap                         # elemento fora da diagonal (i,j)
        Ybarra[para][de]= - y/tap                       # elemento fora da diagonal (j,i)