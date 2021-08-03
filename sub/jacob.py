import numpy as np

def jacob(vmag, vang, G, B, z, zloc, nbus, nP, nq, neq):

    J = np.zeros((neq,2*nbus))
    J1 = np.zeros((neq,neq)) #True Jacobian
    ivar = -1
    # Montagem matriz Jacobiano (expressoes livro Abur)
    # Varre equacao por equacao e verifica em quais colunas serao geradas derivadas (em funcao da barra ser PQ ou PV)  
    # Para cada equacao sao calculadas duas derivadas (para barras PQ) ou uma derivada (para barras PV)  
    for i in range(nP+1):
        iloc = int(zloc[i])
        for j in range(nbus):
            j1 = j + nbus
            if (z["Tipo"][j+1] == 0):
                if(j == iloc):
                    J[i,j] = 0
                    J[i,j1] = 0
                    for k in range(nbus): 
                        J[i,j] = J[i,j] \
                        +vmag[j]*vmag[k]*(-G[j,k]*np.sin(vang[j]-vang[k])+B[j,k]*np.cos(vang[j]-vang[k])) #derivada dPi/dtetai (diag - Hii)
                        
                        J[i,j1] = J[i,j1]+vmag[k]*(G[j,k]*np.cos(vang[j]-vang[k])+B[j,k]*np.sin(vang[j]-vang[k])) #derivada dPi/dVi (diag - Nii))
          
                    J[i,j] = J[i,j]-(vmag[j]**2)*B[j,j]
                    J[i,j1] = J[i,j1]+vmag[j]*G[j,j]
                else:
                    J[i,j] = vmag[iloc]*vmag[j]*(G[iloc,j]*np.sin(vang[iloc]-vang[j])-B[iloc,j]*np.cos(vang[iloc]-vang[j])); # derivada dPi/dtetaj (fora diag - Hij)
                    J[i,j1] = vmag[iloc]*(G[iloc,j]*np.cos(vang[iloc]-vang[j])+B[iloc,j]*np.sin(vang[iloc]-vang[j])); # derivada dPi/dVj (fora diag - Nij)
       
            elif (z["Tipo"][j+1] == 1): # barra j eh PV - apenas tetaj eh variavel do problema
        
                if (j == iloc):
                    J[i,j]=0
                    J[i,j1]=0
                    for k in range(nbus): 
                        J[i,j] = J[i,j]+vmag[j]*vmag[k]*(-G[j,k]*np.sin(vang[j]-vang[k])+B[j,k]*np.cos(vang[j]-vang[k])) # derivada dPi/dtetai (diag - Hii)
                    J[i,j] = J[i,j]-(vmag[j]**2)*B[j,j]
                else:
                    J[i,j] = vmag[iloc]*vmag[j]*(G[iloc,j]*np.sin(vang[iloc]-vang[j])-B[iloc,j]*np.cos(vang[iloc]-vang[j])) # derivada dPi/dtetaj (fora diag - Hij) 
           
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------                       
    for i in range(nP+1,neq): # derivadas da injecao de potencia reativa ## conferir neq ##
        iloc = int(zloc[i])
        for j in range(nbus):
            j1 = j + nbus
            if(z["Tipo"][j+1] == 0): # barra j eh PQ - tetaj e Vj sao variaveis do problema        
                if(j == iloc):
                    J[i,j]=0
                    J[i,j1]=0
                    for k in range(nbus): 
                        J[i,j] = J[i,j]+vmag[j]*vmag[k]*(G[j,k]*np.cos(vang[j]-vang[k])+B[j,k]*np.sin(vang[j]-vang[k])) # derivada dQi/dtetai (diag - Mii)
                        J[i,j1] = J[i,j1]+vmag[k]*(G[j,k]*np.sin(vang[j]-vang[k])-B[j,k]*np.sin(vang[j]-vang[k])) # derivada dQi/dVi (diag - Lii)

                    J[i,j] = J[i,j]-(vmag[j]**2)*G[j,j]
                    J[i,j1] = J[i,j1]-vmag[j]*B[j,j]
                else:
                    J[i,j] = vmag[iloc]*vmag[j]*(-G[iloc,j]*np.cos(vang[iloc]-vang[j])-B[iloc,j]*np.sin(vang[iloc]-vang[j])) # derivada dQi/dtetaj (fora diag - Mij)
                    J[i,j1] = vmag[iloc]*(G[iloc,j]*np.sin(vang[iloc]-vang[j])-B[iloc,j]*np.cos(vang[iloc]-vang[j])) # derivada dQi/dVj (fora diag - Lij)
                
            elif(z["Tipo"][j+1] == 1): # barra j eh PV - apenas tetaj eh variavel do problema        
     
                if (j == iloc):
                    J[i,j]=0
                    J[i,j1]=0
                    for k in range(nbus): 
                        J[i,j] = J[i,j]+vmag[j]*vmag[k]*(G[j,k]*np.cos(vang[j]-vang[k])+B[j,k]*np.sin(vang[j]-vang[k])) # derivada dQi/dtetai (diag - Mii)
    
                    J[i,j] = J[i,j]-(vmag[j]**2)*G[j,j]
                else:
                    J[i,j] = vmag[iloc]*vmag[j]*(-G[iloc,j]*np.cos(vang[iloc]-vang[j])-B[iloc,j]*np.sin(vang[iloc]-vang[j])) # derivada dQi/dtetaj (fora diag - Mij)

    
    for j in range(2*nbus):

        if (j <= nbus - 1):
            j1 = j
            if (z["Tipo"][j1 + 1] == 1):
                ivar = ivar + 1
                J1[:, ivar] = J[:, j]
        else:
            j1=j-nbus
        if (z["Tipo"][j1 + 1] == 0):
            ivar = ivar + 1
            J1[:, ivar] = J[:, j]

    J = J1

    return J