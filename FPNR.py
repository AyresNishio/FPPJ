from sub.ler_dados import *
from sub.output import *
from sub.deltaz import *
from sub.solve import *
import pandas as pd

#teste


# def FPNR(network_file,load_file):
#     #Leitura de dados da rede e das cargas (subrotina ler_dados)

Lines, Z, nbus, nlin = lerdados('NetData14Bus.xlsx','LoadData14Bus.xlsx')
    # Lines, Z, nbus, nlin = lerdados(network_file,load_file)
Vmag,Vang,FPA_dp, FPA_pd, FPR_dp, FPR_pd, IPA, IPR, Ibarra, FC, FC1 = solver(Lines, Z, nbus, nlin)
tipo = ['V']*len(Vmag) + ['P']*len(FPA_dp) + ['P']*len(FPA_pd) + ['P']*len(IPA) + ['Q']*len(FPR_dp) + ['Q']*len(FPR_pd) + ['Q']*len(IPR)
De = list(range(1,nbus+1)) + list(Lines['De']) + list(Lines['Para']) + list(range(1,nbus+1)) + list(Lines['De']) + list(Lines['Para']) + list(range(1,nbus+1))
Para = ['-']*len(Vmag) + list(Lines['Para']) + list(Lines['De']) + ['-']*len(Vmag) + list(Lines['Para']) + list(Lines['De']) + ['-']*len(Vmag) 
Valor = Vmag.tolist() + FPA_dp.tolist() + FPA_pd.tolist() + FPR_dp.tolist() + FPR_pd.tolist() + IPA.tolist() + IPR.tolist()

# Desvios baseados em:
                    # IEEE TRANS. INSTRUMENTATION AND MEASUREMENT
                    # , VOL. 66, NO. 8, pp. 2036-2045, AUG. 2017,
                    #  “A Cubature Kalman Filter Based Power System Dynamic State Estimator”,
                    #  A. Sharma, S. C. Srivastava, and S. Chakrabarti
desvio_padrao = []
for i in tipo:
    if i == 'V':
        desvio_padrao.append(0.006)
    if i == 'P': 
        desvio_padrao.append(0.01)
    if i == 'Q': 
        desvio_padrao.append(0.01)
    #if i == 'Ang':
    #    desvio_padrao.append(0.018)

#Montagem do DataFrame        
valores = {'Tipo':tipo,'De':De,'Para': Para,'Valor': Valor,'Desvio Padrão':desvio_padrao}

df = pd.DataFrame(valores)
print(df)
# if __name__ == "FPNR":
#     FPNR()
