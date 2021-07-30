from sub.ler_dados import *
from sub.output import *
from sub.solve import *
import pandas as pd

#teste


# def FPNR():
#     #Leitura de dados da rede e das cargas (subrotina ler_dados)

Lines, Z, nbus, nlin = lerdados('NetData14Bus.xlsx','LoadData14Bus.xlsx')


solver(Lines, Z, nbus, nlin)

# if __name__ == "FPNR":
#     FPNR()
print("andreiteste")
print("andreiteste2 qualquer coisa")