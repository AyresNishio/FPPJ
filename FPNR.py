from sub.deltaz import *
from sub.flowres import *
from sub.jacob import *
from sub.ler_dados import *
from sub.output import *
from sub.solve import *
import pandas as pd

#teste


# def FPNR():
#     #Leitura de dados da rede e das cargas (subrotina ler_dados)

line, z, nbus, nlin = lerdados('NetData14Bus.xlsx','LoadData14Bus.xlsx')

print("hahaha")
print(line)

# if __name__ == "FPNR":
#     FPNR()