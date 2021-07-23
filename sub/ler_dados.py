import pandas as pd

def lerdados(dadoslinha, dadoscarga):
    line = pd.read_excel(dadoslinha) 
    z = pd.read_excel(dadoscarga)
    nlin = len(line)
    nbus = len(z)
    return  line, z, nbus, nlin