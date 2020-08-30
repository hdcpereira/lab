#!/usr/bin/env python
# coding: utf-8
#%%
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import h,k,c
#%%

#calculando a temperatura
tensao_res = 0.621
tensao_lamp = 8.01
res_resistor = 1
temp0 = 21.1+273

corrente = tensao_res/res_resistor
res_lampada = tensao_lamp /corrente

temp = (temp0 * np.power((res_lampada/res_resistor),(1/1.24)))+273
#%%

#lendo os dados
data = pd.read_csv("CN_Eq1_m6_1.txt",skiprows = 2,header = None, delimiter = '\t', names=['voltasMP','intensidade'])
data = data.apply(lambda x: x.str.replace(',','.'))

data['voltasMP'] = pd.to_numeric(data['voltasMP'], errors = 'coerce')
data['intensidade'] = pd.to_numeric(data['intensidade'], errors = 'coerce')

mask = data['voltasMP']>400

x = data['voltasMP'][mask]
y = data['intensidade'][mask]

#%%
#correções nos dados

x= (x/60)*(np.pi/180)
x= (0.001/300)*np.sin(x)

fundo = y.min()
y = y - fundo
#%%
plt.plot(x,y, 'x')
plt.xlabel('comprimento (m)')
# plt.ylabel("Comprimento de Onda ()")

plt.savefig("experimento1.jpeg")

#%%
def planck_func(x,norm,t):
    return norm*((8*np.pi*h*c)/np.power(x,5))*(1/(np.exp((h*c)/(x*k*t))-1))


ans, cov = curve_fit(planck_func, x, y, p0=[4,temp])

y_fit = planck_func(x,*ans)


#%%


fig, ax1 = plt.subplots()

ax1.plot(x,y,'x')

ax2 = ax1.twinx()

ax2.plot(x,y_fit, '--', color ='red')

plt.xlabel("Comprimento de Onda")
fig.tight_layout()

plt.savefig("exp1.jpeg")
plt.show()
# %%
