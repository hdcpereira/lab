#!/usr/bin/env python
# coding: utf-8
#%%
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import h,k,c
from matplotlib.offsetbox import AnchoredText
#%%
#defining constants that will be used in the analysis.
tensao_res = 0.621
tensao_lamp = 8.01
res_resistor = 1
temp0 = 21.1+273
wien = 2.898e-3
#%%
#defining planck function that will be used to fit the data.
def planck_func(x,norm,t,fund,planck_const,boltzmann_const):
    return norm*((8*np.pi*planck_const*c)/np.power(x,5))*(1/(np.exp((planck_const*c)/(x*boltzmann_const*t))-1))+fund
#%%
#filtering data so then we are able to proceed to curve fitting.

corrente = tensao_res/res_resistor
res_lampada = tensao_lamp /corrente

temp = (temp0 * np.power((res_lampada/res_resistor),(1/1.24)))+273

data = pd.read_csv("dataset/CN_Eq1_m6_2.txt",skiprows = 2,header = None, delimiter = '\t', names=['voltasMP','intensidade'])
data = data.apply(lambda x: x.str.replace(',','.'))

data['voltasMP'] = pd.to_numeric(data['voltasMP'], errors = 'coerce')
data['intensidade'] = pd.to_numeric(data['intensidade'], errors = 'coerce')

mask = data['voltasMP']>400

x = data['voltasMP'][mask]
y = data['intensidade'][mask]

x_not_filtered = data['voltasMP']
y_not_filtered = data['intensidade']

x= (x/60)*(np.pi/180)
x= (0.001/300)*np.sin(x)

x_not_filtered = (x_not_filtered/60)*(np.pi/180)
x_not_filtered = (0.001/300)*np.sin(x_not_filtered)

fundo = y.min()

sigma_y = np.sqrt(np.abs(y))

#%%
## ploting raw set of data, without any filter.

plt.plot(x_not_filtered,y_not_filtered, '.')
plt.xlabel('Comprimento de Onda $\lambda$ (m)')
plt.ylabel('Intensidade $R_T(\lambda)$ (u.a.)')
plt.grid()
plt.savefig("raw_dataset1.jpeg",dpi = 300)
plt.clf()

#%%
ans, cov = curve_fit(planck_func, x, y, p0=[0,temp, fundo,h,k])

y_fit = planck_func(x,*ans)

fitted_df = pd.DataFrame({'wavelenght':x, 'R':y, 'fit':y_fit})
peak = fitted_df.loc[fitted_df['fit'] == fitted_df['fit'].max()]
temperatura = np.divide(wien,peak['wavelenght'])

#%%
##plotting only one set of data.
fig, ax1 = plt.subplots()

ax1.errorbar(x,y,yerr = sigma_y , fmt = '.',elinewidth=0.5, color = 'green',alpha = 0.5, label = 'dataset #2')

parameters_text = AnchoredText('norm = {:.4E} \ntemp = {:.2f} \nfundo = {:.4E} \nplanck = {:.4E} \nboltzmann = {:.4E}'.format(ans[0],ans[1],ans[2],ans[3],ans[4]), loc = 'upper right')

ax1.plot(x,y_fit, '-', color ='red', label = 'fit')

text = '{:.2f}K'.format(temperatura.values[0])
# ax1.annotate(text,(peak['wavelenght'],peak['R']))

ax1.legend(loc='center left', bbox_to_anchor=(0.5, 1.05),
          ncol=1, fancybox=True, shadow=True)

ax1.add_artist(parameters_text) 

plt.grid()
plt.xlabel("Comprimento de Onda")
plt.savefig("exp1_dataset2.jpeg", dpi = 300)
plt.show()
# %%
