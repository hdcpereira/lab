#%%
import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.constants import h,k,c

#%%
#calculando a temperatura
tensao_res = [0.621, 0.675, 0.705, 0.747, 0.786]
tensao_lamp = [8.01, 9, 9.89, 11.12, 12.31]
res_resistor = 1
temp0 = 21.1+273

corrente = np.divide(tensao_res,res_resistor) #tensao_res/res_resistor
res_lampada = np.divide(tensao_lamp,corrente) #tensao_lamp /corrente

temp = (temp0 * np.power((res_lampada/res_resistor),(1/1.24)))+273


#%%
# defining plank function

def planck_func(x,norm,t,fund):
    return norm*((8*np.pi*h*c)/np.power(x,5))*(1/(np.exp((h*c)/(x*k*t))-1))-fund

#%%
#lendo os dados
data = [None]*5
mask = [None]*5
x = [None]*5
y = [None]*5
fundo = [None]*5
ans = [None]*5
cov = [None]*5
y_fit = [None]*5

ini_guess = [4,40000,10000,10000,11999]

for i in range(1,6):
    filename = "CN_Eq1_m6_{}.txt".format(i)
    data[i-1] = pd.read_csv(filename,skiprows = 2,header = None, delimiter = '\t', names=['voltasMP','intensidade'])

    data[i-1] = data[i-1].apply(lambda x: x.str.replace(',','.'))

    data[i-1]['voltasMP'] = pd.to_numeric(data[i-1]['voltasMP'], errors = 'coerce')
    data[i-1]['intensidade'] = pd.to_numeric(data[i-1]['intensidade'], errors = 'coerce')

    mask[i-1] = data[i-1]['voltasMP']>400

    x[i-1] = data[i-1]['voltasMP'][mask[i-1]]
    y[i-1] = data[i-1]['intensidade'][mask[i-1]]
#correção nos dados 
    x[i-1]= (x[i-1]/60)*(np.pi/180)
    x[i-1]= (0.001/300)*np.sin(x[i-1])

    fundo[i-1] = y[i-1].min()
    # y[i-1] = y[i-1]

    ans[i-1], cov[i-1] = curve_fit(planck_func, x[i-1], y[i-1], p0=[ini_guess[i-1],temp[i-1],fundo[i-1]])

    y_fit[i-1] = planck_func(x[i-1],*ans[i-1])


#%%
fig, ax1 = plt.subplots()


ax1.plot(x[0],y[0],'x')
ax1.plot(x[1],y[1],'x')
ax1.plot(x[2],y[2],'x')
ax1.plot(x[3],y[3],'x')
ax1.plot(x[4],y[4],'x')

ax1.plot(x[0],y_fit[0], '--', color ='red')
ax1.plot(x[1],y_fit[1], '--', color ='red')
ax1.plot(x[2],y_fit[2], '--', color ='red')
ax1.plot(x[3],y_fit[3], '--', color ='red')
ax1.plot(x[4],y_fit[4], '--', color ='red')

plt.xlabel("Comprimento de Onda")
plt.ylabel("R")

fig.tight_layout()

plt.savefig("exp1.jpeg")
plt.show()

# %%
