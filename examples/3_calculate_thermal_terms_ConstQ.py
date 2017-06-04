
# coding: utf-8

# In[1]:

get_ipython().magic('cat 0Source_Citation.txt')


# In[2]:

get_ipython().magic('matplotlib inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[3]:

get_ipython().magic("config InlineBackend.figure_format = 'retina'")


# # 0. General note

# This notebook shows how to calculate thermal pressure and associated terms in the constant q approach.

# In[4]:

import numpy as np
import matplotlib.pyplot as plt
import uncertainties as uct
from uncertainties import unumpy as unp
import pandas as pd
import pytheos as eos


# # 1. Calculate Debye energy with uncertainties

# Assign uncertainties to `x`.  

# In[5]:

x = unp.uarray(np.linspace(0.01,15.,20), np.ones(20)*0.5) # 0.1,7.25


# In[6]:

energy = eos.debye_E(x)


# In[7]:

plt.plot(unp.nominal_values(x), unp.nominal_values(energy))
plt.xlabel('x'); plt.ylabel('Energy')
plt.errorbar(unp.nominal_values(x), unp.nominal_values(energy), 
             xerr = unp.std_devs(x), yerr = unp.std_devs(energy));


# # 2. Calculate Gruneisen parameter

# You may get some help on how to call the function using `help()` command.  `constq_grun` calculates Gruneisen $(\gamma)$ parameter with error propagation based on the following relation:
# 
# $$\dfrac{\gamma}{\gamma_0} = \left( \dfrac{V}{V_0} \right)^q$$
# 
# where $\gamma_0$ is the Gruneisen parameter at reference conditions and $V$ is the volume.  $q$ is the logarithmic volume dependence of Gruneisen parameter.

# In[8]:

help(eos.constq_grun)


# Calculate Gruneisen parameter without error bar.

# In[9]:

v0 = 162.3
v = np.linspace(v0, v0*0.8, 20)
grun = eos.constq_grun(v, v0, 1.5, 2)


# In[10]:

plt.plot(v, grun)
plt.xlabel('Unit-cell volume ($\mathrm{\AA}^3$)'); plt.ylabel('$\gamma$');


# The cell below shows how to do error propagation.

# In[11]:

s_v = np.random.uniform(0., 0.1, 20)
v_u = unp.uarray(v, s_v)
gamma = eos.constq_grun(v_u, uct.ufloat(v0, 0.01), uct.ufloat(1.5, 0.1), uct.ufloat(2.,0.5))
gamma


# If you need a pretty table.

# In[12]:

df = pd.DataFrame()
df['volume'] = v_u
df['gamma'] = gamma
print(df.to_string(index=False))


# In[13]:

plt.errorbar(unp.nominal_values(v_u), unp.nominal_values(gamma), xerr=unp.std_devs(v_u), yerr=unp.std_devs(gamma))
plt.xlabel('Unit-cell volume ($\mathrm{\AA}^3$)'); plt.ylabel('$\gamma$');


# You do not need to provide uncertainties for all the parameters.  The cell below shows a case where we do not have error bars for the parameters.  In this case, we have uncertainties for volume.

# In[14]:

eos.constq_grun(v_u, v0, 1.5, 2.)


# # 3. Calculate Debye temperature and thermal pressure

# You can get the Debye temperatures with error bars.

# In[15]:

help(eos.constq_debyetemp)


# In[16]:

eos.constq_debyetemp(v_u, v0, 1.5, 2., 1000.)


# You can get thermal pressures with error bars.

# In[17]:

help(eos.constq_pth)


# In[18]:

p_th = eos.constq_pth(v_u, unp.uarray(np.ones_like(v)*2000., np.ones_like(v)*100), v0, 1.5, 2., 1000., 5, 4)
p_th


# In[19]:

plt.errorbar(unp.nominal_values(v_u), unp.nominal_values(p_th),
            xerr=unp.std_devs(v_u), yerr=unp.std_devs(p_th))
plt.xlabel('Unit-cell volume ($\mathrm{\AA}^3$)'); plt.ylabel('Thermal pressure (GPa)');

