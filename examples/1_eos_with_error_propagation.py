
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

# - This notebook shows how to propagate uncertainties to get reasonable error bars for pressure.
# - We use MgO as an example.

# # 1. General setup

# In[4]:

import numpy as np
import matplotlib.pyplot as plt
import uncertainties as uct
from uncertainties import unumpy as unp
import pandas as pd
import pytheos as eos


# # 2. Assign uncertainties to the EOS parameters

# Uncertainties of parameters can be defined using the `uncertainties` package.  The parameter values used in this example are just for demonstration purpose.  For more accurate values, please refer to the recent literature.

# In[5]:

v0 = uct.ufloat(74.698, 0.004)
k0 = uct.ufloat(160., 3.)
k0p = uct.ufloat(4.0, 0.3)


# We make an `numpy` array for volume at high pressure.

# In[6]:

n_pts = 20 
vv0 = np.linspace(1.,0.8, n_pts)
v = vv0 * v0


# Calculate pressure from `pytheos`.

# In[7]:

p = eos.bm3_p(v, v0, k0, k0p)


# ### How to get help

# You may get help by using `help(function_name)`.

# In[8]:

help(eos.bm3_p)


# Now you can see that error bars for the EOS parameters are used in error propagation calculation to pressure value.  Note that the uncertainties in the EOS parameters are correctly applied for propagating uncertainties to both molar volume and pressure.

# In[9]:

df = pd.DataFrame()
df['unit-cell volume'] = v
df['pressure'] = p
print(df.to_string(index=False))


# Unfortunately to plot with `matplotlib`, you need to separate nominal values from standard deviation.

# In[10]:

f = plt.figure()
plt.errorbar(unp.nominal_values(p), unp.nominal_values(v), fmt='ko',              xerr=unp.std_devs(p), yerr=unp.std_devs(v))
plt.xlabel('Pressure (GPa)'); plt.ylabel('Unit-cell volume ($\mathrm{\AA}^3$)');


# # 3. Calculate volume from pressure using pytheos

# `Pytheos` provides functions to calculate volume at a given pressure with error propagation.

# In[11]:

v_cal = eos.bm3_v(p, v0, k0, k0p)


# In[12]:

df = pd.DataFrame()
df['pressure'] = p
df['unit-cell volume'] = v_cal
print(df.to_string(index=False))


# Compare this table with the one we showed above for accuracy check.

# # 4. Bulk modulus at high pressure

# You can also propagate uncertainties in bulk modulus calculation.

# In[13]:

k = eos.bm3_k(p, v0, k0, k0p)


# In[14]:

df = pd.DataFrame()
df['pressure'] = p
df['bulk modulus'] = k
print(df.to_string(index=False))


# In[15]:

f = plt.figure()
plt.errorbar( unp.nominal_values(p), unp.nominal_values(k),              xerr=unp.std_devs(p), yerr=unp.std_devs(k), fmt='o')
plt.xlabel('Pressure (GPa)'); plt.ylabel('Bulk modulus (GPa)');


# # 5. High temperature equation of state

# The constant q assumption has been used widely for the thermal part of the mantle phases.  Below, we assign uncertainties to the thermal parameters of MgO.

# In[16]:

gamma0 = uct.ufloat(1.45, 0.02)
q = uct.ufloat(0.8, 0.3)
theta0 = uct.ufloat(800., 0.)


# We will use `constq_pth` for calculating thermal pressure part of the EOS.  Below, I demonstrate how to get help for the function.

# In[17]:

help(eos.constq_pth)


# We calculate total pressure at 2000 K below.  `eos.constq_pth` requires input of volume and temperature with the same number of elements.  For 2000-K isotherm, we generate a temperature array with 2000 for all elements.

# In[18]:

p_hT = eos.bm3_p(v, v0, k0, k0p) + eos.constq_pth(v, np.ones_like(v)*2000., v0, gamma0, q, theta0, 2, 4)


# In[19]:

df = pd.DataFrame()
df['unit-cell volume'] = v
df['pressure@300K'] = p
df['pressure@2000K'] = p_hT
print(df.to_string(index=False))


# In[ ]:



