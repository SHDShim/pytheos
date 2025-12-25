#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")


# In[3]:


get_ipython().run_line_magic('cat', '0Source_Citation.txt')


# # 0. General note

# - Even for the same data, Birch-Murnaghan equation and Vinet equation yield different $K_0$ and $K_0'$ values. 
# 
# - This notebook shows how to convert $K_0$ and $K'_0$ from the Vinet equation to the Birch-Murnaghan equation.

# # 1. General setup

# In[4]:


import numpy as np
import matplotlib.pyplot as plt
import pytheos as eos
import uncertainties as uct
from uncertainties import unumpy as unp


# # 2. Vinet EOS

# * source: Ye 2017, JGR
# 
# * Pt: $V_0 = 3.9231^3$ A$^{3}$, $K_0 = 277.3$ GPa, $K_0' = 5.226\pm 0.033$.  
# 
# * Au: $V_0 = 4.07860^3$ A$^{3}$, $K_0 = 167.0$ GPa, $K_0' = 5.813\pm 0.022$.  
# 
# * MgO: $V_0 = 74.698$ A$^{3}$, $K_0 = 160.3$ GPa, $K_0' = 4.109\pm 0.022$.  

# In[5]:


v0 = {'Pt': 3.9231**3, 'Au': 4.07860**3, 'MgO': 74.698}
k0 = {'Pt': 277.3, 'Au': 167.0, 'MgO': 160.3}
k0p = {'Pt': uct.ufloat(5.226, 0.033), 'Au': uct.ufloat(5.813, 0.022), 'MgO': uct.ufloat(4.109, 0.022)}


# Set pressure range and number of data points.

# In[6]:


p_max = 150.
n_pts = 100
p = np.linspace(0.,p_max, n_pts)
standard = 'MgO'


# Calculate volume of MgO at different pressures.

# In[7]:


v = eos.vinet_v(p, v0[standard], k0[standard], k0p[standard])


# In[8]:


plt.plot(p,unp.nominal_values(v))
plt.xlabel('Pressure (GPa)')
plt.ylabel('Unit-cell volume ($\mathrm{\AA}^3$)');


# # 3. What if you use those numbers for a wrong equation

# In[9]:


v_bm3 = eos.bm3_v(p, v0[standard], k0[standard], k0p[standard])


# In[10]:


plt.plot(p,unp.nominal_values(v - v_bm3))
plt.xlabel('Pressure (GPa)')
plt.ylabel('$\Delta$ P (GPa)');


# You will systematically underestimate pressure by using $K_0$ and $K_0'$ obtained from Vinet equation in BM3 equation.  

# # 4. Can we reduce the problem?

# Fit the synthetic data to get $K_0'$ for BM equation.  First, setup bm3 model.

# In[11]:


model_bm3 = eos.BM3Model()


# Generate parameters.  Note that `k0p` is uncertainty formatted and therefore, you need to take nominal value only by `.n`. 

# In[12]:


params = model_bm3.make_params(v0=v0[standard], k0=k0[standard], k0p=k0p[standard].n)


# Fix parameters.

# In[13]:


params['v0'].vary = False
params['k0'].vary = False


# Run fitting

# In[14]:


fitresult_bm3 = model_bm3.fit(p, params, v=unp.nominal_values(v))
print(fitresult_bm3.fit_report())


# In[15]:


eos.plot.static_fit_result(fitresult_bm3)


# As you can see above, there is a fundamental difference between BM3 and Vinet.  Therefore, even if you convert $K_0$ and $K_0'$ between BM3 and Vinet, it will create differences for $\pm 0.5$ GPa up to 150 GPa, for example.
