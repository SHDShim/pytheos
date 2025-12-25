#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('cat', '0Source_Citation.txt')


# In[2]:


get_ipython().run_line_magic('matplotlib', 'inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[3]:


get_ipython().run_line_magic('config', "InlineBackend.figure_format = 'retina'")


# # 0. General note

# This example compares pressure calculated from `pytheos` and original publication for the NaCl scales.

# # 1. Global setup

# In[4]:


import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 2. Compare

# In[20]:


eta = np.linspace(1., 0.65, 8)
print(eta)


# In[21]:


dorogokupets2007 = eos.sodium_chloride.Dorogokupets2007()


# In[22]:


help(eos.sodium_chloride.Dorogokupets2007)


# In[23]:


dorogokupets2007.print_equations()


# In[24]:


dorogokupets2007.print_parameters()


# In[25]:


v0 = 179.44


# In[26]:


dorogokupets2007.three_r


# In[33]:


v = v0 * (eta) 
temp = 2000.


# In[34]:


p = dorogokupets2007.cal_p(v, temp * np.ones_like(v))


# <img src='./tables/Dorogokupets2007_NaCl.png' width=500>

# In[36]:


print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f} ".format(eta_i, p_i))


# Max difference is about 0.02 GPa

# In[37]:


v = dorogokupets2007.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print((v/v0))

