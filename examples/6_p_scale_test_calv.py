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

# This example compares pressure calculated from `pytheos` and original publication for the MgO scale by Dorogokupets 2015.

# # 1. Global setup

# In[4]:


import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[5]:


eta = np.linspace(1., 0.6, 9)
print(eta)


# In[6]:


test_EOS = eos.periclase.Zha2000()


# In[7]:


help(eos.periclase.Zha2000)


# In[8]:


test_EOS.print_equations()


# In[9]:


test_EOS.print_equations()


# In[10]:


test_EOS.print_parameters()


# In[11]:


v0 = 74.698


# In[12]:


p = np.linspace(0, 10000, 11) 
temp = 3000.


# In[13]:


v = test_EOS.cal_v(p, temp * np.ones_like(p), min_strain=0.1)
v_300 = test_EOS.cal_v(p, 300. * np.ones_like(p))


# Table is not given in this publication.

# In[14]:


print('for T = ', temp)
for eta_i, v_i in zip(p, v):
    print("{0: .3f} {1: .2f} {2: .2f}".format(eta_i, v_i, test_EOS.cal_p(v_i, 300.)))


# In[15]:


import pytheos as eos

p = 100000
v0 = 74.698
k0 = 160.2
k0p = 4.03

#74.698+/-0.001, 'k0': 160.2+/-0, 'k0p': 4.03+/-0

v_spline = eos.bm3_v_single_spline(p, v0, k0, k0p)
v_inversion = eos.bm3_v_single(p, v0, k0, k0p)

print(v_spline, v_inversion, (v_spline - v_inversion)/v_inversion)


# In[ ]:




