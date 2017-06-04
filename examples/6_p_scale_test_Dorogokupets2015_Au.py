
# coding: utf-8

# ## Source and citation
# 
# - This notebook is a part of the `pytheos` package ([Github](http://github.com/SHDShim/pytheos)). 
# 
# - __[Citation]__ S.-H. Shim (2017) Pytheos - python equations of state tools. doi:

# In[1]:

get_ipython().magic('matplotlib inline')
# %matplotlib notebook # for interactive


# For high dpi displays.

# In[2]:

get_ipython().magic("config InlineBackend.figure_format = 'retina'")


# # 0. General note

# This example compares pressure calculated from `pytheos` and original publication for the gold scale by Dorogokupets 2015.

# # 1. Global setup

# In[3]:

import matplotlib.pyplot as plt
import numpy as np
from uncertainties import unumpy as unp
import pytheos as eos


# # 3. Compare

# In[4]:

eta = np.linspace(1., 0.65, 8)
print(eta)


# In[5]:

dorogokupets2015_au = eos.gold.Dorogokupets2015()


# In[6]:

help(dorogokupets2015_au)


# In[7]:

dorogokupets2015_au.print_equations()


# In[8]:

dorogokupets2015_au.print_equations()


# In[9]:

dorogokupets2015_au.print_parameters()


# In[10]:

v0 = 67.84742110765599


# In[11]:

dorogokupets2015_au.three_r


# In[12]:

v = v0 * (eta) 
temp = 2500.


# In[13]:

p = dorogokupets2015_au.cal_p(v, temp * np.ones_like(v))


# In[14]:

print('for T = ', temp)
for eta_i, p_i in zip(eta, p):
    print("{0: .3f} {1: .2f}".format(eta_i, p_i))


# Table is not given for this publication.

# In[16]:

v = dorogokupets2015_au.cal_v(p, temp * np.ones_like(p), min_strain=0.6)
print((v/v0))

